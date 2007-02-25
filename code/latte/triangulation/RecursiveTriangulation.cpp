/* RecursiveTriangulation.cpp -- Recursive triangulation of a cone given by rays
	       
   Copyright 2006 Matthias Koeppe

   This file is part of LattE.
   
   LattE is free software; you can redistribute it and/or modify it
   under the terms of the version 2 of the GNU General Public License
   as published by the Free Software Foundation.

   LattE is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with LattE; if not, write to the Free Software Foundation,
   Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
*/

#include <cassert>
#include <vector>
#include "triangulation/RecursiveTriangulation.h"
#include "latte_random.h"
#include "print.h"
#include "dual.h"
#include "triangulation/RegularTriangulationWithCdd.h"
#include "Irrational.h"

using namespace std;

static void
check_rays(listCone *cone, int numOfVars)
{
  listVector *ray;
  for (ray = cone->rays; ray!=NULL; ray=ray->rest) {
    if (ray->first[numOfVars - 1] == 0) {
      cerr << "The following dualized-back cone has bad rays." << endl;
      printCone(cone, numOfVars);
      abort();
    }
  }
}

/* Consume CONE. */
static listCone *
triangulate_recursively
(listCone *cone, BarvinokParameters *Parameters, int level = 0)
{
#if 1
  dualizeCone(cone, Parameters->Number_of_Variables, Parameters);
  //  printCone(cone, Parameters->Number_of_Variables);
  check_rays(cone, Parameters->Number_of_Variables);
  dualizeCone(cone, Parameters->Number_of_Variables, Parameters);
#endif
  int numOfVars = Parameters->Number_of_Variables;
  int num_rays = lengthListVector(cone->rays);
  if (num_rays <= numOfVars) {
    /* FIXME: Check full-dimensional */
    assert(num_rays == numOfVars);
    return cone;
  }
  /* We pick (numOfVars - 1) random rays that define
     a hyperplane with a permissible normal vector that
     partitions the cone into left and right subcone.
  */
#if 1
  cerr << endl << "Level " << level << ": "
       << "Cone with " << num_rays << " rays:  " << flush;
#endif
  vector<listVector*> partition_rays(numOfVars - 1);
  int num_tries;
  const int max_tries = 1000;
  for (num_tries = 0; num_tries<max_tries; num_tries++) {
    vector<bool> taken(num_rays);
    int i;
    for (i = 0; i<num_rays; i++) taken[i] = false;
    /* Pick numOfVars - 1 random rays. */
    int num_remaining_rays = num_rays;
    for (i = 0; i<numOfVars - 1; i++) {
      int index = uniform_random_number(0, num_remaining_rays - 1);
      int j;     /* Counts rays */
      int k = 0; /* Counts non-taken rays */
      listVector *ray;
      for (j = 0, ray = cone->rays; k!=index; ray = ray->rest, j++) {
	if (!taken[j]) k++;
      }
      taken[j] = true;
      partition_rays[i] = ray;
      num_remaining_rays--;
#if 0
      cerr << j << " ";
#endif
    }
    /* Compute the normal vector of partition_rays. */
    mat_ZZ ray_matrix;
    vec_ZZ rhs;
    vec_ZZ normal;
    ray_matrix.SetDims(numOfVars, numOfVars);
    rhs.SetLength(numOfVars);
    normal.SetLength(numOfVars);
    for (i = 0; i<numOfVars - 1; i++) {
      int j;
      for (j = 0; j<numOfVars; j++) 
	ray_matrix[j][i] = partition_rays[i]->first[j];
      rhs[i] = 0;
    }
    /* Unit vector. */
    ray_matrix[numOfVars - 1][numOfVars - 1] = 1;
    rhs[numOfVars - 1] = 1;

    ZZ d;
    solve(d, normal, ray_matrix, rhs);
    if (d == 0) {
#if 0
      cerr << "Picked linearly dependent rays, trying again." << endl;
#endif
      continue;
    }
    cerr << "Normal vector: " << normal << " ";
    /* Check subspace condition. */
    if (normal[numOfVars - 1] == 0) {
#if 0
      cerr << "Normal vector lies in the forbidden subspace, trying again."
	   << endl;
#endif
      continue;
    }
    /* Check partition. */
    int num_left = 0;
    int num_right = 0;
    listVector *ray;
    cerr << "Orientation: ";
    for (i = 0, ray = cone->rays; i<num_rays; i++, ray = ray->rest) {
      ZZ ip;
      InnerProduct(ip, normal, ray->first);
      if (ip > 0) {
	num_right++;
	cerr << i << "(+) ";
      }
      else if (ip < 0) {
	num_left++;
	cerr << i << "(-) ";
      }
      else {
	cerr << i << "(0) ";
      }
    }
    if (num_left == 0 || num_right == 0) {
      cerr << "Picked rays that do not give a partition, trying again." << endl;
      continue;
    }
    /* Do partition. */
    listVector *left_rays = NULL;
    listVector *right_rays = NULL;
    for (i = 0, ray = cone->rays; i<num_rays; i++, ray = ray->rest) {
      ZZ ip;
      InnerProduct(ip, normal, ray->first);
      if (ip >= 0)
	right_rays = appendVectorToListVector(ray->first, right_rays);
      if (ip <= 0)
	left_rays = appendVectorToListVector(ray->first, left_rays);
    }
    listCone *left_cone = createListCone();
    listCone *right_cone = createListCone();
    left_cone->vertex = new Vertex(*cone->vertex);
    left_cone->rays = left_rays;
    right_cone->vertex = new Vertex(*cone->vertex);
    right_cone->rays = right_rays;
    freeCone(cone);
    /* Recurse. */
    listCone *left_triang
      = triangulate_recursively(left_cone, Parameters, level + 1);
    listCone *right_triang
      = triangulate_recursively(right_cone, Parameters, level + 1);
    return appendListCones(left_triang, right_triang);
  };
  cerr << "No success after " << num_tries << " tries, dualizing back and irrationalizing." << endl;
  cerr << "The cone: " << endl;
  printCone(cone, numOfVars);
  dualizeCone(cone, numOfVars, Parameters);
  /* Check that rays do not lie in the forbidden subspace. */
  check_rays(cone, numOfVars);
  irrationalizeCone(cone, numOfVars);
  listCone *triang = triangulate_cone_with_cdd(cone, Parameters);
  freeCone(cone);
  triang = dualizeBackCones(triang, numOfVars);
  return triang;
}

/* Does not consume CONE. */
listCone *
triangulate_cone_recursively_with_subspace_avoiding_facets
(listCone *cone, BarvinokParameters *Parameters)
{
  return triangulate_recursively(copyCone(cone), Parameters);
}
