/* BoundaryTriangulation.cpp -- Boundary triangulation
	       
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
#include "triangulation/BoundaryTriangulation.h"
#include "triangulation/RegularTriangulationWithCddlib.h"
#include "latte_cddlib.h"
#include "latte_random.h"
#include "print.h"
#include "dual.h"

using namespace std;

static bool
rays_ok(listCone *cone, int numOfVars)
{
  listVector *ray;
  for (ray = cone->rays; ray!=NULL; ray=ray->rest) {
    if (ray->first[numOfVars - 1] == 0)
      return false;
  }
  return true;
}

static void
check_rays(listCone *cone, int numOfVars)
{
  if (!rays_ok(cone, numOfVars)) {
    cerr << "The following dualized-back cone has bad rays." << endl;
    printCone(cone, numOfVars);
    abort();
  }
}

/* Construct a vector whose entries are the numbers alpha_k(w) for w =
   e_1, ..., e_d.
*/
static vec_ZZ
alpha_basis(listCone *cone, int numOfVars)
{
  int num_rays = lengthListVector(cone->rays);
  assert(num_rays == numOfVars - 1);
  /* We first complete the rays to a basis of the space. */
  mat_ZZ ray_matrix;
  ray_matrix.SetDims(numOfVars, numOfVars);
  int i;
  listVector *ray;
  for (i = 0, ray = cone->rays; ray!=NULL; i++, ray = ray->rest) {
    ray_matrix[i] = ray->first;
  }
  /* Simple and stupid method: Try all unit vectors. */
  int j;
  for (j = 0; j<numOfVars; j++) {
    ray_matrix[i][j] = 1;
    if (determinant(ray_matrix) != 0) break;
    ray_matrix[i][j] = 0;
  }
  assert(j<numOfVars);
  mat_ZZ inverse;
  ZZ d;
  inv(d, inverse, ray_matrix);

  //cout << "Completed rays: " << ray_matrix << "Inverse: " << inverse << endl;
  
  vec_ZZ result;
  result.SetLength(numOfVars);
  for (j = 0; j<numOfVars; j++)
    result[j] = inverse[j][numOfVars - 1];
  return result;
}

static vec_ZZ 
construct_interior_vector(listCone *boundary_triangulation, int numOfVars, vec_ZZ &det_vector)
{
  /* Create a matrix whose rows generate a d-dimensional lattice. */
  int num_cones = lengthListCone(boundary_triangulation);
  mat_ZZ alpha;
  alpha.SetDims(numOfVars, num_cones);
  listCone *cone;
  int k;
  for (k = 0, cone = boundary_triangulation; cone!=NULL; k++, cone=cone->rest) {
    vec_ZZ alpha_basis_k = alpha_basis(cone, numOfVars);
    int j;
    for (j = 0; j<numOfVars; j++) {
      alpha[j][k] = alpha_basis_k[j];
    }
  }
  cout << "Auxiliary lattice for constructing an interior vector: " << endl
       << alpha << endl;
  /* Find a short vector in this lattice. */
  ZZ det2;
  mat_ZZ U;
  U.SetDims(numOfVars, numOfVars);
  long rank = LLL(det2, alpha, U, 1, 1);
  cout << "Rank: " << rank << " Variables: " << numOfVars << endl;
  assert(rank == numOfVars);
  /* We simply choose the first vector of the reduced basis. */
  det_vector = alpha[0];
  return U[0];
}  

void
boundary_triangulation_of_cone_with_subspace_avoiding_facets
(listCone *cone, BarvinokParameters *Parameters, ConeConsumer &consumer)
{
  int numOfVars = Parameters->Number_of_Variables;
  listCone *boundary_triangulation = NULL;
  vector<listVector *> rays = ray_array(cone);
  dd_PolyhedraPtr poly
    = cone_to_cddlib_polyhedron(cone, Parameters->Number_of_Variables);
  dd_MatrixPtr inequalities = dd_CopyInequalities(poly);
  assert(inequalities->representation == dd_Inequality);
  int num_inequalities = inequalities->rowsize;
  /* For each computed facet, obtain the set of input rays that are
     incident with the facet. */
  dd_SetFamilyPtr incidence = dd_CopyIncidence(poly);
  assert(incidence->setsize == lengthListVector(cone->rays));
  assert(incidence->famsize == num_inequalities);
  int i;
  for (i = 0; i<num_inequalities; i++) {
    cout << "Facet " << i << ": ";
    int j;
    for (j = 0; j<incidence->setsize; j++)
      if (set_member(j + 1, incidence->set[i])) {
	cout << j + 1 << " ";
      }
    cout << "(cardinality " << set_card(incidence->set[i]) << ")" << endl;
    /* Compute a triangulation of that facet. */
    listCone *facet_cone
      = cone_from_ray_set(rays, incidence->set[i], cone->vertex);
    CollectingConeConsumer ccc;
    triangulate_cone_with_cddlib(facet_cone, Parameters,
				 delone_height, NULL,
				 Parameters->Number_of_Variables - 1,
				 ccc);
    listCone *facet_triangulation = ccc.Collected_Cones;
    cout << "Triangulation of facet cone: " << lengthListCone(facet_triangulation)
	 << " simplicial cones." << endl;
    boundary_triangulation
      = appendListCones(facet_triangulation, boundary_triangulation);
  }
  /* Complete the cones with an interior ray. */
  vec_ZZ det_vector;
  vec_ZZ interior_ray_vector
    = construct_interior_vector(boundary_triangulation, Parameters->Number_of_Variables, det_vector);
  cout << "Interior ray vector: " << interior_ray_vector << endl;
  cout << "Reuslting determinants: " << det_vector << endl;
  listCone *resulting_triangulation = NULL;

  listCone *simplicial_cone, *next_simplicial_cone;
  for (simplicial_cone = boundary_triangulation, i = 0;
       simplicial_cone!=NULL;
       simplicial_cone = next_simplicial_cone, i++) {
    next_simplicial_cone = simplicial_cone->rest;
    if (det_vector[i] != 0) {
      /* A full-dimensional cone is created. */
      simplicial_cone->rays
	= appendVectorToListVector(interior_ray_vector, simplicial_cone->rays);
      simplicial_cone->rest = resulting_triangulation;
      resulting_triangulation = simplicial_cone;
    }
    else {
      /* A lower-dimensional cone is created -- discard it. */
      freeCone(simplicial_cone);
    }
  }

  resulting_triangulation
    = dualizeBackCones(resulting_triangulation, Parameters->Number_of_Variables);
  for (cone = resulting_triangulation; cone!=NULL; cone=cone->rest) {
    if (!rays_ok(cone, Parameters->Number_of_Variables)) {
      cerr << "Note: The following dualized-back cone has bad rays." << endl;
      printCone(cone, numOfVars);
    }
  }
  resulting_triangulation
    = dualizeBackCones(resulting_triangulation, Parameters->Number_of_Variables);

  for (cone = resulting_triangulation; cone!=NULL; cone=cone->rest)
    consumer.ConsumeCone(cone);
}
