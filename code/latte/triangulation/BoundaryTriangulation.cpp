/* BoundaryTriangulation.cpp -- Boundary triangulation
	       
   Copyright 2006, 2007 Matthias Koeppe

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
#include "triangulation/RegularTriangulationWith4ti2.h"
#include "latte_cddlib.h"
#include "latte_random.h"
#include "print.h"
#include "dual.h"
#include "genFunction/IntCombEnum.h"

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

   Also construct a matrix F whose columns collect the last coordinates
   of all facet normals created, when the extra vector w runs through
   e_1, ..., e_d. 
*/
static void
alpha_basis(listCone *cone, int numOfVars, 
	    mat_ZZ &alpha, int alpha_column_index,
	    mat_ZZ &F, int F_start_column_index)
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
  
  for (j = 0; j<numOfVars; j++)
    alpha[j][alpha_column_index] = inverse[j][numOfVars - 1];
  
  int s;
  for (s = 0; s<numOfVars; s++) { // Consider unit vector w = e_s.
    int l;
    for (l = 0; l<numOfVars - 1; l++) {
      // Compute the (d+1)-st coordinate of $\tilde b^*_{l,d+1}(e_s)$
      
      F[s][l + F_start_column_index]
	= inverse[s][numOfVars - 1] * inverse[numOfVars - 1][l]
	- inverse[s][l] * inverse[numOfVars - 1][numOfVars - 1];
      
    }
  }
}

static bool
singularity_avoiding(int num_cones, int numOfVars,
		     const vec_ZZ &alpha, const vec_ZZ &f)
{
  int k;
  for (k = 0; k<num_cones; k++) {
    if (alpha[k] != 0) {
      //cout << "Cone " << k << " full-dimensional: ";
      // full-dimensional cone
      int l;
      for (l = 0; l<numOfVars - 1; l++) {
	//cout << f[k * (numOfVars - 1) + l] << " ";
	if (f[k * (numOfVars - 1) + l] == 0) {
	  //cout << "Not singularity-avoiding" << endl;
	  return false;
	}
      }
      //cout << endl;
    }
  }
  return true;
}

static vec_ZZ 
construct_interior_vector(listCone *boundary_triangulation, int numOfVars, vec_ZZ &det_vector)
{
  /* Create a matrix whose rows generate a d-dimensional lattice. */
  int num_cones = lengthListCone(boundary_triangulation);
  mat_ZZ alpha;
  alpha.SetDims(numOfVars, num_cones);
  mat_ZZ F;
  F.SetDims(numOfVars, num_cones * (numOfVars - 1));
  listCone *cone;
  int k;
  for (k = 0, cone = boundary_triangulation; cone!=NULL; k++, cone=cone->rest) {
    alpha_basis(cone, numOfVars, alpha, k, F, k * (numOfVars - 1));
  }
  cout << "Auxiliary lattice for constructing an interior vector: " << endl
       << alpha << endl;
  cout << "Need to avoid zeros here: " << endl
       << F << endl;
  /* Find a short vector in this lattice. */
  ZZ det2;
  mat_ZZ U;
  U.SetDims(numOfVars, numOfVars);
  long rank = LLL(det2, alpha, U, 1, 1);
  cout << "Rank: " << rank << " Variables: " << numOfVars << endl;
  assert(rank == numOfVars);

  /* Determine if the LLL basis vectors avoid the zeros. */
  mat_ZZ UF;
  UF = U * F;
  cout << "L3 basis: " << endl
       << U << endl;
  cout << "In auxiliary space:" << endl
       << alpha << endl;
  cout << "In zero-avoidance space:" << endl
       << UF;
  
  /* We simply choose the first vector of the reduced basis. */
  /* Check basis vectors. */
  int i;
  for (i = 0; i<numOfVars; i++) {
    cout << "# Basis vector " << i << ":" << endl;
    if (singularity_avoiding(num_cones, numOfVars, alpha[i], UF[i])) {
      cout << "### Found vector" << endl;
    }
  }
  /* Try integer combinations of the basis vectors */
  int max_coeff = 5;
  int *max_comb = new int[numOfVars];
  for (i = 0; i<numOfVars; i++)
    max_comb[i] = 2 * max_coeff;
  IntCombEnum iter_comb(max_comb, numOfVars);
  int *next;
  vec_ZZ multi;
  multi.SetLength(numOfVars);
  while((next = iter_comb.getNext())) {
    int j;
    for (j = 0; j<numOfVars; j++)
      multi[j] = next[j] - max_coeff;
    vec_ZZ alpha_comb;
    alpha_comb = multi * alpha;
    vec_ZZ F_comb;
    F_comb = multi * UF;
    if (singularity_avoiding(num_cones, numOfVars, alpha_comb, F_comb)) {
      cout << "### Found vector:" << endl
	   << multi * U << endl
	   << alpha_comb << endl
	   << F_comb << endl;
      det_vector = alpha_comb;
      return multi * U;
    }
  }  
  cout << "No suitable vector found." << endl;
  exit(1);
}  

void
compute_triangulation_of_boundary
(listCone *cone, BarvinokParameters *Parameters, ConeConsumer &consumer)
{
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
    cout << "Facet " << i+1 << "/" << num_inequalities << ": ";
    int j;
    for (j = 0; j<incidence->setsize; j++)
      if (set_member(j + 1, incidence->set[i])) {
	cout << j + 1 << " ";
      }
    cout << "(cardinality " << set_card(incidence->set[i]) << ")" << endl;
    /* Compute a triangulation of that facet. */
    listCone *facet_cone
      = cone_from_ray_set(rays, incidence->set[i], cone->vertex);
//     triangulate_cone_with_cddlib(facet_cone, Parameters,
// 				 delone_height, NULL,
// 				 Parameters->Number_of_Variables - 1,
// 				 consumer);
    triangulate_cone_with_4ti2(facet_cone, Parameters,
			       random_height, &Parameters->triangulation_max_height,
			       Parameters->Number_of_Variables - 1,
			       consumer);
//     cout << "Triangulation of facet cone: " << lengthListCone(facet_triangulation)
// 	 << " simplicial cones." << endl;
  }
}

void
complete_boundary_triangulation_of_cone_with_subspace_avoiding_facets
(listCone *boundary_triangulation, BarvinokParameters *Parameters, ConeConsumer &consumer)
{
  int numOfVars = Parameters->Number_of_Variables;
  /* Complete the cones with an interior ray. */
  vec_ZZ det_vector;
  vec_ZZ interior_ray_vector
    = construct_interior_vector(boundary_triangulation, Parameters->Number_of_Variables, det_vector);
  cout << "Interior ray vector: " << interior_ray_vector << endl;
  cout << "Reuslting determinants: " << det_vector << endl;
  listCone *resulting_triangulation = NULL;

  listCone *simplicial_cone, *next_simplicial_cone;
  int i;
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
  listCone *cone;
  for (cone = resulting_triangulation; cone!=NULL; cone=cone->rest) {
    if (!rays_ok(cone, Parameters->Number_of_Variables)) {
      cerr << "Note: The following dualized-back cone has bad rays." << endl;
      printCone(cone, numOfVars);
    }
  }
  resulting_triangulation
    = dualizeBackCones(resulting_triangulation, Parameters->Number_of_Variables);

  {
    listCone *next = NULL;
    for (cone = resulting_triangulation; cone!=NULL; cone=next) {
      next = cone->rest;
      cone->rest = NULL;
      consumer.ConsumeCone(cone);
    }
  }
}

void
boundary_triangulation_of_cone_with_subspace_avoiding_facets
(listCone *cone, BarvinokParameters *Parameters, ConeConsumer &consumer)
{
  CollectingConeConsumer ccc;
  compute_triangulation_of_boundary(cone, Parameters, ccc);
  complete_boundary_triangulation_of_cone_with_subspace_avoiding_facets
    (ccc.Collected_Cones, Parameters, consumer);
}
