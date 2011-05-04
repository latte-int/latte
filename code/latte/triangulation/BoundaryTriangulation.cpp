/* BoundaryTriangulation.cpp -- Boundary triangulation
	       
   Copyright 2006, 2007 Matthias Koeppe
   Copyright 2011       Christof Soeger

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

#include "./simplify_helpers.cpp"
#include <list>
#include <algorithm>

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

  //cerr << "Completed rays: " << ray_matrix << "Inverse: " << inverse << endl;
  
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
      //cerr << "Cone " << k << " full-dimensional: ";
      // full-dimensional cone
      int l;
      for (l = 0; l<numOfVars - 1; l++) {
	//cerr << f[k * (numOfVars - 1) + l] << " ";
	if (f[k * (numOfVars - 1) + l] == 0) {
	  //cerr << "Not singularity-avoiding" << endl;
	  return false;
	}
      }
      //cerr << endl;
    }
  }
  return true;
}

static void 
prepare_interior_vector_data(listCone *boundary_triangulation, int numOfVars,
  mat_ZZ &alpha, mat_ZZ &U, mat_ZZ &UF)
{
  /* Create a matrix whose rows generate a d-dimensional lattice. */
  int num_cones = lengthListCone(boundary_triangulation);
  alpha.SetDims(numOfVars, num_cones);
  mat_ZZ F; //need to avoid zeros here (for standart basis)
  F.SetDims(numOfVars, num_cones * (numOfVars - 1));
  listCone *cone;
  int k;
  for (k = 0, cone = boundary_triangulation; cone!=NULL; k++, cone=cone->rest) {
    alpha_basis(cone, numOfVars, alpha, k, F, k * (numOfVars - 1));
  }
  cerr << "Auxiliary lattice for constructing an interior vector: " << endl
       << alpha << endl;
  cerr << "Need to avoid zeros here: " << endl
       << F << endl;
  /* Find a short vector in this lattice. */
  ZZ det2;
  U.SetDims(numOfVars, numOfVars);
  cerr << "Computing LLL basis: " << endl;
  long rank = LLL(det2, alpha, U, 1, 1);
  cerr << "Rank: " << rank << " Variables: " << numOfVars << endl;
  assert(rank == numOfVars);

  /* avoid zeros here (in L3 basis) */
  UF = U * F;
}

static vec_ZZ 
construct_interior_vector(listCone *boundary_triangulation, int numOfVars, vec_ZZ &det_vector)
{
  int num_cones = lengthListCone(boundary_triangulation);
  mat_ZZ alpha, U, UF;
  prepare_interior_vector_data(boundary_triangulation, numOfVars, alpha, U, UF);
  cerr << "L3 basis: " << endl
       << U << endl;
  cerr << "In auxiliary space:" << endl
       << alpha << endl;
  cerr << "In zero-avoidance space:" << endl
       << UF;
  
  /* Determine if the LLL basis vectors avoid the zeros. */
  /* We simply choose the first vector of the reduced basis. */
  /* Check basis vectors. */
  int i;
  for (i = 0; i<numOfVars; i++) {
    cerr << "# Basis vector " << i << ":" << endl;
    if (singularity_avoiding(num_cones, numOfVars, alpha[i], UF[i])) {
      cerr << "### Found vector" << endl;
    }
  }
  /* Try integer combinations of the basis vectors */
  int *max_comb = new int[numOfVars];
  int max_coeff;
  for (max_coeff = 1; ; max_coeff++) {
    cerr << "# Trying combinations with maximum coefficient " << max_coeff << endl;
    for (i = 0; i<numOfVars; i++)
      max_comb[i] = 2 * max_coeff;
    IntCombEnum iter_comb(max_comb, numOfVars);
    int *next;
    vec_ZZ multi;
    multi.SetLength(numOfVars);
    while((next = iter_comb.getNext())) {
      int j;
      bool have_max = false;
      for (j = 0; j<numOfVars; j++) {
	multi[j] = next[j] - max_coeff;
	if (abs(multi[j]) == max_coeff) have_max = true; 
      }
      if (have_max) {
	vec_ZZ alpha_comb;
	alpha_comb = multi * alpha;
	vec_ZZ F_comb;
	F_comb = multi * UF;
	if (singularity_avoiding(num_cones, numOfVars, alpha_comb, F_comb)) {
	  cerr << "### Found vector:" << endl
	       << " multipliers: " << multi << endl
	       << " vector:      " << multi * U << endl
	       << " alpha_comb:  " << alpha_comb << endl
	       << " F_comb:      " << F_comb << endl;
	  det_vector = alpha_comb;
	  return multi * U;
	}
      }
    }
  }
  cerr << "No suitable vector found." << endl;
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
    cerr << "Facet " << i+1 << "/" << num_inequalities << ": ";
    int j;
    for (j = 0; j<incidence->setsize; j++)
      if (set_member(j + 1, incidence->set[i])) {
	cerr << j + 1 << " ";
      }
    cerr << "(cardinality " << set_card(incidence->set[i]) << ")" << endl;
    /* Compute a triangulation of that facet. */
    listCone *facet_cone
      = cone_from_ray_set(rays, incidence->set[i], cone->vertex);
//     triangulate_cone_with_cddlib(facet_cone, Parameters,
// 				 delone_height, NULL,
// 				 Parameters->Number_of_Variables - 1,
// 				 consumer);
    triangulate_cone_with_4ti2(facet_cone, Parameters,
			       random_height, &Parameters->triangulation_max_height,
			       consumer);
//     cerr << "Triangulation of facet cone: " << lengthListCone(facet_triangulation)
// 	 << " simplicial cones." << endl;
  }
}


void
generate_interior_vector_simplified_data(listCone *boundary_triangulation, int numOfVars,
  mat_ZZ &alpha, mat_ZZ &basis, mat_ZZ &F)
{
  prepare_interior_vector_data(boundary_triangulation, numOfVars, alpha, basis, F);
  cerr << "Simplify data" << endl;
  // transpose and work on the rows
  alpha = transpose(alpha); //the determinants
  F = transpose(F);         //the last components of the facet normal vectors
  list<vec_ZZ> a_list, F_list;
  for (int i=0; i<alpha.NumRows(); i++) {
    a_list.push_back(alpha[i]);
  }
  for (int i=0; i<F.NumRows(); i++) {
    F_list.push_back(F[i]);
  }
  //make the first row entry positive by row*(-1) if needed
  for_each(a_list.begin(), a_list.end(), make_first_entry_positive);
  for_each(F_list.begin(), F_list.end(), make_first_entry_positive);
  //make the rows in F prime
  for_each(F_list.begin(), F_list.end(), make_coprime);
  //sort
  a_list.sort(vec_ZZ_is_less);
  F_list.sort(vec_ZZ_is_less);
  //unique
  a_list.unique();
  F_list.unique();
  //remove zero vector if existent (it has to be on the first position)
  if(IsZero(a_list.front())!=0) a_list.pop_front();
  //remove rows of alpha that have multiples
  list<vec_ZZ>::iterator it = a_list.begin();
  list<vec_ZZ>::const_iterator cit;
  int fnz; //first non zero position
  ZZ factor; 
  while (it!=a_list.end()) {
    fnz = first_non_zero_pos(*it);
    cit = it; ++cit;
    while (cit!=a_list.end() && fnz == first_non_zero_pos(*cit)) {
      factor = (*cit)[fnz]/(*it)[fnz];
      if ( factor>1 && (*cit)[fnz]%(*it)[fnz] == 0 && //first position is multiple
           (*it)*factor == (*cit) ) {
        fnz = -1; //to mark success
        it = a_list.erase(it);
		break;
      }
	  ++cit;
	}
    if (fnz >= 0) ++it;
  }
  //transpose back
  alpha.SetDims(a_list.size(), numOfVars);
  F.SetDims(F_list.size(), numOfVars);
  for (int i=0; i<alpha.NumRows(); i++) {
    alpha[i] = a_list.front();
	a_list.pop_front();
  }
  for (int i=0; i<F.NumRows(); i++) {
    F[i] = F_list.front();
	F_list.pop_front();
  }
  alpha = transpose(alpha);
  F = transpose(F);
}


void
complete_boundary_triangulation_of_cone_with_subspace_avoiding_facets
(listCone *boundary_triangulation, BarvinokParameters *Parameters, ConeConsumer &consumer)
{
  int numOfVars = Parameters->Number_of_Variables;
  /* Complete the cones with an interior ray. */
  vec_ZZ det_vector;
  vec_ZZ extra_ray_vector
    = construct_interior_vector(boundary_triangulation, Parameters->Number_of_Variables, det_vector);
  cerr << "Interior ray vector: " << extra_ray_vector << endl;
  cerr << "Resulting determinants: " << det_vector << endl;
  complete_boundary_triangulation_of_cone_with_subspace_avoiding_facets
    (boundary_triangulation, Parameters, extra_ray_vector, consumer);
}

ZZ determinant(listVector* rows, int numOfVars) {
  mat_ZZ row_matrix;
  row_matrix.SetDims(numOfVars, numOfVars);
  for (int i = 0; rows!=NULL; i++, rows = rows->rest) {
    row_matrix[i] = rows->first;
  }
  return determinant(row_matrix);
}


void
complete_boundary_triangulation_of_cone_with_subspace_avoiding_facets
(listCone *boundary_triangulation, BarvinokParameters *Parameters, vec_ZZ &extra_ray_vector, ConeConsumer &consumer)
{
  int numOfVars = Parameters->Number_of_Variables;
  /* Complete the cones with an interior ray. */
  ZZ det, scalar;
  cerr << "Extra ray vector: " << extra_ray_vector << endl;
  vec_ZZ interior_vector;
  interior_vector.SetLength(numOfVars);
  vec_ZZ seperating_hyperplane;
  seperating_hyperplane.SetLength(numOfVars);
  listCone *resulting_triangulation = NULL;

  listCone *simplicial_cone, *next_simplicial_cone;
  int i;
  listVector* rays;
  for (simplicial_cone = boundary_triangulation, i = 0;
   simplicial_cone!=NULL;
   simplicial_cone = next_simplicial_cone, i++) {
    next_simplicial_cone = simplicial_cone->rest;
    
    rays = simplicial_cone->rays;
    for (; rays!=NULL; rays = rays->rest) {
      interior_vector += rays->first;
    }

    simplicial_cone->rays
      = appendVectorToListVector(extra_ray_vector, simplicial_cone->rays);
    det = determinant(simplicial_cone->rays, numOfVars);
    if (det != 0) {
      /* A full-dimensional cone is created. */
      simplicial_cone->rest = resulting_triangulation;
      resulting_triangulation = simplicial_cone;
    }
    else {
      /* A lower-dimensional cone is created -- discard it. */
      freeCone(simplicial_cone);
    }
  }

  dualizeCones(resulting_triangulation, Parameters->Number_of_Variables, Parameters);
  listCone *cone;
  for (cone = resulting_triangulation; cone!=NULL; cone=cone->rest) {
    if (!rays_ok(cone, Parameters->Number_of_Variables)) {
      cerr << "Note: The following dualized-back cone has bad rays." << endl;
      printCone(cone, numOfVars);
    }
  }
  dualizeCones(resulting_triangulation, Parameters->Number_of_Variables, Parameters);
  cerr << "Interior vector: " << interior_vector << endl;

  {
    listCone *next = NULL;
    for (cone = resulting_triangulation; cone!=NULL; cone=next) {
      next = cone->rest;
      cone->rest = NULL;
		
      //update Coefficient
      InnerProduct(scalar, interior_vector, cone->facets->first);
      //scalar always non-zero because interior_vector is not in the boundary
      assert(scalar != 0);
      //cerr<<scalar<<" ";
      cone->coefficient = (scalar>0)?-1:1;
      //add up to get a seperating hyperplane
      seperating_hyperplane += cone->coefficient * cone->facets->first;

      consumer.ConsumeCone(cone);
    }
    cerr << "Seperating hyperplane: " << seperating_hyperplane << endl;
    InnerProduct(scalar, extra_ray_vector, seperating_hyperplane);
    if (scalar>0) {
      cerr << "Resulting Cone may be not pointed, change the sign of the extra vector!" << endl;
      exit(1);
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
