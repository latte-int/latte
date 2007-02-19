/* RegularTriangulationWithCddlib.cpp -- Regular triangulation using CDDLIB
	       
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
#include "latte_cddlib.h"
#include "latte_gmp.h"
#include "latte_random.h"
#include "RegularTriangulationWithCddlib.h"

using namespace std;

listCone *
cone_from_ray_set(vector<listVector *> &rays,
		  set_type ray_set,
		  Vertex *vertex)
{
  listCone *c = createListCone();
  c->vertex = new Vertex(*vertex);
  vector<listVector *>::iterator i;
  int j;
  for (i = rays.begin(), j = 1; i!=rays.end(); ++i, j++) {
    if (set_member(j, ray_set)) {
      c->rays = new listVector((*i)->first, c->rays);
    }
  }
  return c;
}

listCone *
triangulate_cone_with_cddlib(listCone *cone,
			     BarvinokParameters *Parameters,
			     height_function_type height_function,
			     void *height_function_data,
			     int cone_dimension)
{
  listCone *triangulation = NULL;
  // Copy rays into an array, so we can index them.
  vector<listVector *> rays = ray_array(cone);
  /* Create a matrix from the rays, with 2 extra coordinates in
     front:  The 0th for homogenization, the 1st for lifting. */
  dd_MatrixPtr matrix
    = rays_to_cddlib_matrix(cone->rays, Parameters->Number_of_Variables,
			    /* num_homogenization_vars: */ 2,
			    /* num_extra_rows: */ 1);
  int num_rays = lengthListVector(cone->rays);
  assert(num_rays + 1 == matrix->rowsize);
  /* Extra row: `vertical' ray -- This kills all upper facets.
     See Verdoolaege, Woods, Bruynooghe, Cools (2005). */
  dd_set_si(matrix->matrix[num_rays][1], 1);
  do {
    /* Compute random lifting. */
    int i;
    for (i = 0; i<num_rays; i++) {
      height_function(matrix->matrix[i][1], rays[i]->first, height_function_data);
    }
    /* Output of the file -- for debugging. */
    {
      FILE *mf = fopen("lifted_cone_for_triangulation", "w");
      dd_WriteMatrix(mf, matrix);
      fclose(mf);
    }
    /* Compute facets by double-description method. */
    dd_ErrorType error;
    dd_PolyhedraPtr poly = dd_DDMatrix2Poly(matrix, &error);
    check_cddlib_error(error, "cone_to_cddlib_polyhedron");
    dd_MatrixPtr inequalities = dd_CopyInequalities(poly);
    assert(inequalities->representation == dd_Inequality);
    int num_inequalities = inequalities->rowsize;
    /* For each computed facet, obtain the set of input rays that are
       incident with the facet. */
    dd_SetFamilyPtr incidence = dd_CopyIncidence(poly);
    assert(incidence->setsize == num_rays + 1);
    assert(incidence->famsize == num_inequalities);
    /* Walk through all facets. */
    for (i = 0; i<num_inequalities; i++) {
      if (set_member(i + 1, inequalities->linset)) {
	/* We ignore equalities. */
      }
      else {
	/* If the extra ray is incident, we have a facet that is a wall
	   of the half-infinite cylinder over the cone -- ignore it.
	   Otherwise, the facet gives a simplicial cone of the
	   triangulation. */
	if (!set_member(num_rays + 1 /* 1-based */,
			incidence->set[i])) {
	  listCone *c = cone_from_ray_set(rays, incidence->set[i], cone->vertex);
	  /* Is a cone of the triangulation -- check it is simplicial */
	  int c_num_rays = set_card(incidence->set[i]);
	  if (c_num_rays > cone_dimension && !Parameters->nonsimplicial_subdivision) {
	    cerr << "Found non-simplicial cone (" << c_num_rays << "rays) "
		 << "in purported triangulation, triangulating it recursively." << endl;
	    /* In the refinement step, always fall back to using a
	       random height vector. */
	    listCone *ct = triangulate_cone_with_cddlib(c, Parameters,
							random_height, &Parameters->triangulation_max_height,
							cone_dimension);
	    freeCone(c);
	    triangulation = appendListCones(ct, triangulation);
	  }
	  else if (c_num_rays < cone_dimension) {
	    cerr << "Lower-dimensional cone in purported triangulation, should not happen."
		 << endl;
	    abort();
	  }
	  else {
	    c->rest = triangulation;
	    triangulation = c;
	  }
	}
      }
    }
    dd_FreeMatrix(inequalities);
    dd_FreeSetFamily(incidence);
    dd_FreeMatrix(matrix);
    dd_FreePolyhedra(poly);
    return triangulation;
  NEW_LIFTING:
    dd_FreeMatrix(inequalities);
    dd_FreeSetFamily(incidence);
    dd_FreePolyhedra(poly);
    freeListCone(triangulation); triangulation = NULL;
  } while (1);
}

listCone *
random_regular_triangulation_with_cddlib(listCone *cone,
					 BarvinokParameters *Parameters)
{
  return triangulate_cone_with_cddlib(cone, Parameters, random_height, &Parameters->triangulation_max_height,
				      Parameters->Number_of_Variables);
}

listCone *
refined_delone_triangulation_with_cddlib(listCone *cone,
					 BarvinokParameters *Parameters)
{
  return triangulate_cone_with_cddlib(cone, Parameters, delone_height, NULL,
				      Parameters->Number_of_Variables);
}
