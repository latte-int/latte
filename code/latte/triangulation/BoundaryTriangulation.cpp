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

listCone *
boundary_triangulation_of_cone_with_subspace_avoiding_facets
(listCone *cone, BarvinokParameters *Parameters)
{
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
    listCone *facet_triangulation
      = triangulate_cone_with_cddlib(facet_cone, Parameters,
				     delone_height, NULL,
				     Parameters->Number_of_Variables - 1);
    cout << "Triangulation of facet cone: " << lengthListCone(facet_triangulation)
	 << " simplicial cones." << endl;
    boundary_triangulation
      = appendListCones(facet_triangulation, boundary_triangulation);
  }
  /* Complete the cones with an interior ray. */
  vec_ZZ interior_ray_vector;
  interior_ray_vector.SetLength(Parameters->Number_of_Variables);
  for (i = 0; i<rays.size(); i++) {
    interior_ray_vector += rays[i]->first;
  }
  listCone *simplicial_cone;
  for (simplicial_cone = boundary_triangulation;
       simplicial_cone!=NULL;
       simplicial_cone = simplicial_cone->rest) {
    simplicial_cone->rays
      = appendVectorToListVector(interior_ray_vector, simplicial_cone->rays);
  }
  bool triangulation_ok = true;
  do {
    /* Pick a random interior ray vector. */
    int j;
    for (j = 0; j<Parameters->Number_of_Variables; j++) {
      interior_ray_vector[j] = uniform_random_number(1, 100);
    }
    cout << "Interior ray vector: " << interior_ray_vector << endl;
    /* Replace it in all cones and clear all facets. */
    for (simplicial_cone = boundary_triangulation;
	 simplicial_cone!=NULL;
	 simplicial_cone = simplicial_cone->rest) {
      simplicial_cone->rays->first = interior_ray_vector;
      if (simplicial_cone->facets != NULL) {
	freeListVector(simplicial_cone->facets);
	simplicial_cone->facets = NULL;
      }
    }
    /* Compute and check the facets. */
    boundary_triangulation
      = dualizeBackCones(boundary_triangulation, Parameters->Number_of_Variables);
    for (simplicial_cone = boundary_triangulation;
	 simplicial_cone!=NULL;
	 simplicial_cone = simplicial_cone->rest) {
      if (!rays_ok(simplicial_cone, Parameters->Number_of_Variables)) {
	triangulation_ok = false;
	break;
      }
    }
    boundary_triangulation
      = dualizeBackCones(boundary_triangulation, Parameters->Number_of_Variables);
  } while (!triangulation_ok);

  return boundary_triangulation;
}
