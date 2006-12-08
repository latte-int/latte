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
#include "latte_cddlib.h"
#include "print.h"

using namespace std;

static void
check_rays(listCone *cone, int numOfVars)
{
  listVector *ray;
  for (ray = cone->rays; ray!=NULL; ray=ray->rest) {
    if (ray->first[numOfVars - 1] == 0) {
      cerr << "The following dualized-back cone has bad rays." << endl;
      printListCone(cone, numOfVars);
      abort();
    }
  }
}

listCone *
boundary_triangulation_of_cone_with_subspace_avoiding_facets
(listCone *cone, BarvinokParameters *Parameters)
{
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
  }
  cerr << "Code unfinished." << endl;
  abort();
}
