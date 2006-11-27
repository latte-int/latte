/* triangulate.cpp -- Compute triangulations.

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

#include <cstdlib>
#include "config.h"
#include "triangulation/triangulate.h"
#include "triangulation/RegularTriangulationWithCdd.h"
#ifdef HAVE_EXPERIMENTS
#include "triangulation/RecursiveTriangulation.h"
#endif

listCone *
triangulateCone(listCone *cone, int numOfVars,
		BarvinokParameters *params)
{
  listCone *result;
  params->triangulate_time.start();
  switch(params->triangulation) {
  case BarvinokParameters::RegularTriangulationWithCdd:
    result = triangulate_cone_with_cdd(cone, params);
    break;
  case BarvinokParameters::SubspaceAvoidingRecursiveTriangulation:
#ifdef HAVE_EXPERIMENTS
    result = triangulate_cone_recursively_with_subspace_avoiding_facets
      (cone, numOfVars);
#else
    cerr << "SubspaceAvoidingRecursiveTriangulation not compiled in, sorry."
	 << endl;
    exit(1);
#endif
    break;
  default:
    cerr << "Unknown triangulation method." << endl;
    exit(1);
  }
  params->triangulate_time.stop();
  return result;
}
