/* triangulate.cpp -- Compute triangulations.

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

#include <cstdlib>
#include "config.h"
#include "triangulation/triangulate.h"
#include "triangulation/RegularTriangulationWithCdd.h"
#ifdef HAVE_CDDLIB
#include "triangulation/RegularTriangulationWithCddlib.h"
#endif
#ifdef HAVE_EXPERIMENTS
#  include "triangulation/BoundaryTriangulation.h"
#endif
#if defined(HAVE_TOPCOM_LIB) || defined(HAVE_TOPCOM_BIN)
#  include "triangulation/TriangulationWithTOPCOM.h"
#endif
#if defined(HAVE_FORTYTWO_LIB)
#  include "triangulation/RegularTriangulationWith4ti2.h"
#endif
#include "print.h"

BarvinokParameters::TriangulationType
triangulation_type_from_name(const char *name)
{
  if (strcmp(name, "cdd") == 0) return BarvinokParameters::RegularTriangulationWithCdd;
  else if (strcmp(name, "cddlib") == 0) return BarvinokParameters::RegularTriangulationWithCddlib;
  else if (strcmp(name, "delone") == 0 || strcmp(name, "delaunay") == 0)
    return BarvinokParameters::DeloneTriangulationWithCddlib;
  else if (strcmp(name, "topcom") == 0) return BarvinokParameters::PlacingTriangulationWithTOPCOM;
  else if (strcmp(name, "boundary") == 0) return BarvinokParameters::SubspaceAvoidingBoundaryTriangulation;
  else if (strcmp(name, "4ti2") == 0) return BarvinokParameters::RegularTriangulationWith4ti2;
  else {
    cerr << "Unknown triangulation type name: " << name << endl;
    exit(1);
  }
}

bool
parse_standard_triangulation_option(const char *arg,
				    BarvinokParameters *params)
{
  if (strncmp(arg, "--triangulation=", 16) == 0) {
    params->triangulation = triangulation_type_from_name(arg + 16);
  }
  else if (strncmp(arg, "--triangulation-max-height=", 27) == 0) {
    params->triangulation_max_height = atoi(arg + 27);
  }
  else if (strncmp(arg, "--nonsimplicial-subdivision", 9) == 0) {
    params->nonsimplicial_subdivision = true;
  }
  else if (strncmp(arg, "--triangulation-bias=", 21) == 0) {
    params->triangulation_bias = atoi(arg + 21);
  }
  else return false;
  return true;
}

void
show_standard_triangulation_options(ostream &stream)
{
  stream << "Triangulation options:" << endl
	 << "  --triangulation={cddlib,4ti2,topcom,...}" << endl
	 << "  --triangulation-max-height=HEIGHT        Use a uniform distribution of height from 1 to HEIGHT." << endl
         << "  --triangulation-bias=PERCENTAGE          Use a non-uniform distribution of heights 1 and 2." << endl;
}

listCone *
triangulateCone(listCone *cone, int numOfVars,
		BarvinokParameters *params)
{
  cout << "Triangulating cone... " << flush;
  params->triangulate_time.start();
  CollectingConeConsumer ccc;
  triangulateCone(cone, numOfVars, params, ccc);
  cout << "done." << endl;
  params->triangulate_time.stop();
  return ccc.Collected_Cones;
}

void
triangulateCone(listCone *cone, int numOfVars,
		BarvinokParameters *params,
		ConeConsumer &consumer)
{
  if (numOfVars == lengthListVector(cone->rays)) {
    // Already simplicial.
    consumer.ConsumeCone(copyCone(cone));
    return;
  }
  switch(params->triangulation) {
  case BarvinokParameters::RegularTriangulationWithCdd:
    triangulate_cone_with_cdd(cone, params, consumer);
    break;
  case BarvinokParameters::RegularTriangulationWithCddlib:
#ifdef HAVE_CDDLIB
    random_regular_triangulation_with_cddlib(cone, params, consumer);
#else
    cerr << "RegularTriangulationWithCddlib not compiled in, sorry."
	 << endl;
    exit(1);
#endif
    break;
  case BarvinokParameters::DeloneTriangulationWithCddlib:
#ifdef HAVE_CDDLIB
    refined_delone_triangulation_with_cddlib(cone, params, consumer);
#else
    cerr << "DeloneTriangulationWithCddlib not compiled in, sorry."
	 << endl;
    exit(1);
#endif
    break;
  case BarvinokParameters::SubspaceAvoidingBoundaryTriangulation:
#ifdef HAVE_EXPERIMENTS
    boundary_triangulation_of_cone_with_subspace_avoiding_facets
      (cone, params, consumer);
#else
    cerr << "SubspaceAvoidingBoundaryTriangulation not compiled in, sorry."
	 << endl;
    exit(1);
#endif
    break;
  case BarvinokParameters::PlacingTriangulationWithTOPCOM:
#if defined(HAVE_TOPCOM_LIB) || defined(HAVE_TOPCOM_BIN)
    triangulate_cone_with_TOPCOM(cone, numOfVars, consumer);
#else
    cerr << "PlacingTriangulationWithTOPCOM not compiled in, sorry."
	 << endl;
    exit(1);
#endif
    break;
  case BarvinokParameters::RegularTriangulationWith4ti2:
#if defined(HAVE_FORTYTWO_LIB)
    random_regular_triangulation_with_4ti2(cone, params, consumer);
#else
    cerr << "RegularTriangulationWith4ti2 not compiled in, sorry."
	 << endl;
    exit(1);
#endif
    break;
  default:
    cerr << "Unknown triangulation method." << endl;
    exit(1);
  }
}


TriangulatingConeTransducer::TriangulatingConeTransducer(BarvinokParameters *some_parameters)
  : params(some_parameters)
{
}

int TriangulatingConeTransducer::ConsumeCone(listCone *cone)
{
  triangulateCone(cone, cone->rays->first.length(), params, 
		  *consumer);
  freeCone(cone);
  return 1; /* OK */
}

