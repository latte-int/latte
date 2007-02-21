/* test-triangulation.cpp -- Stand-alone program for triangulation
	       
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

#include <iostream>
#include <string>
#include <cctype>
#include "print.h"
#include "triangulate.h"
#include "dual.h"
#include "TriangulationWithTOPCOM.h"
#include "latte_cddlib.h"

using namespace std;

static void
print_triangulation(listCone *triang, int Number_of_Variables)
{
#if 0
  listCone *t;
  for (t = triang; t!=NULL; t = t->rest) {
    computeDetAndFacetsOfSimplicialCone(t, Number_of_Variables);
  }
#endif
//   cout << "*** Triangulation:" << endl;
//   printListCone(triang, Number_of_Variables);
  cout << "Triangulation (with " << lengthListCone(triang)
       << " cones) printed to file `triangulation'." << endl;
  printListConeToFile("triangulation", triang, Number_of_Variables);
}

int main(int argc, char **argv)
{
  if (argc < 2) {
    cerr << "usage: triangulate [OPTIONS] [LATTE-CONE-FILE | CDD-EXT-FILE.ext ] " << endl;
    cerr << "Options are: --triangulation={cddlib,4ti2,topcom,...}" << endl
	 << "             --triangulation-max-height=HEIGHT" << endl
	 << "             --nonsimplicial-subdivision" << endl;
    exit(1);
  }
  BarvinokParameters params;
  listCone *cone;

  params.triangulation = BarvinokParameters::RegularTriangulationWithCddlib;

  {
    int i;
    for (i = 1; i<argc-1; i++) {
      if (parse_standard_triangulation_option(argv[i], &params)) {}
      else {
	cerr << "Unknown option " << argv[i] << endl;
	exit(1);
      }
    }
  }
  
  if (strlen(argv[argc-1]) > 4 && strcmp(argv[argc-1] + strlen(argv[argc-1]) - 4, ".ext") == 0) {
    /* Input in CDD format. */
    FILE *in = fopen(argv[argc-1], "r");
    if (in == NULL) {
      cerr << "triangulate: Unable to open CDD-style input file " << argv[argc-1] << endl;
      exit(1);
    }
    dd_MatrixPtr M;
    dd_ErrorType err=dd_NoError;
    M = dd_PolyFile2Matrix(in, &err);
    if (err!=dd_NoError) {
      cerr << "triangulate: Parse error in CDD-style input file " << argv[argc-1] << endl;
      exit(1);
    }
    cone = cddlib_matrix_to_cone(M);
  }
  else {
    /* Input in LattE cone format. */
    ifstream in(argv[argc-1]);
    if (in.bad()) {
      cerr << "triangulate: Unable to open " << argv[argc-1] << endl;
      exit(1);
    }
    cone = readConeFromFile(in);
    if (!cone) {
      cerr << "triangulate: Parse error in file." << endl;
      exit(1);
    }
  }

  if (cone->rays == NULL) {
    cerr << "triangulate: No rays." << endl;
    exit(2);
  }
  params.Number_of_Variables = cone->rays->first.length();
  cone->vertex = new Vertex(new rationalVector(params.Number_of_Variables));

#if 0
  // Compute facets.
  cone = dualizeCones(cone, params.Number_of_Variables);
  cone = dualizeBackCones(cone, params.Number_of_Variables);
#endif
  cout << "*** Input cone:" << endl;
  printListCone(cone, params.Number_of_Variables);

#if 0
  list<listCone *>all_triangulations
    = all_triangulations_of_cone_with_TOPCOM(cone, params.Number_of_Variables);
  for (list<listCone *>::iterator i = all_triangulations.begin();
       i != all_triangulations.end();
       ++i)
    print_triangulation(*i, params.Number_of_Variables);
#else

  listCone *triang
    = triangulateCone(cone, params.Number_of_Variables, &params);
  print_triangulation(triang, params.Number_of_Variables);

#endif
  return 0;
}
