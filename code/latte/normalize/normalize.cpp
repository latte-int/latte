/* normalize.cpp -- Re-implementation of NORMALIZ
	       
   Copyright 2007 Matthias Koeppe

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

#include <cstdio>
#include <iostream>
#include <string>
#include <cctype>
#include "print.h"
#include "triangulation/triangulate.h"
#include "dual.h"
#include "latte_cddlib.h"
#include "latte_4ti2.h"
#include "latte_4ti2_zsolve.h"
#include "timing.h"
// from 4ti2:
#include "Globals.h"

using namespace std;


FILE *output;
string filename;
ofstream stats;
BarvinokParameters params;
int max_facets = INT_MAX;

static void
handle_list_cone(listCone *triang, int level);

static void
handle_cone(listCone *t, int t_count, int t_total, int level)
{
  int num_rays = lengthListVector(t->rays);
  int num_facets;
  int dimension = t->rays->first.length();
  Timer dualization_time("dualization", /*start_timer:*/ false);
  Timer zsolve_time("zsolve", /*start_timer:*/ false);
    
  cout << "### " << "Level " << level << ": "
       << "Cone " << t_count << " of " << t_total << ": "
       << num_rays << " rays "
       << "(dim " << dimension << ")";
  cout.flush();

  //printCone(t, params.Number_of_Variables);

  /* Compute the facets of the cone. */
  dualization_time.start();
  dualizeCone(t, params.Number_of_Variables, &params); // computes and swaps
  dualizeCone(t, params.Number_of_Variables, &params); // just swaps back
  dualization_time.stop();
  num_facets = lengthListVector(t->facets);
  cout << ", " << num_facets << " facets; "
       << dualization_time;
    
#if 0
  string current_cone_filename = filename + ".current_cone.triang";
  {
    ofstream current_cone_file(current_cone_filename.c_str());
    printConeToFile(current_cone_file, t, params.Number_of_Variables);
  }
#endif

  stats << level << "\t" << t_count << "\t"
	<< num_rays << "\t" << num_facets << "\t"
	<< dualization_time.get_seconds() << "\t";

  if (num_facets < max_facets) {
    // Use zsolve to compute the Hilbert basis.
    
    //printCone(t, params.Number_of_Variables);
    LinearSystem ls
      = facets_to_4ti2_zsolve_LinearSystem(t->facets, params.Number_of_Variables);

    //printLinearSystem(ls);

    ZSolveContext ctx
      = createZSolveContextFromSystem(ls, NULL/*LogFile*/, 0/*OLogging*/, -1/*OVerbose*/,
				      /* logcallback: */ NULL, NULL/*backupEvent*/);
    deleteLinearSystem(ls);
    zsolve_time.start();
    zsolveSystem(ctx, /*appendnegatives:*/ true);
    zsolve_time.stop();

    int num_hilberts = ctx->Homs->Size;
    cout << num_hilberts << " Hilbert basis elements; "
	 << zsolve_time;
    
    if (num_hilberts < params.Number_of_Variables) {
      // Sanity check.
      cerr << "Too few Hilbert basis elements " << endl;
      printCone(t, params.Number_of_Variables);
      LinearSystem ls
	= facets_to_4ti2_zsolve_LinearSystem(t->facets, params.Number_of_Variables);
      printLinearSystem(ls);
      fprintf(stdout, "%d %d\n\n", ctx->Homs->Size + ctx->Frees->Size, ctx->Homs->Variables);
      fprintVectorArray(stdout, ctx->Homs, false);
      fprintVectorArray(stdout, ctx->Frees, false);
      abort();
    }

    fprintVectorArray(output, ctx->Homs, false);
    fprintVectorArray(output, ctx->Frees, false);
    deleteZSolveContext(ctx, true);
    
    stats << zsolve_time.get_seconds() << "\t"
	  << num_hilberts << endl;
  }
  else {
    stats << endl;
    cout << "Too many facets, subdividing..." << endl;
    listCone *triang = triangulateCone(t, params.Number_of_Variables, &params);
    handle_list_cone(triang, level + 1);
    freeListCone(triang);
  }
}

static void
handle_list_cone(listCone *triang, int level)
{
  listCone *t;
  int t_count = 1;
  int t_total = lengthListCone(triang);
  cout << "Subdivision has " << lengthListCone(triang) << " cones." << endl;
  for (t = triang; t!=NULL; t = t->rest, t_count++)
    handle_cone(t, t_count, t_total, level);
  cout << "Done." << endl;
}
  
int main(int argc, char **argv)
{
  if (argc < 2) {
    cerr << "usage: normaliz [OPTIONS] { CDD-EXT-FILE.ext | LATTE-TRIANG-FILE.triang } " << endl;
    cerr << "Options are: --triangulation={cddlib,4ti2,topcom,...}" << endl
	 << "             --triangulation-max-height=HEIGHT" << endl
	 << "             --nonsimplicial-subdivision" << endl
	 << "             --dualization={cdd,4ti2}" << endl
	 << "             --max-facets=N                Subdivide further if more than N facets" << endl
      ;
    exit(1);
  }

  cout << "This is pre-NORMALIZ that actually handles the advertised options." << endl;
  
  listCone *cone;
  listCone *triang;
  
  params.triangulation = BarvinokParameters::RegularTriangulationWith4ti2;
  params.dualization = BarvinokParameters::DualizationWith4ti2;

  {
    int i;
    for (i = 1; i<argc-1; i++) {
      if (parse_standard_triangulation_option(argv[i], &params)) {}
      else if (parse_standard_dualization_option(argv[i], &params)) {}
      else if (strncmp(argv[i], "--max-facets=", 13) == 0) {
	max_facets = atoi(argv[i] + 13);
      }
      else {
	cerr << "Unknown option " << argv[i] << endl;
	exit(1);
      }
    }
  }

  filename = argv[argc-1];

  if (strlen(argv[argc-1]) > 7 && strcmp(argv[argc-1] + strlen(argv[argc-1]) - 7, ".triang") == 0) {
    ifstream tf(argv[argc-1]);
    triang = readListConeFromFile(tf);
    if (triang == NULL) {
      cerr << "Failed to read a triangulation." << endl;
      exit(1);
    }
    params.Number_of_Variables = triang->rays->first.length();
  }
  else {
    // Read a cone and triangulate it.
    if (strlen(argv[argc-1]) > 4 && strcmp(argv[argc-1] + strlen(argv[argc-1]) - 4, ".ext") == 0) {
      /* Input in CDD format. */
      FILE *in = fopen(argv[argc-1], "r");
      if (in == NULL) {
	cerr << "normaliz: Unable to open CDD-style input file " << argv[argc-1] << endl;
	exit(1);
      }
      dd_MatrixPtr M;
      dd_ErrorType err=dd_NoError;
      M = dd_PolyFile2Matrix(in, &err);
      if (err!=dd_NoError) {
	cerr << "normaliz: Parse error in CDD-style input file " << argv[argc-1] << endl;
	exit(1);
      }
      cone = cddlib_matrix_to_cone(M);
      dd_FreeMatrix(M);
      params.Number_of_Variables = cone->rays->first.length();
    }
    else {
      cerr << "normaliz: Want a .ext file." << endl;
      exit(1);
    }
    triang = triangulateCone(cone, params.Number_of_Variables, &params);
    freeCone(cone);
    string triang_filename = filename + ".triang";
    printListConeToFile(triang_filename.c_str(), triang, params.Number_of_Variables);
    cout << "Printed triangulation to file `" << triang_filename << "'." << endl;
  }
  
  string output_filename = filename + ".hil";
  cout << "Output goes to file `" << output_filename << "'..." << endl;
  output = fopen(output_filename.c_str(), "w");

  string stats_filename = filename + ".stats";
  cout << "Cone statistics go to file `" << stats_filename << "'..." << endl;
  stats.open(stats_filename.c_str());
  stats << "# Level\tIndex\tRays\tFacets\tDualize\tZSolve\tHilberts" << endl;

  // Redirect 4ti2/qsolve output.
  string fortytwolog_filename = filename + ".4ti2log";
  ofstream fortytwolog(fortytwolog_filename.c_str());
  _4ti2_::out = &fortytwolog;

  handle_list_cone(triang, 1);
    
  freeListCone(triang);
  return 0;
}

