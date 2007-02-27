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
#include <vector>
#include <set>
#include <functional>

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
int verbosity = 1;

// Keeping track of the Hilbert basis candidates, to avoid duplicates

vector<int> *v = NULL;

std::set < vector<int> > known_hilbert_vectors;

static bool insert_hilbert_basis_element(int *vec)
{
  int numOfVars = params.Number_of_Variables;
  if (v == NULL) v = new vector<int>(numOfVars);

  int i;
  for (i = 0; i<numOfVars; i++) {
    (*v)[i] = vec[i];
  }
  
  std::set<vector<int> >::const_iterator where = known_hilbert_vectors.find(*v);
  if (where == known_hilbert_vectors.end()) {
    // Not known yet, so add it and print it. 
    known_hilbert_vectors.insert(*v);
    fprintVector(output, vec, numOfVars);
    fputs("\n", output);
    return true;
  }
  return false;
}

// The recursive decomposition and Hilbert basis computation.

class RecursiveNormalizer : public ConeConsumer {
public:
  int t_count;
  int t_level;
  int t_total;
  RecursiveNormalizer(int level) : t_count(0), t_total(0), t_level(level) {};
  int ConsumeCone(listCone *cone);
  void SetNumCones(size_t num_cones) { t_total = num_cones; }
};

static void
handle_cone(listCone *t, int t_count, int t_total, int level);

int
RecursiveNormalizer::ConsumeCone(listCone *cone)
{
  t_count++;
  handle_cone(cone, t_count, t_total, t_level);
  freeCone(cone);
  return 1; /* OK */
}

bool
cone_unimodular(listCone *cone, int numOfVars)
{
  int i;
  listVector *rays;
  mat_ZZ Mat;
  Mat.SetDims(numOfVars, numOfVars);
  rays=cone->rays;
  for(i = 0; i < numOfVars; i++) {
    Mat[i] = rays->first;
    rays = rays -> rest;
  }
  ZZ d = determinant(Mat);
  return abs(d) == 1;
}

static void
handle_cone(listCone *t, int t_count, int t_total, int level)
{
  params.Number_of_Variables = t->rays->first.length();
  
  int num_rays = lengthListVector(t->rays);

  if (num_rays == params.Number_of_Variables
      && cone_unimodular(t, params.Number_of_Variables)) return;

  int num_facets;
  int dimension = t->rays->first.length();
  Timer dualization_time("dualization", /*start_timer:*/ false);
  Timer zsolve_time("zsolve", /*start_timer:*/ false);

  if (verbosity > 0) {
    cout << "### " << "Level " << level << ": "
	 << "Cone " << t_count << " of at most " << t_total << ": "
	 << num_rays << " rays "
	 << "(dim " << dimension << ")";
    cout.flush();
  }

  //printCone(t, params.Number_of_Variables);

  /* Compute the facets of the cone. */
  dualization_time.start();
  dualizeCone(t, params.Number_of_Variables, &params); // computes and swaps
  dualizeCone(t, params.Number_of_Variables, &params); // just swaps back
  dualization_time.stop();
  num_facets = lengthListVector(t->facets);
  if (verbosity > 0) {
    cout << ", " << num_facets << " facets; "
	 << dualization_time;
  }
    
#if 0
  string current_cone_filename = filename + ".current_cone.triang";
  {
    ofstream current_cone_file(current_cone_filename.c_str());
    printConeToFile(current_cone_file, t, params.Number_of_Variables);
  }
#endif

  stats << level << "\t" << t_count << "\t"
	<< num_rays << "\t" << num_facets << "\t"
	<< t->determinant << "\t"
	<< dualization_time.get_seconds() << "\t";

  if (abs(t->determinant) == 1) {
    // simplicial, unimodular cone: Do nothing.
    stats << endl;
  }
  else if (num_facets < max_facets) {
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
    if (verbosity > 0) {
      cout << num_hilberts << " Hilbert basis elements; "
	   << zsolve_time;
    }
    
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
    
    //fprintVectorArray(output, ctx->Homs, false);
    //fprintVectorArray(output, ctx->Frees, false);
    int i;
    bool any_new = false;
    for (i = 0; i<ctx->Homs->Size; i++) {
      if (insert_hilbert_basis_element(ctx->Homs->Data[i])) any_new = true;
    }
    if (any_new) fflush(output);
    
    deleteZSolveContext(ctx, true);
    
    stats << zsolve_time.get_seconds() << "\t"
	  << num_hilberts << endl;
  }
  else {
    stats << endl;
    if (verbosity > 0) {
      cout << "Too many facets, subdividing..." << endl;
    }
    RecursiveNormalizer rec(level + 1);
    triangulateCone(t, params.Number_of_Variables, &params, rec);
  }
}

  //  cout << "Subdivision has " << lengthListCone(triang) << " cones." << endl;

static void open_output_and_stats()
{
  string output_filename = filename + ".hil";
  cout << "Output goes to file `" << output_filename << "'..." << endl;
  output = fopen(output_filename.c_str(), "w");

  string stats_filename = filename + ".stats";
  cout << "Cone statistics go to file `" << stats_filename << "'..." << endl;
  stats.open(stats_filename.c_str());
  stats << "# Level\tIndex\tRays\tFacets\tDet\tDualize\tZSolve\tHilberts" << endl;

  // Redirect 4ti2/qsolve output.
  string fortytwolog_filename = filename + ".4ti2log";
  static ofstream fortytwolog(fortytwolog_filename.c_str());
  _4ti2_::out = &fortytwolog;
}

static listCone *
read_cone_cdd_format(string &filename)
{
  FILE *in = fopen(filename.c_str(), "r");
  if (in == NULL) {
    cerr << "normaliz: Unable to open CDD-style input file " << filename << endl;
    exit(1);
  }
  dd_MatrixPtr M;
  dd_ErrorType err=dd_NoError;
  M = dd_PolyFile2Matrix(in, &err);
  if (err!=dd_NoError) {
    cerr << "normaliz: Parse error in CDD-style input file " << filename << endl;
    exit(1);
  }
  listCone *cone = cddlib_matrix_to_cone(M);
  dd_FreeMatrix(M);
  return cone;
}

static void
normalize_from_triang_file(string &triang_filename, size_t num_cones = 0)
{
  open_output_and_stats();
  ifstream triang_file(triang_filename.c_str());
  RecursiveNormalizer normalizer(/*level:*/ 1);
  normalizer.t_total = num_cones;
  readListConeFromFile(triang_file, normalizer);
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
	 << "             --quiet                       Do not show much output" << endl
	 << "             --no-triang-file              Do not create a .triang file" << endl
      ;
    exit(1);
  }

  cout << "This is pre-NORMALIZ that actually handles the advertised options." << endl;
  
  listCone *cone;
  listCone *triang;
  bool create_triang_file = true;
  
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
      else if (strcmp(argv[i], "--quiet") == 0) {
	verbosity = 0;
      }
      else if (strcmp(argv[i], "--no-triang-file") == 0) {
	create_triang_file = false;
      }
      else {
	cerr << "Unknown option " << argv[i] << endl;
	exit(1);
      }
    }
  }

  filename = argv[argc-1];

  string triang_filename;

  if (strlen(filename.c_str()) > 7
      && strcmp(filename.c_str() + strlen(filename.c_str()) - 7, ".triang") == 0) {
    
    triang_filename = filename;
    normalize_from_triang_file(triang_filename);
  }
  else {
    // Read a cone and triangulate it.
    if (strlen(filename.c_str()) > 4 && strcmp(filename.c_str() + strlen(filename.c_str()) - 4, ".ext") == 0) {
      /* Input in CDD format. */
      cone = read_cone_cdd_format(filename);
      params.Number_of_Variables = cone->rays->first.length();
    }
    else {
      cerr << "normaliz: Want a .ext file." << endl;
      exit(1);
    }

    if (create_triang_file) {
      // Create a .triang file, then read it back in again one-by-one
      // and feed it to the normalizer.
      triang_filename = filename + ".triang";
      size_t num_cones;
      {
	PrintingConeConsumer triang_file_writer(triang_filename);
	triangulateCone(cone, params.Number_of_Variables, &params, triang_file_writer);
	freeCone(cone);
	num_cones = triang_file_writer.cone_count;
	cout << "Printed triangulation to file `" << triang_filename << "'." << endl;
      }
      normalize_from_triang_file(triang_filename, num_cones);
    }
    else {
      // Triangulate and feed cones one-by-one to the normalizer.
      open_output_and_stats();
      RecursiveNormalizer normalizer(/*level:*/ 1);
      triangulateCone(cone, params.Number_of_Variables, &params, normalizer);
      freeCone(cone);
    }
  }

  return 0;
}

