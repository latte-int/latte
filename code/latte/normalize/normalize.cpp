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
#include "ReadSubcones.h"
#include "ReadLatteStyle.h"

#include "normalize/ReductionTest.h"

// from 4ti2:
#include "groebner/Globals.h"

using namespace std;

string hil_filename;
IncrementalVectorFileWriter *hil_file_writer = NULL;
ReductionTest *reduction_test = NULL;

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
    if (!reduction_test->IsReducible(*v)) {
      // Not known yet, so add it and print it. 
      known_hilbert_vectors.insert(*v);
      hil_file_writer->WriteVector(*v);
      return true;
    }
  }
  return false;
}

// The recursive decomposition and Hilbert basis computation.

class RecursiveNormalizer : public ConeConsumer {
public:
  int t_count;
  int t_level;
  int t_total;
  RecursiveNormalizer(int level) : t_count(0), t_level(level), t_total(0) {};
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

#include <setjmp.h>
int zsolve_time_limit = 0;
jmp_buf timelimit_jmp_buf;

static void check_timelimit(struct zsolvecontext_t *zsolve)
{
  longjmp(timelimit_jmp_buf, 1);
}

static void
handle_cone(listCone *t, int t_count, int t_total, int level)
{
  params.Number_of_Variables = t->rays->first.length();
  if (hil_file_writer == NULL) {
    hil_file_writer
      = new IncrementalVectorFileWriter(hil_filename, params.Number_of_Variables);
    hil_file_writer->UpdateNumVectors();
  }
  
  int num_rays = lengthListVector(t->rays);

  if (num_rays == params.Number_of_Variables
      && cone_unimodular(t, params.Number_of_Variables)) return;

  int num_facets;
  int dimension = t->rays->first.length();
  Timer dualization_time("dualization", /*start_timer:*/ false);
  Timer zsolve_time("zsolve", /*start_timer:*/ false);

  if (verbosity > 0) {
    cerr << "### " << "Level " << level << ": "
	 << "Cone " << t_count << " of at most " << t_total << ": "
	 << num_rays << " rays "
	 << "(dim " << dimension << ")";
    cerr.flush();
  }

  //printCone(t, params.Number_of_Variables);

  /* Compute the facets of the cone. */
  dualization_time.start();
  dualizeCone(t, params.Number_of_Variables, &params); // computes and swaps
  dualizeCone(t, params.Number_of_Variables, &params); // just swaps back
  dualization_time.stop();
  num_facets = lengthListVector(t->facets);
  if (verbosity > 0) {
    cerr << ", " << num_facets << " facets; "
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
				      /* logcallback: */ NULL,
				      /* backupEvent: */ (zsolve_time_limit
							  ? check_timelimit
							  : NULL));
    ctx->BackupTime = zsolve_time_limit;
    deleteLinearSystem(ls);
    zsolve_time.start();
    if (zsolve_time_limit == 0 || setjmp(timelimit_jmp_buf) == 0) {
      zsolveSystem(ctx, /*appendnegatives:*/ true);
      zsolve_time.stop();

      int num_hilberts = ctx->Homs->Size;
      if (verbosity > 0) {
	cerr << num_hilberts << " Hilbert basis elements; "
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
      if (any_new)
	hil_file_writer->UpdateNumVectors();
    
      deleteZSolveContext(ctx, true);
    
      stats << zsolve_time.get_seconds() << "\t"
	    << num_hilberts << endl;
    }
    else {
      /* Longjmp target -- the timelimit was reached. */
      zsolve_time.stop();
      stats << zsolve_time.get_seconds() << endl;
      deleteZSolveContext(ctx, true);
      if (verbosity > 0) {
	cerr << "Spent too much time in zsolve, subdividing..." << endl;
      }
      RecursiveNormalizer rec(level + 1);
      triangulateCone(t, params.Number_of_Variables, &params, rec);
    }
  }
  else {
    stats << endl;
    if (verbosity > 0) {
      cerr << "Too many facets, subdividing..." << endl;
    }
    RecursiveNormalizer rec(level + 1);
    triangulateCone(t, params.Number_of_Variables, &params, rec);
  }
}

static void open_output_and_stats()
{
  hil_filename = filename + ".hil";
  cerr << "Output goes to file `" << hil_filename << "'..." << endl;

  string stats_filename = filename + ".stats";
  cerr << "Cone statistics go to file `" << stats_filename << "'..." << endl;
  stats.open(stats_filename.c_str());
  stats << "# Level\tIndex\tRays\tFacets\tDet\tDualize\tZSolve\tHilberts" << endl;

  // Redirect 4ti2/qsolve output.
  string fortytwolog_filename = "/dev/null"; //filename + ".4ti2log";
  static ofstream fortytwolog(fortytwolog_filename.c_str());
  _4ti2_::out = &fortytwolog;
}

static listCone *
read_cone_cdd_format(const string &filename)
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

static listCone *
read_cone_4ti2_format(const string &filename)
{
  dd_MatrixPtr M = ReadLatteStyleMatrix(filename.c_str(), /*vrep:*/true, /*homogenize:*/true);
  listCone *cone = cddlib_matrix_to_cone(M);
  dd_FreeMatrix(M);
  return cone;
}

ReductionTestFactory reduction_test_factory;

static void
usage()
{
    cerr << "usage: normaliz [OPTIONS] { CDD-EXT-FILE.ext | LATTE-TRIANG-FILE.triang | 4TI2-STYLE-FILE.{rays,tra} } " << endl;
    cerr << "Options are: " << endl
	 << "  --dualization={cdd,4ti2}" << endl;
    show_standard_triangulation_options(cerr);
    cerr << "  --nonsimplicial-subdivision              [Default]" << endl
	 << "  --max-facets=N                           Subdivide further if more than N facets" << endl
         << "  --zsolve-time-limit=SECONDS              Subdivide further if computation of Hilbert" << endl
	 << "                                           basis took longer than this number of seconds." << endl
	 << "  --quiet                                  Do not show much output" << endl
	 << "  --no-triang-file                         Do not create a .triang file" << endl
         << "  --subcones=INPUT-FILE.subcones           Read list of subcone indicators to handle" << endl
         << "  --output-subcones=OUTPUT-FILE.subcones   Write a list of toplevel subcones" << endl
	 << "  --only-triangulate                       Only triangulate, don't normalize" << endl
         << "  --no-initial-triangulation               Don't compute an initial triangulation," << endl
         << "                                           start recursive normalizer on input." << endl
      ;
    reduction_test_factory.show_options(cerr);
}

int main(int argc, char **argv)
{
  if (argc < 2) {
    usage();
    exit(1);
  }

  cerr << "This is the joint LattE/4ti2 almost-NORMALIZ program." << endl;

  listCone *cone;
  bool create_triang_file = true;
  string subcones_filename;
  bool have_subcones = false;
  string output_subcones_filename;
  bool have_output_subcones = false;
  bool normalize = true;
  bool triangulate_toplevel = true;
  
  params.triangulation = BarvinokParameters::RegularTriangulationWith4ti2;
  params.dualization = BarvinokParameters::DualizationWith4ti2;
  params.nonsimplicial_subdivision = true;

  
  {
    int i;
    for (i = 1; i<argc; i++) {
      if (parse_standard_triangulation_option(argv[i], &params)) {}
      else if (parse_standard_dualization_option(argv[i], &params)) {}
      else if (reduction_test_factory.parse_option(argv[i])) {}
      else if (strncmp(argv[i], "--max-facets=", 13) == 0) {
	max_facets = atoi(argv[i] + 13);
      }
      else if (strcmp(argv[i], "--quiet") == 0) {
	verbosity = 0;
      }
      else if (strcmp(argv[i], "--no-triang-file") == 0) {
	create_triang_file = false;
      }
      else if (strncmp(argv[i], "--zsolve-time-limit=", 20) == 0) {
	zsolve_time_limit = atoi(argv[i] + 20);
      }
      else if (strncmp(argv[i], "--subcones=", 11) == 0) {
	subcones_filename = string(argv[i] + 11);
	have_subcones = true;
      }
      else if (strncmp(argv[i], "--output-subcones=", 18) == 0) {
	output_subcones_filename = string(argv[i] + 18);
	have_output_subcones = true;
      }
      else if (strncmp(argv[i], "--only-triangulate", 6) == 0) {
	normalize = false;
      }
      else if (strncmp(argv[i], "--no-initial-triangulation", 12) == 0) {
	triangulate_toplevel = false;
      }
      else if (strncmp(argv[i], "--help", 6) == 0) {
	usage();
	exit(0);
      }
      else if (strncmp(argv[i], "--", 2) == 0) {
	cerr << "Unknown option " << argv[i] << endl;
	exit(1);
      }
      else if (i == argc-1) {
	filename = argv[i];
      }
      else {
	cerr << "Unexpected argument " << argv[i] << endl;
	exit(1);
      }
    }
  }

  if (normalize)
    reduction_test = reduction_test_factory.CreateReductionTest();

  if (filename.length() == 0) {
    cerr << "Missing input filename" << endl;
    exit(1);
  }

  string triang_filename;

  // Input handling.

  ConeProducer *producer = NULL;
  
  if (strlen(filename.c_str()) > 7
      && strcmp(filename.c_str() + strlen(filename.c_str()) - 7, ".triang") == 0) {
    if (have_subcones) {
      cerr << "Cannot use both a triangulation file and a subcones file." << endl;
      exit(1);
    }
    triang_filename = filename;
    producer = new ListConeReadingConeProducer(filename);
    triangulate_toplevel = false;
  }
  else {
    // Read a cone.
    if (strlen(filename.c_str()) > 4 && strcmp(filename.c_str() + strlen(filename.c_str()) - 4, ".ext") == 0) {
      /* Input in CDD format. */
      cone = read_cone_cdd_format(filename);
    }
    else {
      /* Try to read a 4ti2-style file. */
      cerr << "Trying to read `" << filename << "' as a list of rays in 4ti2-style format." << endl;
      cone = read_cone_4ti2_format(filename);
    }
    params.Number_of_Variables = cone->rays->first.length();
    if (have_subcones) {
      // Also a subcones file given.
      producer = new SubconeReadingConeProducer(cone, subcones_filename);
    }
    else {
      producer = new SingletonConeProducer(copyCone(cone));
    }
  }

  if (triangulate_toplevel) {
    // Send input cones through triangulation.
    ConeTransducer *triangulator = new TriangulatingConeTransducer(&params);
    producer = compose(producer, triangulator);
  }

  // Computation handling.   

  if (triangulate_toplevel) {
    if (have_output_subcones) {
      // Create a subcones file, then read it back in again one-by-one
      // and feed it to the normalizer.
      size_t num_cones;
      {
	SubconePrintingConeConsumer subcone_file_writer(cone, 
							output_subcones_filename);
	producer->Produce(subcone_file_writer);
	num_cones = subcone_file_writer.cone_count;
	cerr << "Printed triangulation to subcones file `" << output_subcones_filename << "'." << endl;
      }
      producer = new SubconeReadingConeProducer(cone, output_subcones_filename, num_cones);
    }
    else if (create_triang_file) {
      // Create a .triang file, then read it back in again one-by-one
      // and feed it to the normalizer.
      size_t num_cones;
      triang_filename = filename + ".triang";
      {
	PrintingConeConsumer triang_file_writer(triang_filename);
	producer->Produce(triang_file_writer);
	num_cones = triang_file_writer.cone_count;
	cerr << "Printed triangulation to file `" << triang_filename << "'." << endl;
      }
      producer = new ListConeReadingConeProducer(triang_filename, num_cones);
    }
  }

  if (normalize) {
    // Feed cones one-by-one to the normalizer.
    open_output_and_stats();
    RecursiveNormalizer normalizer(/*level:*/ 1);
    producer->Produce(normalizer);
    if (hil_file_writer) {
      // Update the number of vectors, close the file.
      delete hil_file_writer;
    }
    delete reduction_test;
  }

  return 0;
}

