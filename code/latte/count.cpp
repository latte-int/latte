/* count.cpp -- Master program

   Copyright 2002, 2003 Raymond Hemmecke, Ruriko Yoshida
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

#include <string.h>
#include <stdio.h>
#include <cassert>

#include "config.h"
#include "latte_ntl_integer.h"
#include "barvinok/dec.h"
#include "barvinok/barvinok.h"
#include "barvinok/Triangulation.h"
#include "vertices/cdd.h"
#include "genFunction/maple.h"
#include "genFunction/piped.h"
#include "cone.h"
#include "dual.h"
#include "RudyResNTL.h"
#include "Residue.h"
#include "Grobner.h"
//  #include "jesus.h"
#include "preprocess.h"
#include "print.h"
#include "ramon.h"
#include "rational.h"
#include "timing.h"
#include "flags.h"
//#include "testing.h"
#include "IntegralHull.h"
#include "ExponentialSubst.h"
#include "latte_random.h"
#include "Irrational.h"
#include "ExponentialEhrhart.h"
#include "triangulation/triangulate.h"
#ifdef HAVE_EXPERIMENTS
#include "ExponentialApprox.h"
#include "TrivialSubst.h"
#endif

#include "banner.h"
#include "convert.h"
#include "latte_system.h"
#include "Polyhedron.h"
#include "ReadPolyhedron.h"

static void
usage(const char *progname)
{
  cerr << "usage: " << progname << " [OPTIONS...] " << "INPUTFILE" << endl;
  cerr << "Type `" << progname << " --help' "
       << "for a list of options and input specifications." << endl;
}

/* ----------------------------------------------------------------- */
int
main(int argc, char *argv[])
{
  int i;
  unsigned int flags = 0, print_flag = 0, output_cone = 0;
  char printfile[127],
    Save_Tri[127],
    Load_Tri[127], Print[127], 
    removeFiles[127], command[127];
  bool approx;
  bool ehrhart_polynomial, ehrhart_series, ehrhart_taylor;
  bool triangulation_specified = false;
  double sampling_factor = 1.0;
  long int num_samples = -1;
  ReadPolyhedronData read_polyhedron_data;
  
  struct BarvinokParameters *params = new BarvinokParameters;

  latte_banner(cout);

  if (argc < 2) {
    usage(argv[0]);
    exit(1);
  }
  
  //setbuf(stdout,0);

  cout << "Invocation: ";
  for (i = 0; i<argc; i++) {
    cout << argv[i] << " ";
  }
  cout << endl;

  strcpy(removeFiles,"yes");
  strcpy(printfile,"no");
  strcpy(Save_Tri, "no");
  strcpy(Load_Tri, "no");
  strcpy(Print, "no");
  approx = false;
  ehrhart_polynomial = false;
  params->substitution = BarvinokParameters::PolynomialSubstitution;
  params->decomposition = BarvinokParameters::DualDecomposition;
  params->triangulation = BarvinokParameters::RegularTriangulationWithCdd;
  params->max_determinant = 1;

  for (i=1; i<argc; i++) {
    if(strncmp(argv[i],"nodecom",3)==0
	    || strncmp(argv[i], "--nodecomposition", 5) == 0
	    || strncmp(argv[i], "--no-decomposition", 7) == 0)
      params->max_determinant = 0;
    else if(strncmp(argv[i],"uni",3)==0) strcpy(read_polyhedron_data.assumeUnimodularCones,"yes");
    //else if(strncmp(argv[i],"simp",4)==0) {strcpy(printfile,"yes"); flags |= PRINT;}
    else if(strncmp(argv[i],"file",4)==0) strcpy(read_polyhedron_data.Memory_Save, "no");
    //else if(strncmp(argv[i],"single",6)==0) strcpy(Singlecone,"yes");
    //else if(strncmp(argv[i],"ehrhartsimp",3)==0) strcpy(rationalCone,"yes");
    else if(strncmp(argv[i],"memsave",7)==0) strcpy (read_polyhedron_data.Memory_Save, "yes");
    else if(strncmp(argv[i],"printcones",3)==0) strcpy (Print, "yes");
    //else if(strncmp(argv[i],"hull",3)==0) strcpy (inthull, "yes");
    else if(strncmp(argv[i],"rem",3)==0) {
      strcpy (removeFiles, "no");
      strcpy (read_polyhedron_data.Memory_Save, "no");
    }
    else if(strncmp(argv[i],"trisave",7)==0) {strcpy (Save_Tri, "yes"); flags |= SAVE;}
    else if(strncmp(argv[i],"triload",7)==0) {strcpy (Load_Tri, "yes"); flags |= LOAD;}
    else if (strncmp(argv[i], "--exponential", 5) == 0)
      params->substitution = BarvinokParameters::ExponentialSubstitution;
    else if (strncmp(argv[i], "--polynomial", 6) == 0)
      params->substitution = BarvinokParameters::PolynomialSubstitution;
    else if (strncmp(argv[i], "--maxdet=", 9) == 0)
      params->max_determinant = atoi(argv[i] + 9);
    else if (strncmp(argv[i], "--irrational-all-primal", 14) == 0
	     || strncmp(argv[i], "--all-primal", 5) == 0)
      params->decomposition = BarvinokParameters::IrrationalAllPrimalDecomposition;
    else if (strncmp(argv[i], "--irrational-primal", 5) == 0)
      params->decomposition = BarvinokParameters::IrrationalPrimalDecomposition;
    else if (strcmp(argv[i], "--dual") == 0) // Don't use strncmp to
					     // avoid clash with --dualization=...
      params->decomposition = BarvinokParameters::DualDecomposition;
    else if (strncmp(argv[i], "--count-lattice-points", 7) == 0) {
      // Default.
    }
    else if (strncmp(argv[i], "--ehrhart-polynomial", 11) == 0)
      ehrhart_polynomial = true;
    else if (strncmp(argv[i], "--ehrhart-series", 11) == 0) {
      ehrhart_series = true;
      strcpy(read_polyhedron_data.dualApproach,"yes");
      strcpy(printfile,"yes");
      flags |= PRINT;
    }
    else if (strncmp(argv[i], "--simplified-ehrhart-series", 14) == 0) {
      ehrhart_series = true;
      strcpy(read_polyhedron_data.dualApproach,"yes");
      strcpy(read_polyhedron_data.rationalCone, "yes");
    }
    else if (strncmp(argv[i], "--ehrhart-taylor=", 17) == 0) {
      strcpy(read_polyhedron_data.taylor, "yes");
      read_polyhedron_data.degree = atoi(argv[i] + 17);
      strcpy(read_polyhedron_data.dualApproach,"yes");
    }
    else if (strncmp(argv[i], "--avoid-singularities", 7) == 0) {
      params->shortvector = BarvinokParameters::SubspaceAvoidingLLL;
    }
    else if (parse_standard_triangulation_option(argv[i], params)) {
      if (strncmp(argv[i], "--triangulation=", 16) == 0)
	triangulation_specified = true;
    }
    else if (parse_standard_dualization_option(argv[i], params)) {}
    else if (strncmp(argv[i], "--approximate", 7) == 0)
      approx = true;
    else if (strncmp(argv[i], "--sampling-factor=", 18) == 0)
      sampling_factor = atof(argv[i] + 18);
    else if (strncmp(argv[i], "--num-samples=", 14) == 0)
      num_samples = atol(argv[i] + 14);
    else if (strncmp(argv[i], "--random-seed=", 14) == 0) {
      unsigned int seed = atoi(argv[i] + 14);
      seed_random_generator(seed);
    }
    else if (strcmp(argv[i], "--help") == 0) {
      read_polyhedron_data.show_options(cout);
      cout << "Options that control what to compute:" << endl
	   << "  --count-lattice-points                   Compute the number of lattice points" << endl
	   << "                                           (default)" << endl
	   << "  --ehrhart-polynomial                     Compute an Ehrhart polynomial of an integral polytope" << endl
	   << "  --ehrhart-series                         Compute the unsimplified Ehrhart series" << endl
	   << "                                           as a univariate rational function" << endl
	   << "  --simplified-ehrhart-series              Compute the simplified Ehrhart series" << endl
	   << "                                           as a univariate rational function" << endl
	   << "  --ehrhart-taylor=N                       Compute the first N terms of the Ehrhart series" << endl;
      cout << "Options for the Barvinok algorithm:" << endl
	   << "  --dual                                   Triangulate and signed-decompose in the dual space" << endl
	   << "                                           (traditional method, default)" << endl
	   << "  --irrational-primal                      Triangulate in the dual space, signed-decompose" << endl
	   << "                                           in the primal space using irrationalization" << endl
	   << "  --irrational-all-primal                  Triangulate and signed-decompose in the primal space" << endl
	   << "                                           using irrationalization" << endl
	   << "  --maxdet=N                               Decompose down to an index (determinant) of N" << endl
	   << "                                           instead of index 1 (unimodular cones)" << endl
	   << "  --no-decomposition                       Do not signed-decompose simplicial cones" << endl;
      cout << "Options for specialization:" << endl
	   << "  --polynomial                             Use polynomial substitution for specialization" << endl
	   << "                                           (traditional method, default)" << endl
	   << "  --exponential                            Use exponential substitution for specialization" << endl
	   << "                                           (recommended for maxdet > 1)" << endl;
      cout << "Algorithmic options for subproblems:" << endl;
      show_standard_dualization_option(cout);
      show_standard_triangulation_options(cout);
      exit(0);
    }
    else if (read_polyhedron_data.parse_option(argv[i])) {}
    else {
      cerr << "Unknown command/option " << argv[i] << endl;
      exit(1);
    }
  }

  if (read_polyhedron_data.expect_filename) {
    cerr << "Filename missing" << endl;
    exit(1);
  }
  
  if (params->shortvector == BarvinokParameters::SubspaceAvoidingLLL) {
    if (params->decomposition == BarvinokParameters::IrrationalAllPrimalDecomposition) {
      /* Triangulation will be done in the primal space, so all
	 triangulation methods are fine. */
    }
    else {
      /* Triangulation will be done in the dual space, so we must
	 avoid using facets whose normal vectors lie in the
	 subspace. */
      if (triangulation_specified && params->triangulation != BarvinokParameters::SubspaceAvoidingBoundaryTriangulation) {
	cerr << "Warning: The requested triangulation method is not guaranteed to work with --avoid-singularities."
	     << endl;
      }
      else {
	params->triangulation = BarvinokParameters::SubspaceAvoidingBoundaryTriangulation;
      }
    }
  }

  if (approx) {
    params->substitution = BarvinokParameters::ExponentialSubstitution;
    if (params->decomposition == BarvinokParameters::DualDecomposition) {
      cout << "Exponential approximation not implemented for dual decomposition; switching to irrational primal decomposition." << endl;
      params->decomposition = BarvinokParameters::IrrationalPrimalDecomposition;
    }
  }
  
  if(read_polyhedron_data.minimize[0] == 'y') strcpy(read_polyhedron_data.maximum, "yes");
  if(read_polyhedron_data.grobner[0] == 'y') strcpy(read_polyhedron_data.equationsPresent,"yes");
  if(read_polyhedron_data.maximum[0] == 'y') strcpy(read_polyhedron_data.Memory_Save, "no");
  if(printfile[0] == 'y') strcpy(read_polyhedron_data.Memory_Save, "no");
  if(read_polyhedron_data.rationalCone[0] == 'y') strcpy(read_polyhedron_data.Memory_Save, "no");
  if(printfile[0] == 'y') print_flag = 1;

  if(read_polyhedron_data.rationalCone[0] == 'y'){
    
    //HugInt digit(argv[1]);
    //conv(output_cone, digit.BigInt);
    // User can use only Mode one
    output_cone = 3;
  }

  if(output_cone > 3) output_cone = 0;
  flags |= (output_cone << 1);

  const char *fileName = read_polyhedron_data.filename.c_str();

  if (read_polyhedron_data.dualApproach[0] == 'y') {
    flags |= DUAL_APPROACH;
  }
  
  /* INPUT HANDLING. */

  if(read_polyhedron_data.grobner[0] == 'y'){
    CheckGrobner(fileName, read_polyhedron_data.cddstyle);
    SolveGrobner(fileName,  read_polyhedron_data.nonneg, read_polyhedron_data.dualApproach,
		 read_polyhedron_data.grobner, read_polyhedron_data.equationsPresent, read_polyhedron_data.cddstyle);
    exit(0);
  }

  Polyhedron *Poly = read_polyhedron_data.read_polyhedron(params);

  if (ehrhart_polynomial) {
    /* Translate all cones to the origin, saving the original vertex. */
    listCone *cone;
    for (cone = Poly->cones; cone; cone = cone->rest) {
      ZZ scale_factor;
      cone->vertex->ehrhart_vertex
	= scaleRationalVectorToInteger(cone->vertex->vertex,
				       Poly->numOfVars, scale_factor);
      if (scale_factor != 1) {
	cerr << "Computation of Ehrhart polynomials is only implemented "
	     << "for integral polytopes." << endl
	     << "Use `ehrhart' for computing the Ehrhart series "
	     << "of rational polytopes." << endl;
	exit(1);
      }
      delete cone->vertex->vertex;
      cone->vertex->vertex = new rationalVector(Poly->numOfVars);
    }
  }
    
  params->Flags = flags;
  params->File_Name = (char*) fileName;
  params->Number_of_Variables = Poly->numOfVars;

  switch (params->decomposition) {
  case BarvinokParameters::DualDecomposition:
  case BarvinokParameters::IrrationalPrimalDecomposition:
    if (not Poly->dualized) {
      if (read_polyhedron_data.Vrepresentation[0] != 'y') {
	/* Compute all inequalities tight at the respective vertex.
	   Then dualizeCones just needs to swap rays and facets. */
	computeTightInequalitiesOfCones(Poly->cones, read_polyhedron_data.matrix, Poly->numOfVars);
      }
      dualizeCones(Poly->cones, Poly->numOfVars, params);
      Poly->dualized = true;
    }
    break;
  case BarvinokParameters::IrrationalAllPrimalDecomposition:
    cout << "Irrationalizing polyhedral cones... "; cout.flush();
    if (Poly->dualized) {
      cout << "(First dualizing back... "; cout.flush();
      dualizeCones(Poly->cones, Poly->numOfVars, params);
      cout << "done; sorry for the interruption.) "; cout.flush();
    }
    else {
      if (read_polyhedron_data.Vrepresentation[0] == 'y') {
	cout << "(First computing facets for them... "; cout.flush();
	dualizeCones(Poly->cones, Poly->numOfVars, params);
	dualizeCones(Poly->cones, Poly->numOfVars, params); // just swaps
	cout << "done; sorry for the interruption.) "; cout.flush();
      }      
      else {
	/* Fill in the facets of all cones; we determine them by
	   taking all inequalities tight at the respective vertex. */
	params->dualize_time.start();
	computeTightInequalitiesOfCones(Poly->cones, read_polyhedron_data.matrix, Poly->numOfVars);
	params->dualize_time.stop(); cout << params->dualize_time;
      }
      if (Poly->cones && Poly->cones->rays == NULL) {
	/* Only facets computed, for instance by using the 4ti2
	   method of computing vertex cones.  So dualize twice to
	   compute the rays. */
	cout << "(First computing their rays... "; cout.flush();
	dualizeCones(Poly->cones, Poly->numOfVars, params);
	dualizeCones(Poly->cones, Poly->numOfVars, params); // just swaps
	cout << "done; sorry for the interruption.) "; cout.flush();
      }
    }
    params->irrationalize_time.start();
    {
      listCone *cone;
      for (cone = Poly->cones; cone; cone=cone->rest)
	assert(lengthListVector(cone->facets) >= Poly->numOfVars);
    }
    irrationalizeCones(Poly->cones, Poly->numOfVars);
    params->irrationalize_time.stop();
    cout << params->irrationalize_time;
    break;
  default:
    cerr << "Unknown BarvinokParameters::decomposition" << endl;
    abort();
  }

  try {
    
  switch (params->substitution) {
  case BarvinokParameters::PolynomialSubstitution:
    if (ehrhart_polynomial) {
      cerr << "Computation of Ehrhart polynomials is only implemented "
	   << "for the exponential substitution." << endl;
      exit(1);
    }
    if (read_polyhedron_data.assumeUnimodularCones[0]=='n') {
      if (read_polyhedron_data.Memory_Save[0] == 'n') {
	listCone *decomposed_cones
	  = decomposeCones(Poly->cones, not Poly->dualized,
			   *params);
	freeListCone(Poly->cones);
 	Poly->cones = decomposed_cones;
// 	cout << "Decomposed cones: " << endl;
//	printListCone(Poly->cones, Poly->numOfVars);
	/* Compute points in parallelepipeds */
	computePointsInParallelepipeds(Poly->cones, Poly->numOfVars);
      }
      // Iterator through simplicial cones, DFS
      else {
	Standard_Single_Cone_Parameters *standard_params
	  = new Standard_Single_Cone_Parameters(*params);
	delete params; params = standard_params;
	decomposeAndComputeResidue(Poly->cones, read_polyhedron_data.degree, false,
				   *standard_params);
      }
    }
    break;
  case BarvinokParameters::ExponentialSubstitution:
    if (read_polyhedron_data.dualApproach[0] == 'y') {
      cerr << "Exponential substitution is not yet implemented for the homogenized version."
	   << endl;
      exit(1);
    }
    else {
      if (approx) {
#ifdef HAVE_EXPERIMENTS
	{
	  Write_Exponential_Sample_Formula_Single_Cone_Parameters *write_param
	    = new Write_Exponential_Sample_Formula_Single_Cone_Parameters
	    (*params, "Exponential_Sample_Formula", sampling_factor,
	     num_samples);
	  delete params;
	  params = write_param;
	  decomposeAndWriteExponentialSampleFormula(Poly->cones, *write_param);
	}
#else
	cerr << "Approximation code is not compiled in, sorry." << endl;
	exit(1);
#endif
      }
      else if (ehrhart_polynomial) {
	Exponential_Ehrhart_Parameters *exp_param
	  = new Exponential_Ehrhart_Parameters(*params);
	delete params;
	params = exp_param;
	mpq_vector ehrhart_coefficients
	  = decomposeAndComputeEhrhartPolynomial(Poly->cones, *exp_param);
	cout << endl << "Ehrhart polynomial: ";
	{
	  unsigned int i;
	  for (i = 0; i<ehrhart_coefficients.size(); i++) {
	    if (ehrhart_coefficients[i] > 0) 
	      cout << " + " << ehrhart_coefficients[i] << " * t^" << i;
	    else if (ehrhart_coefficients[i] < 0)
	      cout << " - " << abs(ehrhart_coefficients[i]) << " * t^" << i;
	  }
	}
	cout << endl << endl;
      }
      else {
	Exponential_Single_Cone_Parameters *exp_param
	  = new Exponential_Single_Cone_Parameters(*params);
	delete params;
	params = exp_param;
	Integer number_of_lattice_points
	  = decomposeAndComputeExponentialResidue(Poly->cones, *exp_param);
	cout << endl << "****  The number of lattice points is: "
	     << number_of_lattice_points << "  ****" << endl << endl;
	// FIXME: Centralize this output stuff.
	ofstream out("numOfLatticePoints");
	out << number_of_lattice_points << endl;
      }
    }
    break;
  default:
    cerr << "Unknown BarvinokParameters::substitution" << endl;
    abort();
  }

  if(read_polyhedron_data.grobner[0] == 'y'){

 Poly->cones = ProjectUp(Poly->cones, read_polyhedron_data.oldnumofvars,
			 Poly->numOfVars, read_polyhedron_data.templistVec);
 Poly->numOfVars = read_polyhedron_data.oldnumofvars;

  }
 if(Print[0] == 'y')
  printListCone(Poly->cones,Poly->numOfVars);

 //   printListVector(IntegralHull(Poly->cones,  inequalities, equations, Poly->numOfVars), Poly->numOfVars);
 if(read_polyhedron_data.Memory_Save[0] == 'n')
   {

	if(read_polyhedron_data.dualApproach[0] == 'n'){
	  cout << "Creating generating function.\n"; 
	  //printListVector(templistVec, oldnumofvars); cout << ProjU << endl;
	if(read_polyhedron_data.equationsPresent[0] == 'y'){
	  Poly->cones = ProjectUp2(Poly->cones, read_polyhedron_data.oldnumofvars,
				   Poly->numOfVars, read_polyhedron_data.AA, read_polyhedron_data.bb);
	  Poly->numOfVars = read_polyhedron_data.oldnumofvars;
	}

	  createGeneratingFunctionAsMapleInput(fileName,Poly->cones,Poly->numOfVars);  }
        //printListCone(cones, Poly->numOfVars);

	cout << "Printing decomposed cones to `decomposed_cones'." << endl;
	printListConeToFile("decomposed_cones", Poly->cones, Poly->numOfVars);

	if(read_polyhedron_data.dualApproach[0] == 'n'){
	cout << "Starting final computation.\n";
	cout << endl << "****  The number of lattice points is: " << Residue(Poly->cones,Poly->numOfVars) << "  ****" << endl << endl;}


	if(read_polyhedron_data.dualApproach[0] == 'y')
	{
	  cout << "Starting final computation.\n";
	  //cout << "output_cone: " << output_cone;
	  switch (params->decomposition) {
	  case BarvinokParameters::IrrationalPrimalDecomposition:
	  case BarvinokParameters::IrrationalAllPrimalDecomposition: {
#ifdef HAVE_EXPERIMENTS
	    ofstream out("func.rat");
	    out << "HS := ";
	    TrivialMonomialSubstitutionMapleOutput(out, Poly->cones, Poly->numOfVars);
	    out << ";";
#else
	    cerr << "Trivial monomial substitution not compiled in, sorry." << endl;
#endif
	    break;
	  }
	  case BarvinokParameters::DualDecomposition:
	    ResidueFunction(Poly->cones,Poly->numOfVars, print_flag, read_polyhedron_data.degree, output_cone);
	    // ResidueFunction consumes cones.
	    Poly->cones = NULL;
	    break;
	  default:
	    assert(0);
	  }
	//  Else we have already computed the residue.
	}
   }
  } catch (NotIrrationalException) {
    cerr << "Bug: Irrationalization failed" << endl;
    exit(1);
  };

  freeListVector(read_polyhedron_data.templistVec);
  freeListVector(read_polyhedron_data.matrix);
  delete Poly;

if(read_polyhedron_data.rationalCone[0] == 'y') {
   cout << endl <<"Rational function written to " << argv[argc - 1] << ".rat" << endl << endl;
   strcpy(command, "mv ");
   strcat(command, "simplify.sum ");
   strcat(command, argv[argc - 1]);
   strcat(command, ".rat");
   system_with_error_check(command);
 }

 if(printfile[0] == 'y'){
   cout << endl <<"Rational function written to " << argv[argc - 1] << ".rat" << endl << endl;
   strcpy(command, "mv ");
   strcat(command, "func.rat ");
   strcat(command, argv[argc - 1]);
   strcat(command, ".rat");
   system_with_error_check(command);
 }
 if((removeFiles[0] == 'y') && (read_polyhedron_data.dualApproach[0] == 'n')){
   
  strcpy(command,"rm -f ");
  strcat(command,fileName);
  strcat(command,".ext");
  system_with_error_check(command);
  
  strcpy(command,"rm -f ");
  strcat(command,fileName);
  strcat(command,".cdd");
  system_with_error_check(command); 
  
  if(read_polyhedron_data.Memory_Save[0] == 'n'){
    strcpy(command,"rm -f ");
    strcat(command,fileName);
    strcat(command,".maple");
    system_with_error_check(command); 
  }

  strcpy(command,"rm -f ");
  strcat(command,fileName);
  strcat(command,".ead");
  system_with_error_check(command); 
  
 }

  //cout << "Computation done. " << endl;

  params->total_time.stop();
  cout << params->total_time;
  
  {
    // until we have a more sophisticated test script --mkoeppe
    ofstream totalTime("totalTime");
    totalTime << params->total_time.get_seconds()
	      << " (" << params->read_time.get_seconds() << "r"
	      << ", " << params->vertices_time.get_seconds() << "v"
	      << ", " << params->irrationalize_time.get_seconds() << "i"
	      << ", " << params->dualize_time.get_seconds() << "d"
	      << ", " << params->triangulate_time.get_seconds() << "t"
	      << ", " << params->decompose_time.get_seconds() << "b"
	      << ")" << endl;
    ofstream stats("latte_stats");
    params->print_statistics(stats);
  }

  delete params;
 return(0);
}
/* ----------------------------------------------------------------- */



