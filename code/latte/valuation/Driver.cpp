/*
 * Driver.cpp
 *
 *  Created on: Jun 24, 2010
 *      Author: Brandon Dutra and Gregory Pinto
 *
 *
 *
 *  type --help to print the help menu.
 */




#include <cstdlib>
#include <iostream>
#include <string>
#include <cassert>
#include <vector>







#include "barvinok/barvinok.h"
#include "ReadPolyhedron.h"
#include "triangulation/triangulate.h"
#include "convert.h"
#include "print.h"
#include "gnulib/progname.h"
#include "barvinok/dec.h"
#include "valuation/PolytopeValuation.h"
#include "barvinok/dec.h"
#include "rational.h"
#include <cstdlib>
#include <ctime>
#include "cone.h"
#include <NTL/vec_ZZ.h>
#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <NTL/mat_ZZ.h>
#include "testEhrhart/BuildRandomPolytope.h"
#include "testEhrhart/BuildHypersimplexEdgePolytope.h"

/* EHRHART INCLUDES */
#include <string.h>
#include <stdio.h>

#include "config.h"
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
#include "preprocess.h"
#include "print.h"
#include "ramon.h"
#include "rational.h"
#include "timing.h"
#include "flags.h"
#include "IntegralHull.h"
#include "ReadingFile.h"
#include "binarySearchIP.h"
#include "CheckEmpty.h"
#include "ProjectUp.h"

#include "banner.h"
#include "convert.h"
#include "latte_system.h"

#include "gnulib/progname.h"
/* END EHRHART INCLUDES */


/* START COUNT INCLUDES */

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
#include "preprocess.h"
#include "print.h"
#include "ramon.h"
#include "rational.h"
#include "timing.h"
#include "flags.h"
#include "ExponentialSubst.h"
#include "latte_random.h"
#include "Irrational.h"
#include "ExponentialEhrhart.h"
#include "triangulation/triangulate.h"
#include "genFunction/matrix_ops.h"
#ifdef HAVE_EXPERIMENTS
#include "ExponentialApprox.h"
#include "TrivialSubst.h"
#endif

#include "banner.h"
#include "convert.h"
#include "latte_system.h"
#include "Polyhedron.h"
#include "ReadPolyhedron.h"
#include "ProjectUp.h"

#include "gnulib/progname.h"

/* END COUNT INCLUDES */
using namespace std;




typedef struct
{
	RationalNTL triangulate;
	RationalNTL lawrence;

} VolumesContainer;


VolumesContainer computeVolume(listCone * cones, BarvinokParameters &myParameters,
		const char *valuationType);
VolumesContainer mainValuationDriver(const char *argv[], int argc);
void runOneTest(int ambientDim, int numPoints);
void runTests();

void printVolumeTest(const RationalNTL &correctVolumeAnswer, VolumesContainer volumeAnswer, const string &file, const string &comments);

//BarvinokParameters parameters;
//ReadPolyhedronData read_polyhedron_data;
//string output_filename;


/* ----------------------------------------------------------------- */

/**
 * Computes the volume of the polytope.
 */
VolumesContainer computeVolume(listCone * cones, BarvinokParameters &myParameters,
		const char *valuationAlg, const char * print)
{
	VolumesContainer ans;
	RationalNTL ans1, ans2;
	if (strcmp(valuationAlg, "triangulate") == 0 || strcmp(valuationAlg, "all")
			== 0)
	{
		PolytopeValuation polytopeValuation(cones,
				PolytopeValuation::VertexRayCones,
				myParameters.Number_of_Variables, myParameters);
		ans1 = polytopeValuation.findVolume(PolytopeValuation::DeterminantVolume);

		cout << "(using triangulation-determinant method) VOLUME " << ans1 << endl;
	}//if triangulate

	if (strncmp(valuationAlg, "lawrence", 8) == 0
			|| strcmp(valuationAlg, "all") == 0)
	{

		PolytopeValuation polytopeValuation(cones,
				PolytopeValuation::VertexRayCones,
				myParameters.Number_of_Variables, myParameters);
		ans2 = polytopeValuation.findVolume(PolytopeValuation::LawrenceVolume);

		cout << "(using Lawrence method) VOLUME " << ans2 << endl;

		if (*print == 'y')
			polytopeValuation.printLawrenceVolumeFunction(); //print the lawrence rational function.

	}//if lawrence

	if ( strcmp(valuationAlg, "all") == 0 && ans1 != ans2)
	{
		cout << "Driver.cpp: the two methods are different." << ans1 << "!=" << ans2 << endl;
		exit(1);
	}//if error.


	ans.lawrence = ans2;
	ans.triangulate = ans1;

	return ans;

}//computeVolume


static void usage(const char *progname)
{
	cerr << "usage: " << progname << " [OPTIONS...] " << "INPUTFILE" << endl;
	cerr << "Type `" << progname << " --help' "
			<< "for a list of options and input specifications." << endl;
}


VolumesContainer mainValuationDriver(const char *argv[], int argc)
{
	VolumesContainer volumeAnswers;
	set_program_name(argv[0]);

	int i;
	unsigned int flags = 0, print_flag = 0, output_cone = 0;
	char printfile[127], Save_Tri[127], Load_Tri[127], Print[127],
			removeFiles[127], command[127];
	bool approx;
	bool ehrhart_polynomial, ehrhart_series, ehrhart_taylor;
	bool triangulation_specified = false;
	double sampling_factor = 1.0;
	long int num_samples = -1;
	ReadPolyhedronData read_polyhedron_data;

	struct BarvinokParameters *params = new BarvinokParameters;

	latte_banner(cerr);

	if (argc < 2)
	{
		usage(argv[0]);
		exit(1);
	}

	//setbuf(stdout,0);

	cerr << "Invocation: ";
	for (i = 0; i < argc; i++)
	{
		cerr << argv[i] << " ";
	}
	cerr << endl;

	strcpy(removeFiles, "yes");
	strcpy(printfile, "no");
	strcpy(Save_Tri, "no");
	strcpy(Load_Tri, "no");
	strcpy(Print, "yes");
	approx = false;
	ehrhart_polynomial = false;
	params->substitution = BarvinokParameters::PolynomialSubstitution;
	params->decomposition = BarvinokParameters::DualDecomposition;
	params->triangulation = BarvinokParameters::RegularTriangulationWithCdd;
	params->max_determinant = 1;

	for (i = 1; i < argc; i++)
	{
		if (strncmp(argv[i], "nodecom", 3) == 0 || strncmp(argv[i],
				"--nodecomposition", 5) == 0 || strncmp(argv[i],
				"--no-decomposition", 7) == 0)
			params->max_determinant = 0;
		else if (strncmp(argv[i], "uni", 3) == 0)
			strcpy(read_polyhedron_data.assumeUnimodularCones, "yes");
		//else if(strncmp(argv[i],"simp",4)==0) {strcpy(printfile,"yes"); flags |= PRINT;}
		else if (strncmp(argv[i], "file", 4) == 0)
			strcpy(read_polyhedron_data.Memory_Save, "no");
		//else if(strncmp(argv[i],"single",6)==0) strcpy(Singlecone,"yes");
		//else if(strncmp(argv[i],"ehrhartsimp",3)==0) strcpy(rationalCone,"yes");
		else if (strncmp(argv[i], "memsave", 7) == 0)
			strcpy(read_polyhedron_data.Memory_Save, "yes");
		else if (strncmp(argv[i], "printcones", 3) == 0)
			strcpy(Print, "yes");
		//else if(strncmp(argv[i],"hull",3)==0) strcpy (inthull, "yes");
		else if (strncmp(argv[i], "rem", 3) == 0)
		{
			strcpy(removeFiles, "no");
			strcpy(read_polyhedron_data.Memory_Save, "no");
		} else if (strncmp(argv[i], "trisave", 7) == 0)
		{
			strcpy(Save_Tri, "yes");
			flags |= SAVE;
		} else if (strncmp(argv[i], "triload", 7) == 0)
		{
			strcpy(Load_Tri, "yes");
			flags |= LOAD;
		} else if (strncmp(argv[i], "--exponential", 5) == 0)
			params->substitution = BarvinokParameters::ExponentialSubstitution;
		else if (strncmp(argv[i], "--polynomial", 6) == 0)
			params->substitution = BarvinokParameters::PolynomialSubstitution;
		else if (strncmp(argv[i], "--maxdet=", 9) == 0)
			params->max_determinant = atoi(argv[i] + 9);
		else if (strncmp(argv[i], "--irrational-all-primal", 14) == 0
				|| strncmp(argv[i], "--all-primal", 5) == 0)
			params->decomposition
					= BarvinokParameters::IrrationalAllPrimalDecomposition;
		else if (strncmp(argv[i], "--irrational-primal", 5) == 0)
			params->decomposition
					= BarvinokParameters::IrrationalPrimalDecomposition;
		else if (strcmp(argv[i], "--dual") == 0) // Don't use strncmp to
			// avoid clash with --dualization=...
			params->decomposition = BarvinokParameters::DualDecomposition;
		else if (strncmp(argv[i], "--count-lattice-points", 7) == 0)
		{
			// Default.
		} else if (strncmp(argv[i], "--multivariate-generating-function", 7)
				== 0)
		{
			params->substitution = BarvinokParameters::NoSubstitution;
		} else if (strncmp(argv[i], "--ehrhart-polynomial", 11) == 0)
		{
			ehrhart_polynomial = true;
			params->substitution = BarvinokParameters::ExponentialSubstitution;
		} else if (strncmp(argv[i], "--ehrhart-series", 11) == 0)
		{
			ehrhart_series = true;
			strcpy(read_polyhedron_data.dualApproach, "yes");
			strcpy(printfile, "yes");
			flags |= PRINT;
		} else if (strncmp(argv[i], "--simplified-ehrhart-series", 14) == 0)
		{
			ehrhart_series = true;
			strcpy(read_polyhedron_data.dualApproach, "yes");
			strcpy(read_polyhedron_data.rationalCone, "yes");
		} else if (strncmp(argv[i], "--ehrhart-taylor=", 17) == 0)
		{
			strcpy(read_polyhedron_data.taylor, "yes");
			read_polyhedron_data.degree = atoi(argv[i] + 17);
			strcpy(read_polyhedron_data.dualApproach, "yes");
		} else if (strncmp(argv[i], "--avoid-singularities", 7) == 0)
		{
			params->shortvector = BarvinokParameters::SubspaceAvoidingLLL;
		} else if (parse_standard_triangulation_option(argv[i], params))
		{
			if (strncmp(argv[i], "--triangulation=", 16) == 0)
				triangulation_specified = true;
		} else if (parse_standard_dualization_option(argv[i], params))
		{
		} else if (parse_standard_smith_option(argv[i], params))
		{
		} else if (strncmp(argv[i], "--approximate", 7) == 0)
			approx = true;
		else if (strncmp(argv[i], "--sampling-factor=", 18) == 0)
			sampling_factor = atof(argv[i] + 18);
		else if (strncmp(argv[i], "--num-samples=", 14) == 0)
			num_samples = atol(argv[i] + 14);
		else if (strncmp(argv[i], "--random-seed=", 14) == 0)
		{
			unsigned int seed = atoi(argv[i] + 14);
			seed_random_generator(seed);
		} else if (strcmp(argv[i], "--help") == 0)
		{
			read_polyhedron_data.show_options(cerr);
			cerr << "Options that control what to compute:" << endl
					<< "  --count-lattice-points                   Compute the number of lattice points"
					<< endl
					<< "                                           (default)"
					<< endl
					<< "  --multivariate-generating-function       Compute the multivariate generating function of"
					<< endl
					<< "                                           the set of lattice points of the polyhedron"
					<< endl
					<< "  --ehrhart-polynomial                     Compute an Ehrhart polynomial of an integral polytope"
					<< endl
					<< "  --ehrhart-series                         Compute the unsimplified Ehrhart series"
					<< endl
					<< "                                           as a univariate rational function"
					<< endl
					<< "  --simplified-ehrhart-series              Compute the simplified Ehrhart series"
					<< endl
					<< "                                           as a univariate rational function"
					<< endl
					<< "  --ehrhart-taylor=N                       Compute the first N terms of the Ehrhart series"
					<< endl;
			cerr << "Options for the Barvinok algorithm:" << endl
					<< "  --dual                                   Triangulate and signed-decompose in the dual space"
					<< endl
					<< "                                           (traditional method, default)"
					<< endl
					<< "  --irrational-primal                      Triangulate in the dual space, signed-decompose"
					<< endl
					<< "                                           in the primal space using irrationalization"
					<< endl
					<< "  --irrational-all-primal                  Triangulate and signed-decompose in the primal space"
					<< endl
					<< "                                           using irrationalization"
					<< endl
					<< "  --maxdet=N                               Decompose down to an index (determinant) of N"
					<< endl
					<< "                                           instead of index 1 (unimodular cones)"
					<< endl
					<< "  --no-decomposition                       Do not signed-decompose simplicial cones"
					<< endl;
			cerr << "Options for specialization:" << endl
					<< "  --polynomial                             Use polynomial substitution for specialization"
					<< endl
					<< "                                           (traditional method, default)"
					<< endl
					<< "  --exponential                            Use exponential substitution for specialization"
					<< endl
					<< "                                           (recommended for maxdet > 1)"
					<< endl;
			cerr << "Algorithmic options for subproblems:" << endl;
			show_standard_smith_option(cerr);
			show_standard_dualization_option(cerr);
			show_standard_triangulation_options(cerr);
			exit(0);
		} else if (read_polyhedron_data.parse_option(argv[i]))
		{
		} else
		{
			cerr << "Unknown command/option " << argv[i] << endl;
			exit(1);
		}
	}

	if (read_polyhedron_data.expect_filename)
	{
		cerr << "Filename missing" << endl;
		exit(1);
	}

	if (params->shortvector == BarvinokParameters::SubspaceAvoidingLLL)
	{
		if (params->decomposition
				== BarvinokParameters::IrrationalAllPrimalDecomposition)
		{
			/* Triangulation will be done in the primal space, so all
			 triangulation methods are fine. */
		} else
		{
			/* Triangulation will be done in the dual space, so we must
			 avoid using facets whose normal vectors lie in the
			 subspace. */
			if (triangulation_specified)
			{
				if (params->triangulation
						!= BarvinokParameters::SubspaceAvoidingBoundaryTriangulation
						&& params->triangulation
								!= BarvinokParameters::SubspaceAvoidingSpecialTriangulation)
				{
					cerr
							<< "Warning: The requested triangulation method is not guaranteed to work with --avoid-singularities."
							<< endl;
				}
			} else
			{
				/* Not specified, so choose one that will work. */
				cerr << "Choosing SubspaceAvoidingBoundaryTriangulation method"
						<< endl;
				params->triangulation
						= BarvinokParameters::SubspaceAvoidingBoundaryTriangulation;
			}
		}
	}

	if (approx)
	{
		params->substitution = BarvinokParameters::ExponentialSubstitution;
		if (params->decomposition == BarvinokParameters::DualDecomposition)
		{
			cerr
					<< "Exponential approximation not implemented for dual decomposition; switching to irrational primal decomposition."
					<< endl;
			params->decomposition
					= BarvinokParameters::IrrationalPrimalDecomposition;
		}
	}

	if (read_polyhedron_data.minimize[0] == 'y')
		strcpy(read_polyhedron_data.maximum, "yes");
	if (read_polyhedron_data.grobner[0] == 'y')
		strcpy(read_polyhedron_data.equationsPresent, "yes");
	if (read_polyhedron_data.maximum[0] == 'y')
		strcpy(read_polyhedron_data.Memory_Save, "no");
	if (printfile[0] == 'y')
		strcpy(read_polyhedron_data.Memory_Save, "no");
	if (read_polyhedron_data.rationalCone[0] == 'y')
		strcpy(read_polyhedron_data.Memory_Save, "no");
	if (printfile[0] == 'y')
		print_flag = 1;

	if (read_polyhedron_data.rationalCone[0] == 'y')
	{

		//HugInt digit(argv[1]);
		//conv(output_cone, digit.BigInt);
		// User can use only Mode one
		output_cone = 3;
	}

	if (output_cone > 3)
		output_cone = 0;
	flags |= (output_cone << 1);

	const char *fileName = read_polyhedron_data.filename.c_str();

	if (read_polyhedron_data.dualApproach[0] == 'y')
	{
		flags |= DUAL_APPROACH;
	}

	/* INPUT HANDLING. */

	if (read_polyhedron_data.grobner[0] == 'y')
	{
		CheckGrobner(fileName, read_polyhedron_data.cddstyle);
		SolveGrobner(fileName, read_polyhedron_data.nonneg,
				read_polyhedron_data.dualApproach,
				read_polyhedron_data.grobner,
				read_polyhedron_data.equationsPresent,
				read_polyhedron_data.cddstyle);
		exit(0);
	}

	Polyhedron *Poly = read_polyhedron_data.read_polyhedron(params);

	//VolumesContainer computeVolume(listCone * cones, BarvinokParameters &myParameters,
	//		const char *valuationAlg, const char * print)

	params->Flags = flags;
	params->Number_of_Variables = Poly->numOfVars;
	params->max_determinant = 1;
	//params->File_Name = fileName;
	params->decomposition = BarvinokParameters::DualDecomposition;


	ofstream file;
	file.open("cones_count.txt");
	for(listCone * t = Poly->cones; t; t = t->rest)
		printConeToFile(file, t, Poly->numOfVars);
	file.close();


	PolytopeValuation polytopeValuation(Poly->cones,
					PolytopeValuation::VertexRayCones,
					Poly->numOfVars, *params);

	RationalNTL ans = polytopeValuation.findVolume(PolytopeValuation::DeterminantVolume);
	cout << "VOLUME BY TRIANGULATION: " << ans << endl;



	PolytopeValuation polytopeValuation2(Poly->cones,
					PolytopeValuation::VertexRayCones,
					Poly->numOfVars, *params);
	RationalNTL ans2 = polytopeValuation2.findVolume(PolytopeValuation::LawrenceVolume);

	cout << "VOLUME BY LAWRENCE " << ans2 << endl;



	//computeVolume(Poly->cones, *params, "all", "no");
	exit(1);

	if (Print[0] == 'y')
	{
		cout << "printing on line 512" << endl;
		printListCone(Poly->cones, Poly->numOfVars);
		cout << "end printing on line 512" << endl;
	}


	/* Compute the facet information from tight inequalities if
	 possible.  It is essential that this is done BEFORE translating
	 vertexes to the origin (in Ehrhart polynomial mode) -- otherwise
	 tightness information is wrong. */
	if (not Poly->dualized && Poly->cones != NULL
			&& read_polyhedron_data.matrix != NULL
			&& read_polyhedron_data.Vrepresentation[0] != 'y')
	{
		/* Fill in the facets of all cones; we determine them by
		 taking all inequalities tight at the respective vertex. */
		params->dualize_time.start();
		computeTightInequalitiesOfCones(Poly->cones,
				read_polyhedron_data.matrix, Poly->numOfVars);
		params->dualize_time.stop();
		cerr << params->dualize_time;
	}

	if (ehrhart_polynomial)
	{
		/* Translate all cones to the origin, saving the original vertex. */
		listCone *cone;
		for (cone = Poly->cones; cone; cone = cone->rest)
		{
			ZZ scale_factor;
			cone->vertex->ehrhart_vertex = scaleRationalVectorToInteger(
					cone->vertex->vertex, Poly->numOfVars, scale_factor);
			if (scale_factor != 1)
			{
				cerr
						<< "Computation of Ehrhart (quasi-)polynomials is only implemented "
						<< "for integral polytopes." << endl
						<< "Use `--ehrhart-series' or `--simplfied-ehrhart-series' for computing "
						<< "the Ehrhart series of rational polytopes." << endl;
				exit(1);
			}
			delete cone->vertex->vertex;
			cone->vertex->vertex = new rationalVector(Poly->numOfVars);
		}
	}

	params->Flags = flags;
	params->File_Name = (char*) fileName;
	params->Number_of_Variables = Poly->numOfVars;

	if (Print[0] == 'y')
	{
		cout << "printing on line 560" << endl;
		printListCone(Poly->cones, Poly->numOfVars);
		cout << "end printing on line 560" << endl;
	}
	cout << "calling exit on line 569" << endl;
	exit(1);



	switch (params->decomposition)
	{
		case BarvinokParameters::DualDecomposition:
		case BarvinokParameters::IrrationalPrimalDecomposition:
			if (not Poly->dualized)
			{
				if (read_polyhedron_data.Vrepresentation[0] != 'y')
				{
					/* Compute all inequalities tight at the respective vertex.
					 Then dualizeCones just needs to swap rays and facets. */
					computeTightInequalitiesOfCones(Poly->cones,
							read_polyhedron_data.matrix, Poly->numOfVars);
				}
				dualizeCones(Poly->cones, Poly->numOfVars, params);
				Poly->dualized = true;
			}
			break;
		case BarvinokParameters::IrrationalAllPrimalDecomposition:
			cerr << "Irrationalizing polyhedral cones... ";
			cerr.flush();
			if (Poly->dualized)
			{
				cerr << "(First dualizing back... ";
				cerr.flush();
				dualizeCones(Poly->cones, Poly->numOfVars, params);
				cerr << "done; sorry for the interruption.) ";
				cerr.flush();
			} else if (Poly->cones != NULL)
			{
				if (Poly->cones->facets == NULL)
				{
					cerr << "(First computing facets for them... ";
					cerr.flush();
					dualizeCones(Poly->cones, Poly->numOfVars, params);
					dualizeCones(Poly->cones, Poly->numOfVars, params); // just swaps
					cerr << "done; sorry for the interruption.) ";
					cerr.flush();
				} else if (Poly->cones->rays == NULL)
				{
					/* Only facets computed, for instance by using the 4ti2
					 method of computing vertex cones.  So dualize twice to
					 compute the rays. */
					cerr << "(First computing their rays... ";
					cerr.flush();
					dualizeCones(Poly->cones, Poly->numOfVars, params);
					dualizeCones(Poly->cones, Poly->numOfVars, params); // just swaps
					cerr << "done; sorry for the interruption.) ";
					cerr.flush();
				}
			}
			params->irrationalize_time.start();
			{
				listCone *cone;
				for (cone = Poly->cones; cone; cone = cone->rest)
					assert(lengthListVector(cone->facets) >= Poly->numOfVars);
			}
			irrationalizeCones(Poly->cones, Poly->numOfVars);
			params->irrationalize_time.stop();
			cerr << params->irrationalize_time;
			break;
		default:
			cerr << "Unknown BarvinokParameters::decomposition" << endl;
			abort();
	}

	if (Print[0] == 'y')
	{
		cout << "printing on line 630" << endl;
		printListCone(Poly->cones, Poly->numOfVars);
		cout << "end printing on line 630" << endl;
	}


	try
	{

		switch (params->substitution)
		{
			case BarvinokParameters::NoSubstitution:
			{
				string rat_filename = read_polyhedron_data.filename + ".rat";
				DelegatingSingleConeParameters *write_params =
						new DelegatingSingleConeParameters(*params);
				delete params;
				params = write_params;
				ConeConsumer *writing_consumer =
						new GeneratingFunctionWritingConeConsumer(rat_filename);
				if (Poly->projecting_up_transducer)
					writing_consumer = compose(Poly->projecting_up_transducer,
							writing_consumer);
				writing_consumer = compose(
						new PointsInParallelepipedComputingConeTransducer(
								write_params), writing_consumer);
				write_params->SetConsumer(writing_consumer);
				cerr << "Writing multivariate generating function to `"
						<< rat_filename << "'." << endl;
				listCone *cone;
				for (cone = Poly->cones; cone != NULL; cone = cone->rest)
					barvinokDecomposition_Single(cone, write_params);
				cerr << "Multivariate generating function written to `"
						<< rat_filename << "'." << endl;
				break;
			}
			case BarvinokParameters::PolynomialSubstitution:
				if (ehrhart_polynomial)
				{
					cerr
							<< "Computation of Ehrhart polynomials is only implemented "
							<< "for the exponential substitution (--exp)."
							<< endl;
					exit(1);
				}
				if (Poly->unbounded)
				{
					cerr << "The polyhedron is unbounded." << endl;
					exit(1);
				}
				if (read_polyhedron_data.assumeUnimodularCones[0] == 'n')
				{
					if (read_polyhedron_data.Memory_Save[0] == 'n')
					{
						listCone *decomposed_cones = decomposeCones(
								Poly->cones, not Poly->dualized, *params);
						freeListCone(Poly->cones);
						Poly->cones = decomposed_cones;
						// 	cerr << "Decomposed cones: " << endl;
						//	printListCone(Poly->cones, Poly->numOfVars);
						/* Compute points in parallelepipeds */
						computePointsInParallelepipeds(Poly->cones,
								Poly->numOfVars, params);
					}
					// Iterator through simplicial cones, DFS
					else
					{
						Standard_Single_Cone_Parameters *standard_params =
								new Standard_Single_Cone_Parameters(*params);
						delete params;
						params = standard_params;
						decomposeAndComputeResidue(Poly->cones,
								read_polyhedron_data.degree, false,
								*standard_params);
					}
				}
				break;
			case BarvinokParameters::ExponentialSubstitution:
				if (Poly->unbounded)
				{
					cerr << "The polyhedron is unbounded." << endl;
					exit(1);
				}
				if (read_polyhedron_data.dualApproach[0] == 'y')
				{
					cerr
							<< "Exponential substitution is not yet implemented for the homogenized version."
							<< endl;
					exit(1);
				} else
				{
					if (approx)
					{
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
						cerr << "Approximation code is not compiled in, sorry."
								<< endl;
						exit(1);
#endif
					} else if (ehrhart_polynomial)
					{
						Exponential_Ehrhart_Parameters *exp_param =
								new Exponential_Ehrhart_Parameters(*params);
						delete params;
						params = exp_param;
						mpq_vector ehrhart_coefficients =
								decomposeAndComputeEhrhartPolynomial(
										Poly->cones, *exp_param);

						//PolytopeValuation polytopeValuation(cones,
						//				PolytopeValuation::VertexRayCones,
						//				myParameters.Number_of_Variables, myParameters);
						//		ans1 = polytopeValuation.findVolume(PolytopeValuation::DeterminantVolume);

						for(listCone * t = Poly->cones; t; t = t->rest)
							printConeToFile(cout, Poly->cones, Poly->numOfVars);


						PolytopeValuation polytopeValuation(Poly->cones,
										PolytopeValuation::TriangulatedCones,
										Poly->numOfVars, *params);
						RationalNTL ans = polytopeValuation.findVolume(PolytopeValuation::DeterminantVolume);
						cout << "VOLUME BY DETERMINAT-VOLUME2 " << ans << endl;
						exit(1);
						cerr << endl << "Ehrhart polynomial: ";
						{
							unsigned int i;
							for (i = 0; i < ehrhart_coefficients.size(); i++)
							{
								if (ehrhart_coefficients[i] > 0)
									cout << " + " << ehrhart_coefficients[i]
											<< " * t^" << i;
								else if (ehrhart_coefficients[i] < 0)
									cout << " - " << abs(
											ehrhart_coefficients[i]) << " * t^"
											<< i;
							}
						}
						cout << endl << endl;
					} else
					{
						Exponential_Single_Cone_Parameters *exp_param =
								new Exponential_Single_Cone_Parameters(*params);
						delete params;
						params = exp_param;
						Integer number_of_lattice_points =
								decomposeAndComputeExponentialResidue(
										Poly->cones, *exp_param);
						params->deliver_number_of_lattice_points(
								number_of_lattice_points);
					}
				}
				break;
			default:
				cerr << "Unknown BarvinokParameters::substitution" << endl;
				abort();
		}

		if (read_polyhedron_data.grobner[0] == 'y')
		{

			Poly->cones = ProjectUp(Poly->cones,
					read_polyhedron_data.oldnumofvars, Poly->numOfVars,
					read_polyhedron_data.templistVec);
			Poly->numOfVars = read_polyhedron_data.oldnumofvars;

		}
		if (Print[0] == 'y')
		{
			cout << "printing on line 800" << endl;
			printListCone(Poly->cones, Poly->numOfVars);
			cout << "end printing on line 800" << endl;
		}

		//   printListVector(IntegralHull(Poly->cones,  inequalities, equations, Poly->numOfVars), Poly->numOfVars);
		if (read_polyhedron_data.Memory_Save[0] == 'n')
		{

			if (read_polyhedron_data.dualApproach[0] == 'n')
			{
				cerr << "Creating generating function.\n";
				//printListVector(templistVec, oldnumofvars); cerr << ProjU << endl;
				if (read_polyhedron_data.equationsPresent[0] == 'y')
				{
					Poly->cones = ProjectUp2(Poly->cones,
							read_polyhedron_data.oldnumofvars, Poly->numOfVars,
							read_polyhedron_data.AA, read_polyhedron_data.bb);
					Poly->numOfVars = read_polyhedron_data.oldnumofvars;
					cout << "747::after calling projectup2" << endl;
				}

				createGeneratingFunctionAsMapleInput(fileName, Poly->cones,
						Poly->numOfVars);
			}
			//printListCone(cones, Poly->numOfVars);

			cerr << "Printing decomposed cones to `decomposed_cones'." << endl;
			printListConeToFile("decomposed_cones", Poly->cones,
					Poly->numOfVars);

			if (read_polyhedron_data.dualApproach[0] == 'n')
			{
				cerr << "Starting final computation.\n";
				params->deliver_number_of_lattice_points(Residue(Poly->cones,
						Poly->numOfVars));
			}

			if (read_polyhedron_data.dualApproach[0] == 'y')
			{
				cerr << "Starting final computation.\n";
				//cerr << "output_cone: " << output_cone;
				switch (params->decomposition)
				{
					case BarvinokParameters::IrrationalPrimalDecomposition:
					case BarvinokParameters::IrrationalAllPrimalDecomposition:
					{
#ifdef HAVE_EXPERIMENTS
						ofstream out("func.rat");
						out << "HS := ";
						TrivialMonomialSubstitutionMapleOutput(out, Poly->cones, Poly->numOfVars);
						out << ";";
#else
						cerr
								<< "Trivial monomial substitution not compiled in, sorry."
								<< endl;
#endif
						break;
					}
					case BarvinokParameters::DualDecomposition:
						ResidueFunction(Poly->cones, Poly->numOfVars,
								print_flag, read_polyhedron_data.degree,
								output_cone, params);
						// ResidueFunction consumes cones.
						Poly->cones = NULL;
						break;
					default:
						assert(0);
				}
				//  Else we have already computed the residue.
			}
		}
	} catch (NotIrrationalException)
	{
		cerr << "Bug: Irrationalization failed" << endl;
		exit(1);
	};

	freeListVector(read_polyhedron_data.templistVec);
	freeListVector(read_polyhedron_data.matrix);
	delete Poly;

	if (read_polyhedron_data.rationalCone[0] == 'y')
	{
		cerr << endl << "Rational function written to " << argv[argc - 1]
				<< ".rat" << endl << endl;
		strcpy(command, "mv ");
		strcat(command, "simplify.sum ");
		strcat(command, argv[argc - 1]);
		strcat(command, ".rat");
		system_with_error_check(command);
	}

	if (printfile[0] == 'y')
	{
		cerr << endl << "Rational function written to " << argv[argc - 1]
				<< ".rat" << endl << endl;
		strcpy(command, "mv ");
		strcat(command, "func.rat ");
		strcat(command, argv[argc - 1]);
		strcat(command, ".rat");
		system_with_error_check(command);
	}
	if ((removeFiles[0] == 'y')
			&& (read_polyhedron_data.dualApproach[0] == 'n'))
	{

		strcpy(command, "rm -f ");
		strcat(command, fileName);
		strcat(command, ".ext");
		system_with_error_check(command);

		strcpy(command, "rm -f ");
		strcat(command, fileName);
		strcat(command, ".cdd");
		system_with_error_check(command);

		if (read_polyhedron_data.Memory_Save[0] == 'n')
		{
			strcpy(command, "rm -f ");
			strcat(command, fileName);
			strcat(command, ".maple");
			system_with_error_check(command);
		}

		strcpy(command, "rm -f ");
		strcat(command, fileName);
		strcat(command, ".ead");
		system_with_error_check(command);

	}

	//cerr << "Computation done. " << endl;

	params->total_time.stop();
	cerr << params->total_time;

	{
		// until we have a more sophisticated test script --mkoeppe
		ofstream totalTime("totalTime");
		totalTime << params->total_time.get_seconds() << " ("
				<< params->read_time.get_seconds() << "r" << ", "
				<< params->vertices_time.get_seconds() << "v" << ", "
				<< params->irrationalize_time.get_seconds() << "i" << ", "
				<< params->dualize_time.get_seconds() << "d" << ", "
				<< params->triangulate_time.get_seconds() << "t" << ", "
				<< params->decompose_time.get_seconds() << "b" << ")" << endl;
		ofstream stats("latte_stats");
		params->print_statistics(stats);
	}

	delete params;
	return volumeAnswers;
}//mainValuationDriver



VolumesContainer mainValuationDriverOLD(const char *argv[], int argc)
{

	set_program_name(argv[0]);

	BarvinokParameters params;
	VolumesContainer volumeAns;
#ifdef SUN
	struct tms tms_buf;
#endif
	float z;
	int i, numOfVars, numOfAllVars, degree = 1;
	unsigned int flags = 0, print_flag = 0, output_cone = 0;
	vec_ZZ dim, v, w;
	int oldnumofvars;
	vec_ZZ *generators = 0;
	char fileName[127], invocation[127], decompose[10], equationsPresent[10],
			assumeUnimodularCones[127], dualApproach[127], taylor[127],
			printfile[127], rationalCone[127], nonneg[127], Memory_Save[127],
			Save_Tri[127], Load_Tri[127], Print[127], inthull[127],
			cddstyle[127], grobner[127], removeFiles[127], command[127],
			maximum[127], Singlecone[127], LRS[127], Vrepresentation[127],
			dilation[127], minimize[127], binary[127], interior[127];
	char valuationType[127], valuationAlg[127], printLawrenceFunction[127];
	listVector *matrix = 0, *equations = 0, *inequalities = 0, *rays = 0, *endRays = 0, *tmpRays = 0,
			*matrixTmp = 0;
	vec_ZZ cost;
	listVector *templistVec = 0;
	listCone *cones = 0, *tmp = 0, *tmpcones = 0;

	latte_banner(cerr);

	z = 0;

	strcpy(invocation, "Invocation: ");
	strcat(invocation, argv[0]);
	strcat(invocation, " ");

	strcpy(Vrepresentation, "no");
	strcpy(interior, "no");
	strcpy(dilation, "no");
	strcpy(binary, "no");
	strcpy(Singlecone, "no");
	strcpy(removeFiles, "yes");
	strcpy(grobner, "no");
	strcpy(maximum, "no");
	strcpy(minimize, "no");
	strcpy(decompose, "yes");
	strcpy(dualApproach, "no");
	strcpy(equationsPresent, "no");
	strcpy(assumeUnimodularCones, "no");
	strcpy(printfile, "no");
	strcpy(taylor, "no");
	strcpy(rationalCone, "no");
	strcpy(nonneg, "no");
	strcpy(Memory_Save, "yes");
	strcpy(Save_Tri, "no");
	strcpy(Load_Tri, "no");
	strcpy(Print, "no");
	strcpy(inthull, "no");
	strcpy(cddstyle, "no");
	strcpy(LRS, "no");
	strcpy(valuationType, "volume");
	strcpy(valuationAlg, "all");
	strcpy(printLawrenceFunction, "no");

	flags |= DUAL_APPROACH;


	if ( argc == 1)
	{
		cout << "type --help to print a help menu." << endl;
		exit(1);
	}//if too few parameters.
	if ( argc == 2 && !strcmp(argv[1], "--help"))
	{
		cout << "usage: " << argv[0] << " [valuation  type] [valuation algorithm] <latte file>\n"
		     << "valuation types: volume\n"
		     << "  volume algorithm: [--lawrence  <--printLawrenceFunction> | --triangulate | --all]\n"
		     << "\n"
		     << "Example: " << argv[0] << " --lawrence --printLawrenceFunction file.latte\n"
			 << "         (will print the volume found by the Lawrence method along with the Lawrence rational function.)\n";

		exit(0);
	}//if need help.

	for (i = 1; i < argc - 1; i++) //i = 1...argc-1 because we do not want to process the file name. The file name could be "simp"
	{
		strcat(invocation, argv[i]);
		strcat(invocation, " ");

		if (strncmp(argv[i], "simp", 3) == 0)
			strcpy(rationalCone, "yes");

		if (strncmp(argv[i], "cdd", 3) == 0)
			strcpy(cddstyle, "yes");

		if (strncmp(argv[i], "lrs", 3) == 0)
			strcpy(LRS, "yes");

		if (strncmp(argv[i], "trisave", 7) == 0)
		{
			strcpy(Save_Tri, "yes");
			flags |= SAVE;
		}
		if (strncmp(argv[i], "triload", 7) == 0)
		{
			strcpy(Load_Tri, "yes");
			flags |= LOAD;
		}
		if ( strcmp(argv[i], "keepFiles") == 0)
			strcpy(removeFiles, "no");
		if (strncmp(argv[i], "--lawrence", 8) == 0)
			strcpy(valuationAlg, "--lawrence");
		if (strcmp(argv[i], "--triangulate") == 0)
			strcpy(valuationAlg, "triangulate");
		if (strcmp(argv[i], "--all") == 0)
			strcpy(valuationAlg, "all");
		if ( strcmp(argv[i], "--printLawrenceFunction") == 0)
			strcpy(printLawrenceFunction, "yes");

	}//for i.
#if 0
	/* if cdd and lrs are NOT called. */
	if ((cddstyle[0] == 'n') && (LRS[0] == 'n'))
	{//hit
		if (argc == 2)
		{
			strcpy(printfile, "yes");
			flags |= PRINT;//hit
		} else if (argc == 3)
		{
			if (strncmp(argv[1], "sim", 3) == 0)
				strcpy(rationalCone, "yes");
			else
				strcpy(taylor, "yes");
		}
		//
		//else
		//{
		//	cerr << "Too many arguments.  Check the manual for command line."
		//			<< endl;
		//	exit(1);
		//}
	}


	/* if cdd is called but and lrs is NOT called. */
	else if ((cddstyle[0] == 'y') && (LRS[0] == 'n'))
	{
		if (argc == 3)
		{
			strcpy(printfile, "yes");
			//flags |= PRINT; //hit
			flags &= ~PRINT; //don't print.
		} else if (argc == 4)
		{
			for (i = 1; i < 3; i++)
				if (strncmp(argv[i], "sim", 3) == 0)
					strcpy(rationalCone, "yes");
			if (rationalCone[0] == 'n')
				strcpy(taylor, "yes");
		} else
		{
			cerr << "Too many arguments.  Check the manual for command line."
					<< endl;
			exit(1);
		}
	}

	/* if cdd is NOT called and lrs is called. */
	else if ((cddstyle[0] == 'n') && (LRS[0] == 'y'))
	{
		if (argc == 3)
		{
			strcpy(printfile, "yes");
			flags |= PRINT;
		} else if (argc == 4)
		{
			for (i = 1; i < 3; i++)
				if (strncmp(argv[i], "sim", 3) == 0)
					strcpy(rationalCone, "yes");
			if (rationalCone[0] == 'n')
				strcpy(taylor, "yes");
		} else
		{
			cerr << "Too many arguments.  Check the manual for command line."
					<< endl;
			exit(1);
		}
	}

	/* if cdd and lrs are called. */
	else if ((cddstyle[0] == 'y') && (LRS[0] == 'y'))
	{
		if (argc == 4)
		{
			strcpy(printfile, "yes");
			flags |= PRINT;//hit
		} else if (argc == 5)
		{
			for (i = 1; i < 4; i++)
				if (strncmp(argv[i], "sim", 3) == 0)
					strcpy(rationalCone, "yes");
			// else if(strncmp(argv[3],"sim",3)==0) strcpy(rationalCone,"yes");
			if (rationalCone[0] == 'n')
				strcpy(taylor, "yes");
		} else
		{
			cerr << "Too many arguments.  Check the manual for command line."
					<< endl;
			exit(1);
		}
	}
#endif

	flags |= DUAL_APPROACH;//hit
	if (minimize[0] == 'y')
		strcpy(maximum, "yes");
	if (grobner[0] == 'y')
		strcpy(equationsPresent, "yes");
	if (binary[0] == 'y')
	{
		strcpy(maximum, "yes");
		strcpy(Memory_Save, "no");
	}
	if (maximum[0] == 'y')
		strcpy(Memory_Save, "no");
	if (printfile[0] == 'y')
		strcpy(Memory_Save, "no");
	if (rationalCone[0] == 'y')
		strcpy(Memory_Save, "no");
	if (printfile[0] == 'y')
		print_flag = 1;

	if (rationalCone[0] == 'y')
	{

		//HugInt digit(argv[1]);
		//conv(output_cone, digit.BigInt);
		// User can use only Mode one
		output_cone = 3;
	}
	int dilation_const = 1;

	if (dilation[0] == 'y')
		dilation_const = atoi(argv[argc - 2]);

	if (output_cone > 3)
		output_cone = 0;
	flags |= (output_cone << 1);

	if ((dualApproach[0] == 'y') && (nonneg[0] == 'y') && (equationsPresent[0]
			== 'n'))
	{
		cerr << "You cannot use + and dua at the same time." << endl;
		exit(2);
	}

	if ((cddstyle[0] == 'y') && (Vrepresentation[0] == 'y'))
	{
		cerr << "Use not cdd style and v-representation." << endl;
		exit(2);
	}

	if ((Memory_Save[0] == 'y') && (inthull[0] == 'y'))
	{
		cerr << "You cannot use int and memsave at the same time." << endl;
		exit(3);
	}

	strcat(invocation, argv[argc - 1]);
	strcat(invocation, "\n\n");
	cerr << invocation;
	char costFile[127];

	strcpy(fileName, argv[argc - 1]);

	/* Check input file. */
	if (Vrepresentation[0] == 'n')
	{//hit
		if ((cddstyle[0] == 'n') && (grobner[0] == 'n') && (maximum[0] == 'n')
				&& (minimize[0] == 'n'))
		{ //hit
			CheckInputFile(fileName);
			CheckLength(fileName, equationsPresent);
		}
		if (minimize[0] == 'y')
			strcpy(maximum, "yes");

		if ((cddstyle[0] == 'n') && (grobner[0] == 'n') && (maximum[0] == 'y'))
		{
			CheckInputFile(fileName);
			CheckLength(fileName, equationsPresent);
		}

		if (cddstyle[0] == 'y')
		{
			CheckInputFileCDDRep(argv[argc - 1]);
			CheckInputFileCDDRep1(argv[argc - 1]);
			CheckInputFileCDDRep3(argv[argc - 1]);
			CheckInputFileCDDRep4(argv[argc - 1]);
		}
	} else
		CheckInputFileVrep(fileName);
	CheckEmpty(fileName);//hit

	/* Read problem data. */
	if ((cddstyle[0] == 'n') && (Vrepresentation[0] == 'n'))//hit
		CheckRed(fileName, equationsPresent, maximum, nonneg, interior,
				dilation, dilation_const);
	//file, yes, no, no, no, no, 1

	dilation_const = 1;
	if ((cddstyle[0] == 'n') && (grobner[0] == 'n'))//hit
		readLatteProblem(fileName, &equations, &inequalities, equationsPresent,
				&numOfVars, nonneg, dualApproach, grobner, Vrepresentation);

	if (cddstyle[0] == 'y')
	{
		int tmpoutput;
		CDDstylereadLatteProblem(fileName, &equations, &inequalities,
				equationsPresent, &numOfVars, nonneg, dualApproach, taylor,
				degree, rationalCone, tmpoutput, Memory_Save,
				assumeUnimodularCones, inthull, grobner);
		output_cone = 3;
		if (dualApproach[0] == 'y')
		{
			flags |= DUAL_APPROACH;
		}
	}

	if ((dualApproach[0] == 'y') && (nonneg[0] == 'y') && (equationsPresent[0]
			== 'n'))
	{
		cerr << "You cannot use + and dua at the same time." << endl;
		exit(2);
	}

	if ((Memory_Save[0] == 'y') && (inthull[0] == 'y'))
	{
		cerr << "You cannot use int and memsave at the same time." << endl;
		exit(3);
	}

	if ((Vrepresentation[0] == 'y') && (equationsPresent[0] == 'y'))
	{
		cerr << "You cannot use vrep and equ at the same time." << endl;
		exit(4);
	}

	numOfVars--;
	numOfAllVars = numOfVars;
	mat_ZZ ProjU;
	ProjU.SetDims(numOfVars, numOfVars);
	oldnumofvars = numOfVars;
	generators = createArrayVector(numOfVars);
	if (equationsPresent[0] == 'y') //hit
	{
		matrixTmp = preprocessProblem(equations, inequalities, &generators,
				&numOfVars, cost, ProjU, interior, dilation_const);
		templistVec = transformArrayBigVectorToListVector(ProjU,
				ProjU.NumCols(), ProjU.NumRows());
	} else
	{
		dilateListVector(inequalities, numOfVars, dilation_const);
		matrixTmp = inequalities;
	}

	matrix = matrixTmp;

	/* Now matrix contains the new inequalities. */

	RR LP_OPT;
	cerr << "\nTime: " << GetTime() << " sec\n\n";
	vec_RR Rat_solution, tmp_den, tmp_num;
	mat_RR ProjU_RR;
	ProjU_RR.SetDims(ProjU.NumRows(), ProjU.NumCols());
	for (i = 0; i < ProjU.NumRows(); i++)
		for (int j = 0; j < ProjU.NumCols(); j++)
			conv(ProjU_RR[i][j], ProjU[i][j]);
	//cerr << ProjU << ProjU_RR << endl;
	Rat_solution.SetLength(numOfVars);
	tmp_den.SetLength(numOfVars);
	tmp_num.SetLength(numOfVars);

	/* Compute vertices and edges. */

	rationalVector* LP_vertex;
	if ((dualApproach[0] == 'n') && (Vrepresentation[0] == 'n'))
	{//hit
		if (LRS[0] == 'n')//hit
			tmpcones = computeVertexCones(fileName, matrix, numOfVars);
		else
			tmpcones = computeVertexConesViaLrs(fileName, matrix, numOfVars);

		cones = tmpcones;
		cerr << "\nThe polytope has " << lengthListCone(cones)
				<< " vertices.\n";
		//system_with_error_check("rm -f numOfLatticePoints");
		cerr << endl;
	}

	/* Reading from the vertex representation. */

	if (Vrepresentation[0] == 'y') //not hit
		cones = computeVertexConesFromVrep(fileName, numOfVars);

	/* Compute triangulation or decomposition of each vertex cone. */

	if (dualApproach[0] == 'n') // hit
	{
		if (assumeUnimodularCones[0] == 'n')
		{
			if (decompose[0] == 'y')
			{
				ofstream file;
				file.open("cones_count_oldMain.txt");
				for(listCone * t = cones; t; t = t->rest)
					printConeToFile(file, t, numOfVars);
				file.close();

				BarvinokParameters myParameters;
				myParameters.Flags = flags;
				myParameters.Number_of_Variables = numOfVars;
				myParameters.max_determinant = 1;
				myParameters.File_Name = fileName;
				myParameters.decomposition
						= BarvinokParameters::DualDecomposition;

				if (strcmp(valuationType, "volume") == 0)
				{
					volumeAns = computeVolume(cones, myParameters, valuationAlg, printLawrenceFunction);
				} else
				{
					cout << "Ops, " << valuationType << " is not supported" << endl;
					exit(1);
				}//if not finding volume.

			}//if decompose
		}//assumeUnimodularCones
	}//dualApproach


//	if (rationalCone[0] == 'y')
//	{
//		cerr << endl << "Rational function written to " << argv[argc - 1]
//				<< ".rat" << endl << endl;
//		strcpy(command, "mv ");
//		strcat(command, "simplify.sum ");
//		strcat(command, argv[argc - 1]);
//		strcat(command, ".rat");
//		system_with_error_check(command);
//	}

//	if (printfile[0] == 'y')
//	{
//		cerr << endl << "Rational function written to " << argv[argc - 1]
//				<< ".rat" << endl << endl;
//		strcpy(command, "mv ");
//		strcat(command, "func.rat ");
//		strcat(command, argv[argc - 1]);
//		strcat(command, ".rat");
//		system_with_error_check(command);
//	}

	if (removeFiles[0] == 'y')
	{

		strcpy(command, "rm -f ");
		strcat(command, fileName);
		strcat(command, ".ext");
		system_with_error_check(command);

		strcpy(command, "rm -f ");
		strcat(command, fileName);
		strcat(command, ".cdd");
		system_with_error_check(command);

		if (cddstyle[0] == 'n')
		{
			strcpy(command, "rm -f ");
			strcat(command, fileName);
			strcat(command, ".ead");
			system_with_error_check(command);

			strcpy(command, "rm -f ");
			strcat(command, fileName);
			system_with_error_check(command);
		}
	}// remove the latte_nonredundant_input files.


	//free up memory.
	//if ( matrix) 		freeListVector(matrix);
	if ( equations) 	freeListVector(equations);
	if (inequalities) 	freeListVector(inequalities);
	//if (rays) 			freeListVector(rays);
	//if (endRays) 		freeListVector(endRays);
	//if (tmpRays)		freeListVector(tmpRays);
	//if (matrixTmp)		freeListVector(matrixTmp);
	//if (templistVec)	freeListVector(templistVec);

	//if (generators)		delete generators;


	if (cones) 			freeListCone(cones);
	//if (tmp) 			freeListCone(tmp);
	//if (tmpcones) 		freeListCone(tmpcones);

	return (volumeAns);
}//mainValuationDrivver()


/**
 * Checks to see if the triangulation and lawrence volume equal the expected volume.
 */
void printVolumeTest(const RationalNTL &correctVolumeAnswer, VolumesContainer volumeAnswer, const string &file, const string &comments)
{
	if ( correctVolumeAnswer != volumeAnswer.lawrence
			|| correctVolumeAnswer != volumeAnswer.triangulate)
	{
		cout << "******* ERROR ******" << endl;
		cout << "correct answer: " << correctVolumeAnswer << endl;
		cout << "lawrence: " << volumeAnswer.lawrence << endl;
		cout << "triangulate: " << volumeAnswer.triangulate << endl;
		cout << "see file " << file.c_str() << endl;
		exit(1); //dont' delete the latte file.
	}//if error
	else
		cout << comments.c_str() << " CORRECT!" << endl;
}//printVolumeTest

/**
 * Calls polymake to make a random interger (or rational) vertex polytope, and then makes the latte file.
 * The latte file is then passed into mainValuationDriver() to find the volume
 *
 * We cannot check our volume with polymake for low-dimensional polytopes.
 */
void runOneTest(int ambientDim, int numPoints)
{
	const char * argv[] = {"runTests()", "--all", 0};
	stringstream comments;
	comments << "Making random integer polytope with " << numPoints << " points in R^" << ambientDim<< " for volume testing";

	BuildRandomPolytope buildPolytope(ambientDim);
	buildPolytope.setComments(comments.str().c_str());
	buildPolytope.setIntegerPoints(false); //make random rational points.
	buildPolytope.buildPolymakeFile(numPoints); //make the file
	buildPolytope.callPolymake(); //run polymake
	buildPolytope.findVolumeWithPolymake(); //run polymake for the volume
	buildPolytope.convertFacetEquations(); //fix facet equations
	buildPolytope.printFacetEquationsForLattE(); //make latte file.
	cout << comments.str().c_str();
	//buildPolytope.findEhrhardPolynomial();
	//buildPolytope.findVolumeWithPolymake();

	string file = buildPolytope.getLatteFile();


	char * sFile = new char[file.size() + 1];
	strcpy(sFile, file.c_str());
	argv[2] = sFile;
	mainValuationDriver(argv, 3);
	delete [] sFile;
}//RunOneTest

/**
 * Runs many random tests.
 */
void runTests()
{
	int startAmbientDim = 6, endAmbientDim = 50;
	int pointStepSize = 5;


	for(int ambientDim = startAmbientDim; ambientDim < endAmbientDim; ambientDim = ambientDim + 3)
	{
		for(int numberPoints = startAmbientDim / 2; numberPoints < startAmbientDim/2 + startAmbientDim + 10; numberPoints = numberPoints + pointStepSize)
			runOneTest(ambientDim, numberPoints);

	}//for ambientDim

}//runTests

/**
 * Finds the volume of hypersimplex polytopes and checks for correctness.
 */
void runHyperSimplexTests()
{
	const char * argv[] = {"runHyperSimplexTests()", "--all", 0};
						//   n  k  num/denom
	int  hyperSimplexData[][4] = { /*{4, 1, 1, 6},
							{4, 2, 2, 3},
							{5, 1, 1, 24},
							{5, 2, 11, 24},
							{6, 1, 1, 120},
							{6, 2, 13, 60},
							{6, 3, 11, 20},
							{7, 1, 1, 720},
							{7, 2, 19, 240},
							{7, 3, 151, 360},
							{8, 1, 1, 5040},
							{8, 2, 1, 42},
							{8, 3, 397, 1680},
							{8, 4, 151, 315},
							{9, 1, 1, 40320},
							{9, 2, 247, 40320},
							{9, 3, 477, 4480},
							{9, 4, 15619, 40320},
							{10, 1, 1, 362880},
							{10, 2, 251, 181440},
							{10, 3, 913, 22680},
							{10, 4, 44117, 181440},*/
							{10, 5, 15619, 36288}, //start here
							{11, 1, 1, 3628800},
							{11, 2, 1013, 3628800},
							{11, 3, 299, 22680},
							{11, 4, 56899, 45300},
							{11, 5, 655177, 1814400},
							{12, 1, 1, 39916800},
							{12, 2, 509, 9979200},
							{12, 3, 50879, 13305600},
							{12, 4, 1093, 19800},
							{12, 5, 1623019, 6652800},
							{12, 6, 655177, 1663200}
	};//hyperSimplexData

	int numberTestCases = 34;

	VolumesContainer volumeAnswer;
	for(int i = 0; i < numberTestCases; ++i)
	{
		stringstream comments;
		BuildHypersimplexEdgePolytope hyperSimplex(hyperSimplexData[i][0], hyperSimplexData[i][1]);

		comments << "finding volume of Hypersimplex(" << hyperSimplexData[i][0] << ", " << hyperSimplexData[i][1] << ")";
		hyperSimplex.buildPolymakeFile();

		hyperSimplex.setComments(comments.str().c_str());

		hyperSimplex.buildPolymakeFile(); //make the file
		hyperSimplex.callPolymake(); //run polymake
		//hyperSimplex.findVolumeWithPolymake(); //run polymake for the volume
		hyperSimplex.convertFacetEquations(); //fix facet equations
		hyperSimplex.printFacetEquationsForLattE(); //make latte file.

		string file = hyperSimplex.getLatteFile();

		cout << comments.str().c_str() << endl;

		char * sFile = new char[file.size() + 1];
		strcpy(sFile, file.c_str());
		argv[2] = sFile;
		volumeAnswer = mainValuationDriver(argv, 3);
		delete [] sFile;

		RationalNTL correctVolumeAnswer(hyperSimplexData[i][2], hyperSimplexData[i][3]);
		printVolumeTest(correctVolumeAnswer, volumeAnswer, file, comments.str());
	}//for i.
}//runHyperSimplexTests



/**
 * Finds the volume of Birkhoff polytopes and checks for correctness.
 */
void runBirkhoffTests()
{

	string birkhoff[] = { "../../EXAMPLES/birkhoff/birkhoff-5.latte",
						  "../../EXAMPLES/birkhoff/birkhoff-6.latte",
						  "../../EXAMPLES/birkhoff/birkhoff-7.latte"};
	string birkhoffVolume[][2] = { {"188723", "836911595520"}, //5
								   {"9700106723", "10258736801144832000000"}, //6
								   {"225762910421308831", "4709491654300668677115504230400000000"} //7
								 };
	int numberTestCases = 3;


    VolumesContainer volumeAnswer;
	const char * argv[] = {"runBirkhoffTests()", "--all", 0};


    for(int i = 0; i < numberTestCases; ++i)
    {
		char * sFile = new char[birkhoff[i].length() + 1];
		strcpy(sFile, birkhoff[i].c_str());
		argv[2] = sFile;
		volumeAnswer = mainValuationDriver(argv, 3);
		delete [] sFile;

		RationalNTL correctVolumeAnswer(birkhoffVolume[i][0], birkhoffVolume[i][1]);
		printVolumeTest(correctVolumeAnswer, volumeAnswer, string(birkhoff[i]), string("testing ") + string(birkhoff[i]));
    }//for ever file in the directory
}//runBirkhoffTests



int main(int argc, char *argv[])
{
	mainValuationDriver((const char **) argv, argc);
	//mainValuationDriverOLD((const char **) argv, argc); //remove the file printing later.
	//runHyperSimplexTests();
	//runBirkhoffTests();
	//runTests();
	//runOneTest(atoi(argv[1]), atoi(argv[2]));

	return 0;
}//main()
