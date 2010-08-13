/*
 * valuation.cpp
 *
 *  Created on: Jun 24, 2010
 *      Author: Brandon Dutra and Gregory Pinto
 */

#include "valuation.h"

/**
 * Computes the volume of the polytope.
 */
ValuationContainer Valuation::computeVolume(Polyhedron * poly,
		BarvinokParameters &myParameters, const char *valuationAlg,
		const char * print)
{
	ValuationContainer ans;
	RationalNTL ans1, ans2;
	Timer timer("");
	if (strcmp(valuationAlg, "triangulate") == 0 || strcmp(valuationAlg, "all")
			== 0)
	{
		PolytopeValuation polytopeValuation(poly, myParameters);
		timer.start();
		ans1 = polytopeValuation.findVolume(
				PolytopeValuation::DeterminantVolume);
		timer.stop();
		cout << "(using triangulation-determinant method) VOLUME " << ans1
				<< endl;
		cout << "Time for the triangulation-determinant method" << timer << endl;
	}//if triangulate

	if (strncmp(valuationAlg, "lawrence", 8) == 0
			|| strcmp(valuationAlg, "all") == 0)
	{

		PolytopeValuation polytopeValuation(poly, myParameters);
		timer.start();
		ans2 = polytopeValuation.findVolume(PolytopeValuation::LawrenceVolume);
		timer.stop();
		cout << "(using Lawrence method) VOLUME " << ans2 << endl;
		cout << "Time for the Lawrence method" << timer << endl;

		if (*print == 'y')
			polytopeValuation.printLawrenceVolumeFunction(); //print the lawrence rational function.

	}//if lawrence

	if (strcmp(valuationAlg, "all") == 0 && ans1 != ans2)
	{
		cerr << "Driver.cpp: the two methods are different." << ans1 << "!="
				<< ans2 << endl;
		exit(1);
	}//if error.
	else if (strncmp(valuationAlg, "lawrence", 8) == 0)
		cout << "Decimal: " << ans2.to_RR() << endl; //lawrence.
	else
		cout << "Decimal: " << ans1.to_RR() << endl; //triangulate or all.


	ans.lawrence = ans2;
	ans.triangulate = ans1;

	return ans;

}//computeVolume

/**
 * Computes the integral of a polynomial over the polytope.
 *
 * Currently returns a ValuationContainer object because we might
 * add Lawrence functionality in the future. --Brandon July 2010
 */
ValuationContainer Valuation::computeIntegral(Polyhedron *poly,
		BarvinokParameters &myParameters, const char *valuationAlg,
		const char *printLawrence, const char * polynomialString)
{
	ValuationContainer answer;
	RationalNTL ans1;
	monomialSum originalPolynomial;// polynomial without the updated coefficients.
	Timer timer("Integration time");
	PolytopeValuation polytopeValuation(poly, myParameters);

	loadMonomials(originalPolynomial, polynomialString); //get the polynomial from the string.
	timer.start();
	ans1 = polytopeValuation.integrate(originalPolynomial);
	timer.stop();
	cout << "integrate answer = " << ans1 << endl;
	cout << "Decimal answer = " << ans1.to_RR() << endl;
	cout << timer << endl;

	answer.triangulate = ans1;

	destroyMonomials(originalPolynomial);
	return answer;

}//computeIntegral


static void Valuation::usage(const char *progname)
{
	cerr << "usage: " << progname << " [OPTIONS...] " << "INPUTFILE" << endl;
	cerr << "Type `" << progname << " --help' "
			<< "for a list of options and input specifications." << endl;
}


/**
 * The main function and reads in the polytope and computes a valuation.
 */
ValuationContainer Valuation::mainValuationDriver(const char *argv[], int argc)
{
	ValuationContainer valuationAnswers;
	set_program_name(argv[0]);

	int i;
	unsigned int flags = 0, print_flag = 0, output_cone = 0;
	char printfile[127], Save_Tri[127], Load_Tri[127], Print[127],
			removeFiles[127], command[127];
	char valuationAlg[127], valuationType[127], printLawrence[127],
			polynomialFile[127];
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
	strcpy(Print, "no");
	strcpy(valuationAlg, "all");
	strcpy(printLawrence, "no");
	strcpy(polynomialFile, "");
	//strcpy(valuationType, "volume");
	strcpy(valuationType, "integrate");

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
					<< "  valuation types: --valuation=volume, or --valuation=integrate\n"
					<< "volume algorithms and options:\n"
					<< "  --lawrence  [--printLawrenceFunction]       Computes the volume using the Lawrence formula and\n"
					<< "                                              and prints the Lawrence rational function.\n"
					<< "  --triangulate                               Computes the volume using the triangulation method.\n"
					<< "  --all                                       Computes the volume using all the methods.\n"
					<< "\n" << "integration options\n"
					<< "  --monomials=<file>                           Looks at the first line of file for a polynomial\n"
					<< "                                              encoded in maple-syntax: [ [coef, [exponent vector]], ...]\n"
					<< "                                              If cannot open file, the line is read from std in.\n"
					<< "Example: " << argv[0]
					<< " --valuation=volume --lawrence --printLawrenceFunction file.latte\n"
					<< "         (will print the volume found by the Lawrence method along with the Lawrence rational function.)\n"
					<< "Example: " << argv[0]
					<< "--valuation=integrate --monomials=poly.txt file.latte\n"
					<< "         (will compute the integral of the polynomial in poly.txt over the polytope in file.latte.)\n"
					<< endl;
			exit(0);
		} else if (strcmp(argv[i], "--lawrence") == 0)
			strcpy(valuationAlg, "lawrence");
		else if (strcmp(argv[i], "--triangulate") == 0)
			strcpy(valuationAlg, "triangulate");
		else if (strcmp(argv[i], "--all") == 0)
			strcpy(valuationAlg, "all");
		else if (strcmp(argv[i], "--printLawrenceFunction") == 0)
			strcpy(printLawrence, "yes");
		else if (strcmp(argv[i], "--valuation=integrate") == 0)
			strcpy(valuationType, "integrate");
		else if (strcmp(argv[i], "--valuation=volume") == 0)
			strcpy(valuationType, "volume");
		else if (strncmp(argv[i], "--monomials=", 12) == 0)
		{
			if (strlen(argv[i]) > 127)
			{
				cerr << "polynomial file name is too long" << endl;
				exit(1);
			}
			strncpy(polynomialFile, argv[i] + 12, strlen(argv[i]) - 12 + 1);
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

	//ValuationContainer computeVolume(listCone * cones, BarvinokParameters &myParameters,
	//		const char *valuationAlg, const char * print)

	params->Flags = flags;
	params->Number_of_Variables = Poly->numOfVars;
	params->max_determinant = 1;
	params->File_Name = (char*) fileName;
	//params->File_Name = fileName;
	params->decomposition = BarvinokParameters::DualDecomposition;

	if (strcmp(valuationType, "volume") == 0)
		valuationAnswers = computeVolume(Poly, *params, valuationAlg,
				printLawrence);
	else if (strcmp(valuationType, "integrate") == 0) //add input of polynomial.
	{
		//read the polynomial from the file or from std in.
		ifstream inFile;
		istream inStream(cin.rdbuf());

		if (strcmp(polynomialFile, "") != 0)
		{
			inFile.open(polynomialFile);
			if (!inFile.is_open())
			{
				cerr << "Error: cannot open " << polynomialFile;
				exit(1);
			}
			inStream.rdbuf(inFile.rdbuf());
		}//set the inStream.


		string polynomialLine;
		polynomialLine = "";

		if (inFile.is_open())
		{
			cerr << "Reading polynomial from file " << polynomialFile << endl;
			getline(inStream, polynomialLine, '\n');
			inFile.close();
		} else
		{
			cerr << "Enter Polynomial >";
			getline(inStream, polynomialLine, '\n');
		}//user supplied polynomial in file.

		valuationAnswers = computeIntegral(Poly, *params,
				valuationAlg, printLawrence, polynomialLine.c_str());
	} else
	{
		cerr << "ops, valuation type is not known: " << valuationType << endl;
		exit(1);
	}//else error. This else block should not be reachable!

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

	params->total_time.stop();
	cerr << params->total_time;
	delete params;

	return valuationAnswers;
}//mainValuationDriver


/**
 * Checks to see if the triangulation and lawrence volume equal the expected volume.
 */
void VolumeTests::printVolumeTest(const RationalNTL &correctVolumeAnswer,
		ValuationContainer volumeAnswer, const string &file,
		const string &comments)
{
	if (correctVolumeAnswer != volumeAnswer.lawrence || correctVolumeAnswer
			!= volumeAnswer.triangulate)
	{
		cerr << "******* ERROR ******" << endl;
		cerr << "correct answer: " << correctVolumeAnswer << endl;
		cerr << "lawrence: " << volumeAnswer.lawrence << endl;
		cerr << "triangulate: " << volumeAnswer.triangulate << endl;
		cerr << "see file " << file.c_str() << endl;
		exit(1); //dont' delete the latte file.
	}//if error
	else
		cerr << comments.c_str() << " CORRECT!" << endl;
}//printVolumeTest

/**
 * Calls polymake to make a random interger (or rational) vertex polytope, and then makes the latte file.
 * The latte file is then passed into mainValuationDriver() to find the volume
 *
 * We cannot check our volume with polymake for low-dimensional polytopes.
 */
void VolumeTests::runOneTest(int ambientDim, int numPoints)
{
	const char * argv[] =
	{ "runTests()", "--valuation=volume", "--all", 0 };
	stringstream comments;
	comments << "Making random integer polytope with " << numPoints
			<< " points in R^" << ambientDim << " for volume testing";

	BuildRandomPolytope buildPolytope(ambientDim);
	buildPolytope.setComments(comments.str().c_str());
	buildPolytope.setIntegerPoints(false); //make random rational points.
	buildPolytope.buildPolymakeFile(numPoints); //make the file
	buildPolytope.callPolymake(); //run polymake
	buildPolytope.findVolumeWithPolymake(); //run polymake for the volume
	buildPolytope.convertFacetEquations(); //fix facet equations
	buildPolytope.printFacetEquationsForLattE(); //make latte file.
	cerr << comments.str().c_str();
	//buildPolytope.findEhrhardPolynomial();
	//buildPolytope.findVolumeWithPolymake();

	string file = buildPolytope.getLatteFile();

	char * sFile = new char[file.size() + 1];
	strcpy(sFile, file.c_str());
	argv[3] = sFile;
	Valuation::mainValuationDriver(argv, 4);
	delete[] sFile;
}//RunOneTest

/**
 * Runs many random tests by calling runOneTest
 */
void VolumeTests::runTests()
{
	int startAmbientDim = 6, endAmbientDim = 50;
	int pointStepSize = 5;

	for (int ambientDim = startAmbientDim; ambientDim < endAmbientDim; ambientDim
			= ambientDim + 3)
	{
		for (int numberPoints = startAmbientDim / 2; numberPoints
				< startAmbientDim / 4 + startAmbientDim; numberPoints
				= numberPoints + pointStepSize)
			runOneTest(ambientDim, numberPoints);

	}//for ambientDim

}//runTests

/**
 * Finds the volume of hypersimplex polytopes and checks for correctness.
 */
void VolumeTests::runHyperSimplexTests()
{
	const char * argv[] =
	{ "runHyperSimplexTests()", "--valuation=volume", "--all", 0 };
	//   n  k  num/denom
	int hyperSimplexData[][4] =
	{ {4, 1, 1, 6},
	 {4, 2, 2, 3},
	 /*{5, 1, 1, 24},
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
	{ 10, 5, 15619, 36288 }, //start here

			{ 11, 1, 1, 3628800 },
			{ 11, 2, 1013, 3628800 },
			{ 11, 3, 299, 22680 },
			{ 11, 4, 56899, 45300 },
			{ 11, 5, 655177, 1814400 },
			{ 12, 1, 1, 39916800 },
			{ 12, 2, 509, 9979200 },
			{ 12, 3, 50879, 13305600 },
			{ 12, 4, 1093, 19800 },
			{ 12, 5, 1623019, 6652800 },
			{ 12, 6, 655177, 1663200 } };//hyperSimplexData

	int numberTestCases = 34;

	ValuationContainer volumeAnswer;
	for (int i = 0; i < numberTestCases; ++i)
	{
		stringstream comments;
		BuildHypersimplexEdgePolytope hyperSimplex(hyperSimplexData[i][0],
				hyperSimplexData[i][1]);

		comments << "finding volume of Hypersimplex(" << hyperSimplexData[i][0]
				<< ", " << hyperSimplexData[i][1] << ")";
		hyperSimplex.buildPolymakeFile();

		hyperSimplex.setComments(comments.str().c_str());

		hyperSimplex.buildPolymakeFile(); //make the file
		hyperSimplex.callPolymake(); //run polymake
		//hyperSimplex.findVolumeWithPolymake(); //run polymake for the volume
		hyperSimplex.convertFacetEquations(); //fix facet equations
		hyperSimplex.printFacetEquationsForLattE(); //make latte file.

		string file = hyperSimplex.getLatteFile();

		cerr << comments.str().c_str() << endl;

		char * sFile = new char[file.size() + 1];
		strcpy(sFile, file.c_str());
		argv[3] = sFile;
		volumeAnswer = Valuation::mainValuationDriver(argv, 4);
		delete[] sFile;

		RationalNTL correctVolumeAnswer(hyperSimplexData[i][2],
				hyperSimplexData[i][3]);
		VolumeTests::printVolumeTest(correctVolumeAnswer, volumeAnswer, file, comments.str());
	}//for i.
}//runHyperSimplexTests


/**
 * Finds the volume of Birkhoff polytopes and checks for correctness.
 */
void VolumeTests::runBirkhoffTests()
{

	string birkhoff[] =
	{ "../../../../EXAMPLES/birkhoff/birkhoff-5.latte",
			"../../../../EXAMPLES/birkhoff/birkhoff-6.latte",
			"birkhoff7.latte.vrep" };
	string birkhoffVolume[][2] =
	{
	{ "188723", "836911595520" }, //5

			{ "9700106723", "10258736801144832000000" }, //6

			{ "225762910421308831", "4709491654300668677115504230400000000" } //7
	};
	int numberTestCases = 3;

	ValuationContainer volumeAnswer;
	const char * argv[] =
	{ "runBirkhoffTests()", "--valuation=volume", "--all", 0 };

	for (int i = 0; i < numberTestCases; ++i)
	{
		char * sFile = new char[birkhoff[i].length() + 1];
		strcpy(sFile, birkhoff[i].c_str());
		argv[3] = sFile;
		if ( i != 2)
			volumeAnswer = Valuation::mainValuationDriver(argv, 4);
		else
		{
			const char * argv2[5];
			argv2[0] = argv[0];
			argv2[1] = argv[1];
			argv2[2] = argv[2];
			argv2[3] = "--vrep";
			argv2[4] = sFile;
			volumeAnswer = Valuation::mainValuationDriver(argv2, 5);
		}//need to handel the v-rep file differently.
		delete[] sFile;

		RationalNTL correctVolumeAnswer(birkhoffVolume[i][0],
				birkhoffVolume[i][1]);
		VolumeTests::printVolumeTest(correctVolumeAnswer, volumeAnswer, string(birkhoff[i]),
				string("testing ") + string(birkhoff[i]));
	}//for ever file in the directory
}//runBirkhoffTests

