/*
 * valuation.cpp
 *
 *  Created on: Jun 24, 2010
 *      Author: Brandon Dutra and Gregory Pinto
 */

#include "valuation.h"
#include <iomanip>

//defines the input file/string.

/**
 * Computes the volume of the polytope.
 * The triangulation method will not change the org. polytope, but the
 * lawrence method will dilate it.
 */
Valuation::ValuationContainer Valuation::computeVolume(Polyhedron * poly,
		BarvinokParameters &myParameters, const char *valuationAlg,
		const char * print)
{
	ValuationContainer ans;
	RationalNTL ans1, ans2;

	if (strcmp(valuationAlg, "triangulate") == 0 || strcmp(valuationAlg, "all")
			== 0)
	{
		ValuationData timer_and_result;
		PolytopeValuation polytopeValuation(poly, myParameters);
		timer_and_result.timer.start();
		ans1 = polytopeValuation.findVolume(
				PolytopeValuation::DeterminantVolume);
		timer_and_result.timer.stop();

		timer_and_result.valuationType = ValuationData::volumeTriangulation;
		timer_and_result.answer = ans1;
		ans.add(timer_and_result);
	}//if triangulate. origional polytope should not have changed.

	if (strncmp(valuationAlg, "lawrence", 8) == 0
			|| strcmp(valuationAlg, "all") == 0)
	{
		ValuationData timer_and_result;
		PolytopeValuation polytopeValuation(poly, myParameters);
		timer_and_result.timer.start();
		ans2 = polytopeValuation.findVolume(PolytopeValuation::LawrenceVolume);
		timer_and_result.timer.stop();

		if (*print == 'y')
			polytopeValuation.printLawrenceVolumeFunction(); //print the lawrence rational function.

		timer_and_result.valuationType = ValuationData::volumeLawrence;
		timer_and_result.answer = ans2;
		ans.add(timer_and_result);
	}//if lawrence. Origional polytope is now dilated.

	if (strcmp(valuationAlg, "all") == 0 && ans1 != ans2)
	{
		cerr << "valuation.cpp: the two methods are different." << endl;
		cerr << "Lawrence:      " << ans2 << endl;
		cerr << "Triangulation: " << ans1 << endl;
		exit(1);
	}//if error.

	return ans;
}//computeVolume

/**
 * Computes the integral of a polynomial over the polytope.
 *
 * Both the triangulation and lawrence methods will dilate the org. polytope.
 */
Valuation::ValuationContainer Valuation::computeIntegral(Polyhedron *poly,
		BarvinokParameters &myParameters, const char *valuationAlg,
		const char * integrandString, const IntegrandType integrandType)
{
	ValuationContainer answer;
	ValuationData tiangulate_timer_and_result;
	ValuationData lawrence_timer_and_result;
	RationalNTL ans1, ans2;
	Polyhedron *poly2 = poly;//if doing both methods, make a deep copy of the origional polytopel.

	if (strcmp(valuationAlg, "all") == 0)
	{
		poly2 = new Polyhedron(*poly); //copy org. polytope, because it will be dilated.
	}

	if (strcmp(valuationAlg, "triangulate") == 0 || strcmp(valuationAlg, "all")
			== 0)
	{
		cout << "Going to run the triangulation integration method" << endl;
		PolytopeValuation polytopeValuation(poly, myParameters);

		if (integrandType == inputPolynomial)
		{
			monomialSum originalPolynomial;// polynomial without the updated coefficients.
			loadMonomials(originalPolynomial, integrandString); //get the polynomial from the string.

			tiangulate_timer_and_result.timer.start();
			ans1 = polytopeValuation.findIntegral(originalPolynomial,
					PolytopeValuation::TriangulationIntegration);
			tiangulate_timer_and_result.timer.stop();

			tiangulate_timer_and_result.valuationType
					= ValuationData::integrateTriangulation;
			tiangulate_timer_and_result.answer = ans1;
			answer.add(tiangulate_timer_and_result);

			destroyMonomials(originalPolynomial);
		}
		else
		{
			linFormSum originalLinearForm;
			loadLinForms(originalLinearForm, integrandString);

			tiangulate_timer_and_result.timer.start();
			ans1 = polytopeValuation.findIntegral(originalLinearForm,
					PolytopeValuation::TriangulationIntegration);
			tiangulate_timer_and_result.timer.stop();

			tiangulate_timer_and_result.valuationType
					= ValuationData::integrateTriangulation;
			tiangulate_timer_and_result.answer = ans1;
			answer.add(tiangulate_timer_and_result);

			destroyLinForms(originalLinearForm);

		}
	}//if doing triangulation method.


	if (strncmp(valuationAlg, "lawrence", 8) == 0
			|| strcmp(valuationAlg, "all") == 0)
	{
		cout << "Going to run the cone-decomposition integration method" << endl;
		if (integrandType == inputPolynomial)
		{
			monomialSum originalPolynomial;// polynomial without the updated coefficients.
			PolytopeValuation polytopeValuation(poly2, myParameters);

			loadMonomials(originalPolynomial, integrandString); //get the polynomial from the string.
			lawrence_timer_and_result.timer.start();
			ans2 = polytopeValuation.findIntegral(originalPolynomial,
					PolytopeValuation::LawrenceIntegration);
			lawrence_timer_and_result.timer.stop();

			lawrence_timer_and_result.valuationType
					= ValuationData::integrateLawrence;
			lawrence_timer_and_result.answer = ans2;
			answer.add(lawrence_timer_and_result);

			destroyMonomials(originalPolynomial);
		}
		else
		{
			linFormSum originalLinearForm;// polynomial without the updated coefficients.
			PolytopeValuation polytopeValuation(poly2, myParameters);

			loadLinForms(originalLinearForm, integrandString); //get the polynomial from the string.
			lawrence_timer_and_result.timer.start();
			ans2 = polytopeValuation.findIntegral(originalLinearForm,
					PolytopeValuation::LawrenceIntegration);
			lawrence_timer_and_result.timer.stop();

			lawrence_timer_and_result.valuationType
					= ValuationData::integrateLawrence;
			lawrence_timer_and_result.answer = ans2;
			answer.add(lawrence_timer_and_result);

			destroyLinForms(originalLinearForm);
		}
	}

	if (strcmp(valuationAlg, "all") == 0 && ans1 != ans2)
	{
		cerr << "Valuation.cpp: the two methods are different.\n"
				<< "triangulateion: " << ans1 << "\nlawrence       " << ans2
				<< endl;
		exit(1);
	}//if error.
	if (strcmp(valuationAlg, "all") == 0)
	{
		delete poly2;
	}//delete the copy we made.

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
Valuation::ValuationContainer Valuation::mainValuationDriver(
		const char *argv[], int argc)
{
	ValuationContainer valuationAnswers;
	set_program_name(argv[0]);

	int i;
	unsigned int flags = 0, print_flag = 0, output_cone = 0;
	char printfile[127], Save_Tri[127], Load_Tri[127], Print[127],
			removeFiles[127], command[127];
	char valuationAlg[127], valuationType[127], printLawrence[127],
			integrandFile[127];
	IntegrandType integrandType = nothing;
	bool approx;
	bool ehrhart_polynomial, ehrhart_series, ehrhart_taylor;
	bool triangulation_specified = false;
	bool useStokes = false;
	double sampling_factor = 1.0;
	long int num_samples = -1;
	ReadPolyhedronData read_polyhedron_data;
	//ReadPolyhedronDataRecursive read_polyhedron_data;

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
	strcpy(integrandFile, "");
	//strcpy(valuationType, "volume");
	strcpy(valuationType, "integrate");

	approx = false;
	ehrhart_polynomial = false;
	params->substitution = BarvinokParameters::PolynomialSubstitution;
	//params->decomposition = BarvinokParameters::DualDecomposition;
	params->decomposition
			= BarvinokParameters::IrrationalAllPrimalDecomposition;
	params->triangulation = BarvinokParameters::RegularTriangulationWithCdd;
	params->max_determinant = 1;

	for (i = 1; i < argc; i++)
	{
		//cout << "currently doing " << argv[i] << "." << endl;
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
					<< "  --cone-decompose                            Computes the volume using the Lawrence formula and\n"
					<< "                                              and prints the Lawrence rational function.\n"
					<< "  --triangulate                               Computes the volume using the triangulation method.\n"
					<< "  --all                                       Computes the volume using all the methods.\n"
					<< "\n" << "integration options\n"
					<< "  --monomials=<file>                           Looks at the first line of file for a polynomial\n"
					<< "                                              encoded in maple-syntax: [ [coef, [exponent vector]], ...]\n"
					<< "                                              If cannot open file, the line is read from std in.\n"
					<< "Example: " << argv[0]
					<< " --valuation=volume --cone-decompose --print-cone-decompose-function file.latte\n"
					<< "         (will print the volume found by the cone decomposition method along with the Lawrence rational function.)\n"
					<< "Example: " << argv[0]
					<< " --valuation=integrate --monomials=poly.txt file.latte\n"
					<< "         (will compute the integral of the polynomial in poly.txt over the polytope in file.latte.)\n"
					<< endl;
			exit(0);
		} else if (strcmp(argv[i], "--lawrence") == 0 || strcmp(argv[i], "--cone-decompose") == 0)
		{
			strcpy(valuationAlg, "lawrence");
			strcpy(read_polyhedron_data.dualApproach, "no");
		}
		else if (strcmp(argv[i], "--triangulate") == 0)
		{
			strcpy(valuationAlg, "triangulate");
			strcpy(read_polyhedron_data.dualApproach, "yes");
		}
		else if (strcmp(argv[i], "--all") == 0)
		{
			strcpy(valuationAlg, "all");
			strcpy(read_polyhedron_data.dualApproach, "no");
		}
		else if (strcmp(argv[i], "--print-lawrence-function") == 0)
			strcpy(printLawrence, "yes");
		else if (strcmp(argv[i], "--valuation=integrate") == 0)
			strcpy(valuationType, "integrate");
		else if (strcmp(argv[i], "--valuation=volume") == 0)
			strcpy(valuationType, "volume");
		else if ( strcmp(argv[i], "--stokes") == 0)
			useStokes = true;
		else if (strncmp(argv[i], "--monomials=", 12) == 0)
		{
			if (strlen(argv[i]) > 127)
			{
				cerr << "polynomial file name is too long" << endl;
				exit(1);
			}
			strncpy(integrandFile, argv[i] + 12, strlen(argv[i]) - 12 + 1);
			integrandType = inputPolynomial;
		} else if ( strncmp(argv[i], "--linear-forms=", 15) == 0 )
		{
			if ( strlen(argv[i]) > 127)
			{
				cerr << "linear form file name is too long " << endl;
				exit(1);
			}
			strncpy(integrandFile, argv[i] + 15, strlen(argv[i]) - 15 + 1);
			integrandType = inputLinearForm;
		}
		else if (read_polyhedron_data.parse_option(argv[i]))
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




	if ( useStokes == false)
	{
		Polyhedron *Poly = read_polyhedron_data.read_polyhedron(params);

		params->Flags = flags;
		params->Number_of_Variables = Poly->numOfVars;
		params->max_determinant = 1;
		params->File_Name = (char*) fileName;

		//poly is updated.
		polyhedronToCones(valuationAlg, Poly, params);

		//now the cones of poly are the tangent cones or the lifted cone of the vertices.
		if (strcmp(valuationType, "volume") == 0)
		{
			valuationAnswers = computeVolume(Poly, *params, valuationAlg,
					printLawrence);
		} else if (strcmp(valuationType, "integrate") == 0) //add input of polynomial.
		{
			//read the polynomial from the file or from std in.
			ifstream inFile;
			istream inStream(cin.rdbuf());

			if (integrandType != nothing)
			{
				inFile.open(integrandFile);
				if (!inFile.is_open())
				{
					cerr << "Error: cannot open " << integrandFile;
					exit(1);
				}
				inStream.rdbuf(inFile.rdbuf());
			}//set the inStream.


			string integrandLine;
			integrandLine = "";

			if (inFile.is_open())
			{
				cerr << "Reading " << (integrandType == inputPolynomial ? "polynomial" : "linear forms" ) << " from file " << integrandFile << endl;
				getline(inStream, integrandLine, '\n');
				inFile.close();
			} else
			{
				cerr << "Enter a 'p' or a 'l' if you want to integrate a"
					 << (Poly->homogenized ? Poly->numOfVars - 1 : Poly->numOfVars)
					 << " dimensional "
					 << "\n polynomial or a power of a linear form respectively"
					 << "\n followed by the integrand >";
				char pl = cin.get();
				if (pl == 'p')
					integrandType = inputPolynomial;
				else if (pl = 'l')
					integrandType = inputLinearForm;
				else
				{
					cerr << "The character " << pl << " is not a p or l" << endl;
					exit(1);
				}
				getline(inStream, integrandLine, '\n');
			}//user supplied polynomial in file.

			valuationAnswers = computeIntegral(Poly, *params, valuationAlg,
					integrandLine.c_str(), integrandType);
		} else
		{
			cerr << "ops, valuation type is not known: " << valuationType << endl;
			exit(1);
		}//else error. This else block should not be reachable!
		ValuationData totalValuationTimer;
		totalValuationTimer.valuationType = ValuationData::entireValuation;
		params->total_time.stop();
		totalValuationTimer.timer = params->total_time;
		valuationAnswers.add(totalValuationTimer);

		freeListVector(read_polyhedron_data.templistVec);
		freeListVector(read_polyhedron_data.matrix);
		delete Poly;
	}
	else
	{ //use strokes

		ReadPolyhedronDataRecursive rpdr(read_polyhedron_data);
		rpdr.readHrepMatrixFromFile(read_polyhedron_data.filename, params);

		RecursivePolytopeValuation rpv;
		rpv.setMaxRecursiveLevel(4);
		int d;
		cout << "min dim: ";
		cin >> d;
		if ( d >= 0)
		{
			rpv.setMinDimension(d);
			rpv.findVolume(rpdr, params);
		}
		else
		{
			d *= (-1);
			linFormSum linform;
			Polyhedron *Poly2 = rpdr.findTangentCones();
			vec_ZZ exp;
			exp.SetLength(Poly2->numOfVars);
			for(int i = 0; i < exp.length(); ++i)
				exp[i]=i+1;
			//exp[0]=0;
			//exp[1]=1;
			//exp[2]=2;
			//exp[3]=0;

			linform.termCount = 0;
			linform.varCount = Poly2->numOfVars;
			int powerFactorial = 1;
			for(int i = 1; i <= d; ++i)
				powerFactorial *= i;
			insertLinForm(RationalNTL(powerFactorial,1), d, exp, linform);

			rpdr.latticeInverse();

			PolytopeValuation pv(Poly2, *params);
			if ( rpdr.getFullDimensionCount() < Poly2->numOfVars)
			{
				cout << " getFull dim = " << rpdr.getFullDimensionCount() << endl;
				pv.setLatticeInverse(rpdr.getLatticeInverse(), rpdr.getLatticeInverseDilation());
				pv.setFullDimension(rpdr.getFullDimensionCount());
			}

			RationalNTL ans;
			//ans = pv.findIntegral(linform);
			ans = pv.findVolume(PolytopeValuation::LawrenceVolume);

			cout << "volume non-stokes, not full-dim" << ans << endl;
			destroyLinForms(linform);
		}
		cout << "valuation.cPP exit called" << endl;
		exit(1);

	}//else. use strokes.




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


	cout << params->triangulate_time << endl;
	cout << params->dualize_time << endl;
	valuationAnswers.printResults(cout);

	delete params;
	return valuationAnswers;
}//mainValuationDriver


void Valuation::polyhedronToCones(const char valuationAlg[], Polyhedron *Poly, BarvinokParameters * params)
{
	/*
		cout << "#####valuation.cpp, first printing of the cones" << endl;
		printListCone(Poly->cones, Poly->numOfVars);
		cout << "end valuation.cpp, first printing of the cones" << endl;
	*/

		assert(Poly->cones != NULL);
	//	cout << "poly.homogenized:" << Poly->homogenized << endl;
	//	cout << "poly.dualized:" << Poly->dualized << endl;
	//	cout << "read:dualApproach:" << read_polyhedron_data.dualApproach << endl;
	//	cout << "raed.input_dualized" << read_polyhedron_data.input_dualized << endl;

		if (strcmp(valuationAlg, "lawrence") == 0 || strcmp(valuationAlg, "all") == 0 )
		{
			assert(Poly->homogenized == false);
			if (Poly->dualized)
			{
				cerr << "(First dualizing back... ";
				cerr.flush();
				dualizeCones(Poly->cones, Poly->numOfVars, params);
				cerr << "done.) ";
				cerr.flush();
			}
			if (Poly->cones->rays == NULL)
			{
				// dualize twice tocompute the rays.
				cerr << "(First computing their rays... ";
				cerr.flush();
				dualizeCones(Poly->cones, Poly->numOfVars, params);
				dualizeCones(Poly->cones, Poly->numOfVars, params); // just swaps
				cerr << "done!) ";
				cerr.flush();
			}

		}//find vertex-rays
		else if ( strcmp(valuationAlg, "triangulate") == 0)
		{
			assert(Poly->homogenized == true);
			if (Poly->dualized)
			{
				cerr << "(First dualizing back... ";
				cerr.flush();
				dualizeCones(Poly->cones, Poly->numOfVars, params);
				cerr << "done!) ";
				cerr.flush();
			}
		}//only need vertices
		else
		{
			cerr << "Valuation.cpp:: unknown valuation type: " << valuationAlg << '.' << endl;
			exit(1);
		}//else error.

		/*
			cout << "#####valuation.cpp, 2nd printing of the cones" << endl;
			printListCone(Poly->cones, Poly->numOfVars);
			cout << "end valuation.cpp, 2nd printing of the cones" << endl;
		*/



	/*  //OLD WAY OF SETTING UP THE POLYTOPE...NEW WAY IS ABOVE.
		if ( read_polyhedron_data.Vrepresentation[0] != 'y')
		{
		switch (params->decomposition)
		{
			case BarvinokParameters::DualDecomposition:
			case BarvinokParameters::IrrationalPrimalDecomposition:
				if (not Poly->dualized)
				{
					if (read_polyhedron_data.Vrepresentation[0] != 'y')
					{
						// Compute all inequalities tight at the respective vertex.
						// Then dualizeCones just needs to swap rays and facets.
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
						// Only facets computed, for instance by using the 4ti2
						// method of computing vertex cones.  So dualize twice to
						// compute the rays.
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
				//irrationalizeCones(Poly->cones, Poly->numOfVars);
				//params->irrationalize_time.stop();
				//cerr << params->irrationalize_time;
				break;
			default:
				cerr << "Unknown BarvinokParameters::decomposition" << endl;
				abort();
		}//switch
		}//if lawrence (or all) and h-rep.
	*/

		/*
		cout << "#####valuation.cpp, cones after" << endl;
		printListCone(Poly->cones, Poly->numOfVars);
		cout << "end valuation.cpp, cones after" << endl;
		exit(1);
		//*/
}//polyhedronToCones

Valuation::ValuationData::ValuationData() :
	timer(string(""), false)
{
}

void Valuation::ValuationContainer::add(const ValuationData & d)
{
	answers.push_back(d);
}

void Valuation::ValuationContainer::printResults(ostream & out) const
{
	out << "\n";
	for (int i = 0; i < answers.size(); ++i)
	{
		if (answers[i].valuationType == ValuationData::volumeLawrence)
			out << "Volume (using the Lawrence method)" << endl;
		else if (answers[i].valuationType == ValuationData::volumeTriangulation)
			out << "Volume (using the triangulation-determinant method)"
					<< endl;
		else if (answers[i].valuationType
				== ValuationData::integrateTriangulation)
			out << "Integration (using the triangulation method)" << endl;
		else if (answers[i].valuationType == ValuationData::integrateLawrence)
			out << "Integration (using the Lawrence method)" << endl;
		else if (answers[i].valuationType == ValuationData::entireValuation)
		{
			out
					<< "Computational time (algorithms + processing + program control)"
					<< endl;
			out << "     " << answers[i].timer;
			continue;
		}

		RR decimalAns;
		decimalAns = answers[i].answer.to_RR();
		decimalAns.SetOutputPrecision(32);

		out << "     Answer: " << answers[i].answer << endl;
		out << "     Decimal: " << decimalAns << endl;
		out << "     Time" << answers[i].timer;
	}//for each valuation entry.

}//printResults


