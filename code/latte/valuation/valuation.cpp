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
		BarvinokParameters &myParameters, const IntegrationInput &intInput,
		const char * print)
{
	ValuationContainer ans;
	RationalNTL ans1, ans2;
	Polyhedron * poly2;

	if ( intInput.all)
		poly2 = new Polyhedron(*poly);
	else
		poly2 = poly;

	if (intInput.volumeTriangulation || intInput.all)
	{
		ValuationData timer_and_result;
		PolytopeValuation polytopeValuation(poly, myParameters);
		timer_and_result.timer.start();
		ans1 = polytopeValuation.findVolume(
				PolytopeValuation::volumeTriangulation);
		timer_and_result.timer.stop();

		timer_and_result.valuationType = PolytopeValuation::volumeTriangulation;
		timer_and_result.answer = ans1;
		ans.add(timer_and_result);
	}//if triangulate. origional polytope should not have changed.

	if (intInput.volumeCone || intInput.all)
	{
		ValuationData timer_and_result;
		PolytopeValuation polytopeValuation(poly2, myParameters);
		timer_and_result.timer.start();
		ans2 = polytopeValuation.findVolume(PolytopeValuation::volumeCone);
		timer_and_result.timer.stop();

		if (*print == 'y')
			polytopeValuation.printLawrenceVolumeFunction(); //print the lawrence rational function.

		timer_and_result.valuationType = PolytopeValuation::volumeCone;
		timer_and_result.answer = ans2;
		ans.add(timer_and_result);
	}//if lawrence. Origional polytope is now dilated.

	if (intInput.all && ans1 != ans2)
	{
		cerr << "valuation.cpp: the two volume methods are different." << endl;
		cerr << "Lawrence:      " << ans2 << endl;
		cerr << "Triangulation: " << ans1 << endl;
		exit(1);
	}//if error.
	else if (intInput.all)
		delete poly2;

	return ans;
}//computeVolume

/**
 * Computes the integral over the polytope.
 *
 * Both the triangulation and lawrence methods will dilate/change the org. polytope.
 */
Valuation::ValuationContainer Valuation::computeIntegral(Polyhedron *poly,
		BarvinokParameters &myParameters, const IntegrationInput & intInput)
{
	if ( intInput.integrandType == IntegrationInput::inputPolynomial)
	{
		return computeIntegralPolynomial(poly, myParameters, intInput);
	} else if ( intInput.integrandType == IntegrationInput::inputLinearForm)
	{
		return computeIntegralLinearForm(poly, myParameters, intInput);
	} else if ( intInput.integrandType == IntegrationInput::inputProductLinearForm)
	{
		return computeIntegralProductLinearForm(poly, myParameters, intInput);
	}
	else
	{
		THROW_LATTE_MSG(LattException::bug_Unknown, "integrand type not supported.");
	}
}
Valuation::ValuationContainer Valuation::computeIntegralPolynomial(Polyhedron *poly,
		BarvinokParameters &myParameters, const IntegrationInput & intInput)
{
	ValuationContainer answer;
	ValuationData tiangulate_timer_and_result;
	ValuationData lawrence_timer_and_result;
	ValuationData plf_time_and_result;
	RationalNTL ans1, ans2, ans3, ans4;
	Polyhedron *polyCopy;//if doing more than 1 method, make a deep copy of the origional polytopel.

	assert(intInput.integrandType == IntegrationInput::inputPolynomial);

	if (intInput.integratePolynomialAsLinearFormTriangulation || intInput.all )
	{
		if(intInput.all)
			polyCopy = new Polyhedron(*poly);
		else
			polyCopy = poly;
		cerr << "Going to run the triangulation integration method" << endl;
		PolytopeValuation polytopeValuation(polyCopy, myParameters);

		monomialSum originalPolynomial;// polynomial without the updated coefficients.
		loadMonomials(originalPolynomial, intInput.integrand); //get the polynomial from the string.

		tiangulate_timer_and_result.timer.start();
		ans1 = polytopeValuation.findIntegral(originalPolynomial,
				PolytopeValuation::integratePolynomialAsLinearFormTriangulation);
		tiangulate_timer_and_result.timer.stop();

		tiangulate_timer_and_result.valuationType
					= PolytopeValuation::integratePolynomialAsLinearFormTriangulation;
		tiangulate_timer_and_result.answer = ans1;
		answer.add(tiangulate_timer_and_result);

		destroyMonomials(originalPolynomial);

		if (intInput.all)
			delete polyCopy;
	}//if doing triangulation method.


	if (intInput.integratePolynomialAsLinearFormCone || intInput.all)
	{
		cerr << "Going to run the cone-decomposition integration method" << endl;

		if(intInput.all)
			polyCopy = new Polyhedron(*poly);
		else
			polyCopy = poly;

		monomialSum originalPolynomial;// polynomial without the updated coefficients.
		PolytopeValuation polytopeValuation(polyCopy, myParameters);

		loadMonomials(originalPolynomial, intInput.integrand); //get the polynomial from the string.
		lawrence_timer_and_result.timer.start();
		ans2 = polytopeValuation.findIntegral(originalPolynomial,
				 PolytopeValuation::integratePolynomialAsLinearFormCone);
		lawrence_timer_and_result.timer.stop();

		lawrence_timer_and_result.valuationType
					= PolytopeValuation::integratePolynomialAsLinearFormCone;
		lawrence_timer_and_result.answer = ans2;
		answer.add(lawrence_timer_and_result);

		destroyMonomials(originalPolynomial);

		if(intInput.all)
			delete polyCopy;
	}

	if ( intInput.integratePolynomialAsPLFTriangulation || intInput.all)
	{
		cerr << "Going to run the polynomial to PLF method" << endl;

		if(intInput.all)
			polyCopy = new Polyhedron(*poly);
		else
			polyCopy = poly;

		monomialSum originalPolynomial;// polynomial without the updated coefficients.
		PolytopeValuation polytopeValuation(polyCopy, myParameters);

		loadMonomials(originalPolynomial, intInput.integrand); //get the polynomial from the string.
		plf_time_and_result.timer.start();
		ans3 = polytopeValuation.findIntegral(originalPolynomial,
				 PolytopeValuation::integratePolynomialAsPLFTriangulation);
		plf_time_and_result.timer.stop();

		plf_time_and_result.valuationType
					= PolytopeValuation::integratePolynomialAsPLFTriangulation;
		plf_time_and_result.answer = ans3;
		answer.add(plf_time_and_result);

		destroyMonomials(originalPolynomial);

		if(intInput.all)
			delete polyCopy;
	}//polynomial to PLF.


	if (intInput.all && (ans1 != ans2 || ans1 != ans3 ) )
	{
		cerr << "Valuation.cpp: the methods are different.\n"
				<< "triangulateion    : " << ans1 << "\n"
				<< "cone-decomposition: " << ans2 << "\n"
				<< "plf               : " << ans3 << "\n"
				<< endl;
		THROW_LATTE(LattException::bug_Unknown);
	}//if error.


	return answer;
}//computeIntegral



Valuation::ValuationContainer Valuation::computeIntegralLinearForm(Polyhedron *poly,
		BarvinokParameters &myParameters, const IntegrationInput & intInput)
{
	ValuationContainer answer;
	ValuationData tiangulate_timer_and_result;
	ValuationData lawrence_timer_and_result;
	ValuationData product_time_and_result;
	RationalNTL ans1, ans2;
	Polyhedron *poly2 = poly;//if doing both methods, make a deep copy of the origional polytopel.

	assert(intInput.integrandType == IntegrationInput::inputLinearForm);

	if (intInput.all)
	{
		poly2 = new Polyhedron(*poly); //copy org. polytope, because it will be dilated.
	}

	if (intInput.integrateLinearFormTriangulation  || intInput.all )
	{
		cerr << "Going to run the triangulation integration method on linear forms" << endl;
		PolytopeValuation polytopeValuation(poly, myParameters);


		linFormSum originalLinearForm;
		loadLinForms(originalLinearForm, intInput.integrand);

		tiangulate_timer_and_result.timer.start();
		ans1 = polytopeValuation.findIntegral(originalLinearForm,
					PolytopeValuation::integrateLinearFormTriangulation);
		tiangulate_timer_and_result.timer.stop();

		tiangulate_timer_and_result.valuationType
					= PolytopeValuation::integrateLinearFormTriangulation;
		tiangulate_timer_and_result.answer = ans1;
		answer.add(tiangulate_timer_and_result);

		destroyLinForms(originalLinearForm);
	}//if doing triangulation method.


	if (intInput.integrateLinearFormCone  || intInput.all)
	{
		cerr << "Going to run the cone-decomposition integration method on linear forms" << endl;

		linFormSum originalLinearForm;// polynomial without the updated coefficients.
		PolytopeValuation polytopeValuation(poly2, myParameters);

		loadLinForms(originalLinearForm, intInput.integrand); //get the polynomial from the string.
		lawrence_timer_and_result.timer.start();
		ans2 = polytopeValuation.findIntegral(originalLinearForm,
					PolytopeValuation::integrateLinearFormCone);
		lawrence_timer_and_result.timer.stop();

		lawrence_timer_and_result.valuationType
					= PolytopeValuation::integrateLinearFormCone;
		lawrence_timer_and_result.answer = ans2;
		answer.add(lawrence_timer_and_result);

		destroyLinForms(originalLinearForm);
	}

	if (intInput.all && ans1 != ans2)
	{
		cerr << "computeIntegralLinearForm(): the two methods are different.\n"
				<< "triangulateion: " << ans1 << "\nlawrence       " << ans2
				<< endl;
		THROW_LATTE(LattException::bug_Unknown);
	}//if error.
	if (intInput.all)
	{
		delete poly2;
	}//delete the copy we made.

	return answer;
}//computeIntegralLinearForm

Valuation::ValuationContainer Valuation::computeIntegralProductLinearForm(Polyhedron *poly,
		BarvinokParameters &myParameters, const IntegrationInput & intInput)
{
	ValuationContainer answer;
	ValuationData product_time_and_result;
	RationalNTL ans1;


	assert(intInput.integrandType == IntegrationInput::inputProductLinearForm);

	cerr << "Going to run the product of linear forms method" << endl;
	PolytopeValuation polytopeValuation(poly, myParameters);

	linFormProductSum originalProducts;
	loadLinFormProducts(originalProducts, intInput.integrand);


	//cout << "integrand string:" << intInput.integrand.c_str() << endl;
	//cout << "integrand load  :" << printLinFormProducts(originalProducts).c_str();
	//exit(1);




	product_time_and_result.timer.start();
	ans1 = polytopeValuation.findIntegral(originalProducts, PolytopeValuation::integrateProductLinearFormsTriangulation);
	product_time_and_result.timer.stop();

	product_time_and_result.valuationType = PolytopeValuation::integrateProductLinearFormsTriangulation;
	product_time_and_result.answer = ans1;
	answer.add(product_time_and_result);

	destroyLinFormProducts(originalProducts);
	return answer;
}//computeIntegralProductLinearForm

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
	IntegrationInput integrationInput;
	set_program_name(argv[0]);

	int i;
	unsigned int flags = 0, print_flag = 0, output_cone = 0;
	char printfile[127], Save_Tri[127], Load_Tri[127], Print[127],
			removeFiles[127], command[127];
	char printLawrence[127];
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
	strcpy(printLawrence, "no");
	integrationInput.integrandType = IntegrationInput::nothing;
	integrationInput.all = true;
	integrationInput.numEhrhartCoefficients = -1; // incremental mode
	integrationInput.realDilations = false; // don't output
						// formulas valid for
						// real dilations also

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
			     		<< "  --valuation=volume                          Computes the volume\n"
			     		<< "  --valuation=integrate                       Computes an integral\n"
	       			     	<< "  --valuation=top-ehrhart                     Computes the top weighted Ehrhart coefficients\n"
					<< "volume algorithms and options:\n"
					<< "  --cone-decompose                            Computes the volume using the Lawrence formula and\n"
					<< "                                              and prints the Lawrence rational function.\n"
					<< "  --triangulate                               Computes the volume using the triangulation method.\n"
					<< "  --all                                       Computes the volume using all the methods.\n"
					<< "\n" << "integration and weighted summation options:\n"
					<< "  --monomials=<file>                          Looks at the first line of file for a polynomial\n"
					<< "                                              encoded in maple-syntax: [ [coef, [exponent vector]], ...]\n"
			  /// FIXME: Add further options!
					<< "                                              If cannot open file, the line is read from std in.\n"
					<< "top Ehrhart options:\n"
					<< "  --num-coefficients=K                        Number of highest Ehrhart coefficients to compute\n"
					<< "                                              (default: compute all, incrementally)\n"
					<< "  --real-dilations                            Output formulas valid for real dilations also\n"
					<< "Example: " << argv[0]
					<< " --valuation=volume --cone-decompose --print-cone-decompose-function file.latte\n"
					<< "         (will print the volume found by the cone decomposition method along with the Lawrence rational function.)\n"
					<< "Example: " << argv[0]
					<< " --valuation=integrate --monomials=poly.txt file.latte\n"
					<< "         (will compute the integral of the polynomial in poly.txt over the polytope in file.latte.)\n"
					<< endl;
			exit(0);
		} else if ( strncmp(argv[i], "--valuation-alg=", 15) == 0)
		{
			integrationInput.all = false;
			if ( strcmp(argv[i], "--valuation-alg=volume-cone") == 0)
			{
				integrationInput.volumeCone = true;
				integrationInput.integrandType = IntegrationInput::inputVolume;
				strcpy(read_polyhedron_data.dualApproach, "no");
			}
			else if ( strcmp(argv[i], "--valuation-alg=volume-triangulation") == 0)
			{
				integrationInput.volumeTriangulation = true;
				integrationInput.integrandType = IntegrationInput::inputVolume;
				strcpy(read_polyhedron_data.dualApproach, "yes");
			}
			else if ( strcmp(argv[i], "--valuation-alg=poly-lf-triangulation") == 0)
			{
				integrationInput.integratePolynomialAsLinearFormTriangulation = true;
				strcpy(read_polyhedron_data.dualApproach, "yes");
			}
			else if ( strcmp(argv[i], "--valuation-alg=poly-lf-cone") == 0)
			{
				integrationInput.integratePolynomialAsLinearFormCone = true;
				strcpy(read_polyhedron_data.dualApproach, "no");
			}
			else if ( strcmp(argv[i], "--valuation-alg=poly-plf-triangulation") == 0)
			{
				integrationInput.integratePolynomialAsPLFTriangulation = true;
				strcpy(read_polyhedron_data.dualApproach, "yes");
			}
			else if ( strcmp(argv[i], "--valuation-alg=lf-triangulation") == 0)
			{
				integrationInput.integrateLinearFormTriangulation = true;
				strcpy(read_polyhedron_data.dualApproach, "yes");
			}
			else if ( strcmp(argv[i], "--valuation-alg=lf-cone") == 0)
			{
				integrationInput.integrateLinearFormCone = true;
				strcpy(read_polyhedron_data.dualApproach, "no");
			}
			else if ( strcmp(argv[i], "--valuation-alg=plf-triangulation") == 0)
			{
				integrationInput.integrateProductLinearFormsTriangulation = true;
				strcpy(read_polyhedron_data.dualApproach, "yes");
			}
			else if ( strcmp(argv[i], "--valuation-alg=volume") == 0)
			{
				integrationInput.all = true;
				integrationInput.integrandType = IntegrationInput::inputVolume;
				strcpy(read_polyhedron_data.dualApproach, "no");
			}
			else if ( strcmp(argv[i], "--valuation-alg=integrate") == 0)
			{
				integrationInput.all = true;
				strcpy(read_polyhedron_data.dualApproach, "no");
			}
		}else if (strcmp(argv[i], "--valuation=integrate") == 0)
		{
			cout << "TODO: add the old style integration optinos back." << endl;
			exit(1);
		} else if (strcmp(argv[i], "--valuation=volume") == 0)
		{
			cout << "TODO: add the old style volume optinos back." << endl;
			exit(1);
		}
		else if (strcmp(argv[i], "--valuation=top-ehrhart") == 0) {
		  strcpy(read_polyhedron_data.dualApproach, "yes");  // Want
								     // vertices
								     // of polytope
		  integrationInput.all = false;
		  integrationInput.topEhrhart = true;
		}
		else if (strncmp(argv[i], "--num-coefficients=", 19) == 0) {
		  integrationInput.numEhrhartCoefficients = atoi(argv[i] + 19);
		}
		else if (strncmp(argv[i], "--real-dilations", 6) == 0) {
		  integrationInput.realDilations = true;
		}
		/* else if (strcmp(argv[i], "--lawrence") == 0 || strcmp(argv[i], "--cone-decompose") == 0)
		{
			strcpy(valuationAlg, "lawrence");
			strcpy(read_polyhedron_data.dualApproach, "no");
		}
		else if (strcmp(argv[i], "--triangulate") == 0)
		{
			strcpy(valuationAlg, "triangulate");
			strcpy(read_polyhedron_data.dualApproach, "yes");
		}
		else if ( strcmp(argv[i], "--products-powers") == 0)
		{
			strcpy(valuationAlg, "products-powers");
			strcpy(read_polyhedron_data.dualApproach, "yes");
		}
		else if (strcmp(argv[i], "--all") == 0)
		{
			strcpy(valuationAlg, "all");
			strcpy(read_polyhedron_data.dualApproach, "no");
		}
		else if (strcmp(argv[i], "--print-lawrence-function") == 0 || strcmp(argv[i], "--print-cone-decompose-function") == 0)
			strcpy(printLawrence, "yes");
		else if ( strcmp(argv[i], "--stokes") == 0)
			useStokes = true;
		*/else if (strncmp(argv[i], "--monomials=", 12) == 0)
		{
			integrationInput.integrandType = IntegrationInput::inputPolynomial;
			integrationInput.fileName = (argv[i] + 12);
		} else if ( strncmp(argv[i], "--linear-forms=", 15) == 0 )
		{
			integrationInput.integrandType = IntegrationInput::inputLinearForm;
			integrationInput.fileName = (argv[i] + 15);
		} else if ( strncmp(argv[i], "--product-linear-forms=", 23) == 0)
		{
			integrationInput.integrandType = IntegrationInput::inputProductLinearForm;
			integrationInput.fileName = (argv[i] + 23);
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




	if ( useStokes == false)
	{
		Polyhedron *Poly = read_polyhedron_data.read_polyhedron(params);

		params->Flags = flags;
		params->Number_of_Variables = Poly->numOfVars;
		params->max_determinant = 1;
		params->File_Name = (char*) fileName;

		//poly is updated.
		polyhedronToCones(integrationInput, Poly, params);

		if (integrationInput.topEhrhart) {
		  valuationAnswers = computeTopEhrhart(Poly, *params, integrationInput);
		}
		//now the cones of poly are the tangent cones or the lifted cone of the vertices.
		else if (integrationInput.integrandType == IntegrationInput::inputVolume)
		{
			valuationAnswers = computeVolume(Poly, *params, integrationInput,
					printLawrence);
		} else //integration
		{
			//read the integrand from the file or from std in.
			ifstream inFile;
			istream inStream(cin.rdbuf());

			if (integrationInput.integrandType != IntegrationInput::nothing)
			{
				inFile.open(integrationInput.fileName.c_str());
				if (!inFile.is_open())
				{
					cerr << "Error: cannot open " << integrationInput.fileName;
					exit(1);
				}
				inStream.rdbuf(inFile.rdbuf());
			}//set the inStream.


			integrationInput.integrand = "";

			if (inFile.is_open())
			{
				cerr << "Reading " << (integrationInput.integrandType == IntegrationInput::inputPolynomial ? "polynomial" : "linear forms" ) << " from file " << integrationInput.fileName.c_str() << endl;
				getline(inStream, integrationInput.integrand, '\n');
				inFile.close();
			} else
			{
				cout << "\nEnter an integrand type: \n"
					 << "1) p (for a polynomial)\n"
					 << "2) l (for powers of linear forms)\n"
					 << "3) d (for products of powers of linear forms)\n"
					 << " :> ";
				char pl = cin.get();
				cin.ignore(1000, '\n');
				cout << "Enter the integrand in "
					 << (Poly->homogenized ? Poly->numOfVars - 1 : Poly->numOfVars)
					 << " variables: ";
				if (pl == 'p')
					integrationInput.integrandType = IntegrationInput::inputPolynomial;
				else if (pl == 'l')
					integrationInput.integrandType = IntegrationInput::inputLinearForm;
				else if ( pl == 'd' )
					integrationInput.integrandType = IntegrationInput::inputProductLinearForm;
				else
				{
					cerr << "The character " << pl << " is not a p or l or d" << endl;
					exit(1);
				}

				//char integrandLine[150];
				getline(inStream, integrationInput.integrand, '\n');
				//integrationInput.integrand = integrandLine;
			}//user supplied polynomial in file.


			if (integrationInput.topEhrhart) {
			  valuationAnswers = computeTopEhrhart(Poly, *params, integrationInput);
			}
			else {
			  valuationAnswers = computeIntegral(Poly, *params, integrationInput);
			}
		} //else integration.
		ValuationData totalValuationTimer;
		totalValuationTimer.valuationType = PolytopeValuation::entireValuation;
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
			ans = pv.findVolume(PolytopeValuation::volumeCone);

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


	//cout << params->triangulate_time << endl;
	//cout << params->dualize_time << endl;
	valuationAnswers.printResults(cout);

	delete params;
	return valuationAnswers;
}//mainValuationDriver


void Valuation::polyhedronToCones(const IntegrationInput &intInput, Polyhedron *Poly, BarvinokParameters * params)
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

		if (intInput.volumeCone || intInput.integrateLinearFormCone || intInput.integratePolynomialAsLinearFormCone || intInput.all)
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
		else
		{
			assert(Poly->homogenized == true);
			if (Poly->dualized)
			{
				cerr << "(First dualizing back... ";
				cerr.flush();
				dualizeCones(Poly->cones, Poly->numOfVars, params);
				cerr << "done!) ";
				cerr.flush();
				Poly->dualized = false; // Adjust state
			}
		}//only need vertices


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





Valuation::IntegrationInput::IntegrationInput()
{
	integrandType = nothing;

	volumeCone = false;									//volume using the cone method.
	volumeTriangulation = false;							//volume using triangulation
	integratePolynomialAsLinearFormTriangulation = false; 	//decompose polynomial to LF, use triangulation.
	integratePolynomialAsLinearFormCone = false;			//decompose polynomila to LF, use cone method.
	integratePolynomialAsPLFTriangulation = false; 			//decompose polynomial to PLF
	integrateLinearFormTriangulation = false;				//integrate linear forms using triangulation
	integrateLinearFormCone = false;						//integrate linear forms using cone method
	integrateProductLinearFormsTriangulation = false;		//integrate product of linear forms using triangulation.
	all = false;
}

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
		if (answers[i].valuationType == PolytopeValuation::volumeCone)
			out << "Volume (using the cone decomposition method)" << endl;
		else if (answers[i].valuationType == PolytopeValuation::volumeTriangulation)
			out << "Volume (using the triangulation-determinant method)"
					<< endl;
		else if (answers[i].valuationType
				== PolytopeValuation::integrateLinearFormTriangulation)
			out << "Integration of linear forms (using the triangulation method)" << endl;
		else if (answers[i].valuationType == PolytopeValuation::integrateLinearFormCone)
			out << "Integration of linear forms (using the cone method)" << endl;
		else if ( answers[i].valuationType == PolytopeValuation::integrateProductLinearFormsTriangulation)
			out << "Integration of products of linear forms (using the triangulation method)" << endl;
		else if ( answers[i].valuationType == PolytopeValuation::integratePolynomialAsLinearFormCone)
			out << "Integration of a polynomial as linear forms (using the cone method)" << endl;
		else if ( answers[i].valuationType == PolytopeValuation::integratePolynomialAsLinearFormTriangulation)
			out << "Integration of a polynomial as linear forms (using the triangulation method)" << endl;
		else if ( answers[i].valuationType == PolytopeValuation::integratePolynomialAsPLFTriangulation)
			out << "Integration of a polynomail as products of linear forms (using the triangulation method)" << endl;
		else if (answers[i].valuationType == PolytopeValuation::entireValuation)
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


