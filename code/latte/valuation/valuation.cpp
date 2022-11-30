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

	if ( intInput.volumeCone && intInput.volumeTriangulation)
		poly2 = new Polyhedron(*poly);
	else
		poly2 = poly;

	if (intInput.volumeTriangulation)
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
	}//if triangulate. original polytope should not have changed.

	if (intInput.volumeCone)
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
	}//if lawrence. Original polytope is now dilated.

	if (intInput.volumeCone && intInput.volumeTriangulation && ans1 != ans2)
	{
		cerr << "valuation.cpp: the two volume methods are different." << endl;
		cerr << "Cone-decompose:      " << ans2 << endl;
		cerr << "Triangulation: " << ans1 << endl;
		THROW_LATTE_MSG( LattException::bug_Unknown, 1, "volume computed by both methods are different. Please send bug report" );
	}//if error.
	else if (intInput.volumeCone && intInput.volumeTriangulation)
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
		THROW_LATTE_MSG(LattException::bug_Unknown, 1, "integrand type not supported.");
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
	Polyhedron *polyCopy;//if doing more than 1 method, make a deep copy of the original polytope

	assert(intInput.integrandType == IntegrationInput::inputPolynomial);


	if (intInput.integratePolynomialAsLinearFormTriangulation)
	{
		if(intInput.integratePolynomial)
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

		if (intInput.integratePolynomial)
			delete polyCopy;
	}//if doing triangulation method.


	if (intInput.integratePolynomialAsLinearFormCone)
	{
		cerr << "Going to run the cone-decomposition integration method" << endl;

		if(intInput.integratePolynomial)
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

		if(intInput.integratePolynomial)
			delete polyCopy;
	}

	if ( intInput.integratePolynomialAsPLFTriangulation )
	{
		cerr << "Going to run the polynomial to PLF method" << endl;

		if(intInput.integratePolynomial)
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

		if(intInput.integratePolynomial)
			delete polyCopy;
	}//polynomial to PLF.


	if (intInput.integratePolynomial && (ans1 != ans2 || ans1 != ans3 ) )
	{
		cerr << "Valuation.cpp: the methods are different.\n"
				<< "triangulateion    : " << ans1 << "\n"
				<< "cone-decomposition: " << ans2 << "\n"
				<< "prod linear form  : " << ans3 << "\n"
				<< endl;
		THROW_LATTE_MSG(LattException::bug_Unknown, 1, "The integrals are different. Please send bug report.");
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
	Polyhedron *poly2 = poly;//if doing both methods, make a deep copy of the original polytope.

	assert(intInput.integrandType == IntegrationInput::inputLinearForm);

	if (intInput.integrateLinearFormCone && intInput.integrateLinearFormTriangulation)
	{
		poly2 = new Polyhedron(*poly); //copy org. polytope, because it will be dilated.
	}

	if (intInput.integrateLinearFormTriangulation)
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


	if (intInput.integrateLinearFormCone)
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

	if (intInput.integrateLinearFormTriangulation && intInput.integrateLinearFormCone && ans1 != ans2)
	{
		cerr << "computeIntegralLinearForm(): the two methods are different.\n"
				<< "triangulation: " << ans1 << "\nlawrence       " << ans2
				<< endl;
		THROW_LATTE_MSG(LattException::bug_Unknown, 1, "The integrals are different. Please send bug report");
	}//if error.
	if (intInput.integrateLinearFormTriangulation && intInput.integrateLinearFormCone)
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

void Valuation::computeTopEhrhart(Polyhedron *poly,
		BarvinokParameters &myParameters, const IntegrationInput & intInput)
{
	ValuationContainer answer;
	ValuationData timer_and_result;
	RationalNTL ans;


	if (intInput.integrandType == IntegrationInput::inputPolynomial)
	{
		TopEhrhart topEhrhart(poly, myParameters, intInput.numEhrhartCoefficients, intInput.realDilations, intInput.saveTopEhrhartPolynomial);

		monomialSum originalPolynomial;// polynomial without the updated coefficients.
		loadMonomials(originalPolynomial, intInput.integrand); //get the polynomial from the string.
		topEhrhart.computeTopEhrhartPolynomial(originalPolynomial);

		destroyMonomials(originalPolynomial);
	}//the weight is a polynomial
	else if ( intInput.integrandType == IntegrationInput::inputLinearForm)
	{
		linFormSum originalLinearForm;// polynomial without the updated coefficients.

		TopEhrhart topEhrhart(poly, myParameters, intInput.numEhrhartCoefficients, intInput.realDilations,  intInput.saveTopEhrhartPolynomial);

		loadLinForms(originalLinearForm, intInput.integrand); //get the polynomial from the string.

		topEhrhart.computeTopEhrhartPolynomial(originalLinearForm);

		destroyLinForms(originalLinearForm);
	}//the weight is a power of a linear form.
	else if (intInput.unweightedCounting == true)
	{
		TopEhrhart topEhrhart(poly, myParameters, intInput.numEhrhartCoefficients, intInput.realDilations, intInput.saveTopEhrhartPolynomial);
		topEhrhart.computeTopEhrhartPolynomial();
	}
	else
	{
		THROW_LATTE_MSG(LattException::bug_NotImplementedHere, 1, "integrand type not supported");
	}

}//computeTopEhrhart


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
			removeFiles[127];
	char printLawrence[127];
	bool approx;
	bool ehrhart_polynomial, ehrhart_series, ehrhart_taylor;
	bool triangulation_specified = false;
	bool useStokes = false;
	double sampling_factor = 1.0;
	long int num_samples = -1;

	//options used by top ehrhart
	bool interactiveLatte = false; //only in the case of top ehrhart do we look at this value.

	ReadPolyhedronData read_polyhedron_data;

	struct BarvinokParameters *params = new BarvinokParameters;

	latte_banner(cerr);

	if (argc < 2)
	{
		usage(argv[0]);
		THROW_LATTE( LattException::ue_BadCommandLineOption, 0);
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
	approx = false;
	ehrhart_polynomial = false;
	params->substitution = BarvinokParameters::PolynomialSubstitution;
	//params->decomposition = BarvinokParameters::DualDecomposition;
	params->decomposition
			= BarvinokParameters::IrrationalAllPrimalDecomposition;
	// params->triangulation = BarvinokParameters::RegularTriangulationWithCdd;
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
		else if (strcmp(argv[i], "--polynomial") == 0)
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
			     << "  --valuation=volume                       Computes the volume\n"
			     << "  --valuation=integrate                    Computes an integral\n"
	       		 << "  --valuation=top-ehrhart                  Computes the top weighted Ehrhart coefficients\n"
					<< "volume and integration and options:\n"
					<< "  --cone-decompose                         Computes the valuation using tangent cones\n"
					<< "  --triangulate                            Computes the valuation over a triangulation.\n"
					<< "  --all                                    Computes the using all the two methods above.\n"
					<< "\n"
					<< "integration and weighted summation options:\n"
					<< "  --monomials=FILE                         Looks at the first line of file for a polynomial\n"
					<< "                                           encoded in maple-syntax: [ [coef, [exponent vector]], ...]\n"
					<< "  --linear-forms=FILE                      Looks at the first line of file for a sum of powers of linear forms\n"
					<< "                                           encoded in maple-syntax: [ [coef, [power, [linear form]]], ...]\n"
					<< "  --product-linear-forms=FILE              Looks at the first line of file for a product of linear forms\n"
					<< "                                           in maple-syntax: [ [coef, [[power, [linear form]],\n"
					<< "                                           [power, [linear form]], ...]],  ...] \n"
					/// FIXME: Add further options!
					<< "                                           If cannot open file, the line is read from std in.\n"
					<< "top Ehrhart options:\n"
					<< "  --num-coefficients=K                     Number of highest Ehrhart coefficients to compute\n"
					<< "                                           (default: compute all, incrementally)\n"
					<< "  --real-dilations                         Output formulas valid for real dilations also\n"
					<< "  --top-ehrhart-save=FILE                  Writes the Ehrhart polynomial to the file.\n"
					<< "  --top-ehrhart-unweighted                 The weight function is assumed to be 1 (default)\n"
					<< "  --interactive-mode                       Print the screen asking for an integrand type,\n"
					<< "                                           Only polynomials and linear forms are valid.\n"
					<< "Example: " << argv[0]
					<< " --valuation=volume --cone-decompose file.latte\n"
					<< "         (will compute the volume using the cone-decomposition method over the polytope in file.latte.)\n"
					<< "Example: " << argv[0]
					<< " --valuation=integrate --monomials=poly.txt file.latte\n"
					<< "         (will compute the integral of the polynomial in poly.txt over the polytope in file.latte.)\n"
					<< "Example: " << argv[0]
					<< " --valuation=top-ehrhart --num-coefficients=4 --real-dilations --monomials=poly.txt file.latte\n"
					<< "         (will compute the top 4 Ehrhart coefficients with a polynomial weight in poly.txt over the polytope in file.latte.)\n"

					<< endl;
			exit(0);
		} else if ( strncmp(argv[i], "--valuation=", 11) == 0)
		{
			if ( strcmp(argv[i], "--valuation=volume") == 0)
				integrationInput.valuationVolume = true;
			else if ( strcmp(argv[i], "--valuation=integrate") == 0)
				integrationInput.valuationIntegrate = true;
			else if ( strcmp(argv[i], "--valuation=top-ehrhart") == 0)
			{
				integrationInput.valuationEhrhart = true;
				integrationInput.unweightedCounting = true;
			}
			else
				; //nothing.
		}
		else if (strcmp(argv[i], "--top-ehrhart-unweighted")==0)
		{
			integrationInput.unweightedCounting = true;
			interactiveLatte = false;
		}
		else if (strcmp(argv[i], "--interactive-mode") == 0)
		{
			interactiveLatte = true; //note, when in
		}
		else if (strncmp(argv[i], "--num-coefficients=", 19) == 0) {
			integrationInput.numEhrhartCoefficients = atoi(argv[i] + 19);
		}
		else if (strncmp(argv[i], "--top-ehrhart-save=", 19) == 0)
		{
			integrationInput.saveTopEhrhartPolynomial = (argv[i] + 19);
		}
		else if (strncmp(argv[i], "--real-dilations", 6) == 0) {
			integrationInput.realDilations = true;
		}
		 else if (strcmp(argv[i], "--cone-decompose") == 0)
		{
			integrationInput.useTangentCones = true;
		}
		else if (strcmp(argv[i], "--triangulate") == 0)
		{
			integrationInput.useTriangulation = true;
		}
		else if ( strcmp(argv[i], "--polynomial-as-plf") == 0)
		{
			integrationInput.polynomialAsPLF = true;
		}
/*
//		 else if (strcmp(argv[i], "--all") == 0)
		{
			strcpy(valuationAlg, "all");
			strcpy(read_polyhedron_data.dualApproach, "no");
		}
		else if (strcmp(argv[i], "--print-lawrence-function") == 0 || strcmp(argv[i], "--print-cone-decompose-function") == 0)
			strcpy(printLawrence, "yes");
		else if ( strcmp(argv[i], "--stokes") == 0)
			useStokes = true;
//
*/		else if (strncmp(argv[i], "--monomials=", 12) == 0)
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
			THROW_LATTE(LattException::ue_UnknownCommandLineOption, 0);
		}
	}

	if ( ! integrationInput.valuationEhrhart
			&& !integrationInput.valuationIntegrate
			&& !integrationInput.valuationVolume
			)
	{
		cerr << "--valuation=??? missing. Type --help" << endl;
		THROW_LATTE(LattException::ue_BadCommandLineOption, 0);
	}


	integrationInput.processUserInput(); //sets which algorithms will be used (if the user gave them on the command line)

	if ( integrationInput.volumeCone
			|| integrationInput.integrateLinearFormCone
			|| integrationInput.integratePolynomialAsLinearFormCone
			|| integrationInput.topEhrhart
			|| integrationInput.useTangentCones
			|| (integrationInput.useTangentCones == false && integrationInput.useTriangulation == false)
		)
	{
		strcpy(read_polyhedron_data.dualApproach, "no");
	}//compute tangent cones.
	else
	{
		strcpy(read_polyhedron_data.dualApproach, "yes");
	}//compute vertices.

	if (read_polyhedron_data.expect_filename)
	{
		cerr << "Filename missing" << endl;
		THROW_LATTE( LattException::ue_FileNameMissing,0);
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



	//the stokes recursive method is not finished and experimental.
	// Users should always enter this if statement.
	if ( useStokes == false)
	{
		Polyhedron *Poly = read_polyhedron_data.read_polyhedron(params);

		params->Flags = flags;
		params->Number_of_Variables = Poly->numOfVars;
		params->max_determinant = 1;
		params->File_Name = (char*) fileName;

		//poly is updated.
		polyhedronToCones(integrationInput, Poly, params);
		//now the cones of poly are the tangent cones or the lifted cone of the vertices.

		if (integrationInput.integrandType == IntegrationInput::inputVolume)
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
					THROW_LATTE( LattException::fe_Open, 0);
				}
				inStream.rdbuf(inFile.rdbuf());
			}//set the inStream.


			integrationInput.integrand = "";

			//get the integrand from the file or the user.
			if (inFile.is_open())
			{
				cerr << "Reading " << (integrationInput.integrandType == IntegrationInput::inputPolynomial ? "polynomial" : "linear forms" ) << " from file " << integrationInput.fileName.c_str() << endl;
				getline(inStream, integrationInput.integrand, '\n');
				inFile.close();
			} else if (integrationInput.unweightedCounting == true && interactiveLatte == false)
			{
				; //don't ask the user for the integrand.
			}
			else
			{
				cout << "\nEnter an integrand type: \n"
					 << "1) p (for a polynomial)\n"
					 << "2) l (for powers of linear forms)\n"
					 << "3) d (for products of powers of linear forms)\n"
					 << " :> ";
				char pl = cin.get();
				cin.ignore(1000, '\n');
				cout << "Enter the integrand using "
					 << (Poly->homogenized ? Poly->numOfVars - 1 : Poly->numOfVars)
					 << " variables\n >";
				if (pl == 'p')
				{
					integrationInput.integrandType = IntegrationInput::inputPolynomial;
				}
				else if (pl == 'l')
				{
					integrationInput.integrandType = IntegrationInput::inputLinearForm;
				}
				else if ( pl == 'd' )
					integrationInput.integrandType = IntegrationInput::inputProductLinearForm;
				else
				{
					cerr << "The character " << pl << " is not a p or l or d" << endl;
					THROW_LATTE( LattException::ue_BadCommandLineOption, 0);
				}

				integrationInput.unweightedCounting = false;

				integrationInput.processUserInput();
				//char integrandLine[150];
				getline(inStream, integrationInput.integrand, '\n');
				//integrationInput.integrand = integrandLine;
			}//user supplied integrand in a file.


			if (integrationInput.topEhrhart) {
				computeTopEhrhart(Poly, *params, integrationInput);
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
	{ //use Stokes
		//NOTE:: THIS RECURSIVE METHOD IS NOT DONE.
		//THIS SHOULD NOT WORK.
		THROW_LATTE_MSG( LattException::bug_NotImplementedHere, 1, "The experimental Stokes code has been removed" );

	}//else use stokes.




	if (read_polyhedron_data.rationalCone[0] == 'y')
	{
	  string new_name = string(argv[argc - 1]) + ".rat";
	  cerr << endl << "Rational function written to `" << new_name << "'" << endl << endl;
	  rename_with_error_check("simplify.sum", new_name);
	}

	if (printfile[0] == 'y')
	{
	  string new_name = string(argv[argc - 1]) + ".rat";
	  cerr << endl << "Rational function written to `" << new_name << "'" << endl << endl;
	  rename_with_error_check("func.rat", new_name);
	}
	if ((removeFiles[0] == 'y')
			&& (read_polyhedron_data.dualApproach[0] == 'n'))
	{

		// strcpy(command, "rm -f ");
		// strcat(command, fileName);
		// strcat(command, ".ext");
		// system_with_error_check(command);

		// strcpy(command, "rm -f ");
		// strcat(command, fileName);
		// strcat(command, ".cdd");
		// system_with_error_check(command);

		// if (read_polyhedron_data.Memory_Save[0] == 'n')
		// {
		// 	strcpy(command, "rm -f ");
		// 	strcat(command, fileName);
		// 	strcat(command, ".maple");
		// 	system_with_error_check(command);
		// }

		// strcpy(command, "rm -f ");
		// strcat(command, fileName);
		// strcat(command, ".ead");
		// system_with_error_check(command);

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
	//	cout << "read.input_dualized" << read_polyhedron_data.input_dualized << endl;

		if (intInput.volumeCone
				|| intInput.integrateLinearFormCone
				|| intInput.integratePolynomialAsLinearFormCone
				|| intInput.topEhrhart
				|| (intInput.useTangentCones == false && intInput.useTriangulation == false) //the user didn't give an integrand nor integration method
				|| (intInput.useTangentCones == true) //the user gave --cone-decompose but didn't give an integrand.
				)
		{
			assert(Poly->homogenized == false);
			if (Poly->dualized)
			{
				cerr << "(First dualizing back... ";
				cerr.flush();
				dualizeCones(Poly->cones, Poly->numOfVars, params);
				cerr << "done.) ";
				cerr.flush();
				Poly->dualized = false;
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

	//volume algorithms
	volumeCone = false;									//volume using the cone method.
	volumeTriangulation = false;							//volume using triangulation

	//integration algorithms.
	integratePolynomial = false;
	integratePolynomialAsLinearFormTriangulation = false; 	//decompose polynomial to LF, use triangulation.
	integratePolynomialAsLinearFormCone = false;			//decompose polynomial to LF, use cone method.
	integratePolynomialAsPLFTriangulation = false; 		//decompost polynomial to PLF, use triangulation.
	integrateLinearFormTriangulation = false;				//integrate linear forms using triangulation
	integrateLinearFormCone = false;						//integrate linear forms using cone method
	integrateProductLinearFormsTriangulation = false;		//integrate product of linear forms using triangulation.

	//Ehrhart algorithms.
	topEhrhart = false;					//compute top Ehrhart coefficients using cone method only.
	numEhrhartCoefficients = -1; // incremental mode
	realDilations = false; //formula is only valid for integer dilations, not rational/real ones.
	saveTopEhrhartPolynomial = "-1"; //default is to not save the polynomial.

	//Command line options. These are used by processUserInput()
	//  to set up which algorithms will be used.
	valuationVolume =false;			//--valuation=volume
	valuationIntegrate= false;	//etc
	valuationEhrhart= false;
	useTangentCones= false;		//--cone-decompose
	useTriangulation= false;	//--triangulate
	polynomialAsPLF= false;		//--polynomial-as-plf
}


void Valuation::IntegrationInput::processUserInput()
{
	if ( valuationVolume)
	{
		integrandType = inputVolume; //the integrand is 1
		if ( useTangentCones ) volumeCone = true;
		else if ( useTriangulation) volumeTriangulation = true;
		else
		{
			volumeCone = true;
			volumeTriangulation = true;
		}//default, do both.
	}//if volume
	else if ( valuationIntegrate)
	{
		if (integrandType == inputLinearForm)
		{
			if ( useTangentCones ) integrateLinearFormCone = true;
			else if ( useTriangulation) integrateLinearFormTriangulation = true;
			else
			{
				integrateLinearFormCone = true;
				integrateLinearFormTriangulation = true;
			}//default, do both.
		}//if linear form
		else if ( integrandType == inputProductLinearForm)
		{
			integrateProductLinearFormsTriangulation = true;
		}//if plf
		else if ( integrandType == inputPolynomial)
		{
			if ( useTangentCones ) integratePolynomialAsLinearFormCone = true;
			else if ( useTriangulation) integratePolynomialAsLinearFormTriangulation = true;
			else if ( polynomialAsPLF) integratePolynomialAsPLFTriangulation = true;
			else
			{
				integratePolynomialAsLinearFormCone = true;
				integratePolynomialAsLinearFormTriangulation = true;
				integratePolynomialAsPLFTriangulation = true;
				integratePolynomial = true; //do all of them.
			}//default, do all.
		}//if polynomial.
	}//else integration
	else if (valuationEhrhart)
	{
		topEhrhart = true;
	}
}//processUserInput

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
	for (size_t i = 0; i < answers.size(); ++i)
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


