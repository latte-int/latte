/*
 * nonlinearOptimizationDriver.cpp
 *
 *  Created on: Apr 23, 2014
 *      Author: bedutra
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <climits>
#include <sstream>
#include <getopt.h>


#include "banner.h"
#include "barvinok/barvinok.h"
#include "ReadPolyhedron.h"
#include "LattException.h"
#include "nonlinearOptimization/WeightedExponentialSubs.h"
#include "nonlinearOptimization/BoxOptimization.h"
#include "nonlinearOptimization/BoxOptimizationContinuous.h"
#include "print.h"
#include "dual.h"
#include "rational.h"
#include "integration/burstTrie.h"
#include "integration/PolyTrie.h"

#include "MpqClassLazy.h"

using namespace std;


void printHelpMenu()
{


	cout <<
	"Input files:\n"
	"  --boxFile FILENAME{.box}               Input knapsack file.\n"
	"  --monomials FILENAME{.poly}            Input polynomial file\n"
	"  --linear-forms FILENAME{.lf}           Input powers of linear forms polynomial (used with --count)\n"
    "Options that control what to compute:\n"
    "  --opt-lf                               Optimization via linear forms\n"
	"  --opt-lfc                              Optimization via linear forms with cache\n"
	"  --opt-ns                               Optimization via natural summation\n"
	"  --count                                Finds the weighted average of f(x)\n"
	"  --opt-cont-ns                          Finds lower bound for continuous optimization\n"
    "Other options:\n"
    "  --k, -k INT                            Sets (f(x)+s)^k\n"
    "  --epi                                  If --opt, sets k\n"
    "  --help, -h                             Prints this help message\n"
    "\nExamples:\n"
    "  ./boxOpt --boxFile box20d.box --monomials deg5.poly --opt -k 5\n"
    << endl;

}






void run()
{
	vector<MpqLazy> v;
	MpqLazy a;
	a = 34;

	v.push_back(a);
	cout << "ok" <<endl;
}
int main(int argc, char *argv[]) {


	run();
	//return 0;
	string linFormFileName, polynomialFileName, boxFileName;
	string cmd;
	int userK = -1;
	int k;
	RR epsilon;
	Timer totalTime("Total Time");

	epsilon = 0.10;

	latte_banner(cerr);

	cerr << "Invocation: ";
	for (int i = 0; i < argc; i++) {
		cerr << argv[i] << " ";
	}
	cerr << endl;

	//process each option.
	while (1)
	{
		static struct option long_options[] =
		{
		{ "linear-forms",	  			  required_argument, 0, 0x100 },
		{ "monomials",  				  required_argument, 0, 0x101 },
		{ "boxFile",					  required_argument, 0, 0x102 }, //hex value larger than 1 byte/char
		{ "help",						  no_argument,       0, 'h' },
		{ "k",							  required_argument, 0, 'k' },
		{ "epi",          				  required_argument, 0, 0x103},
		{ "range",						  no_argument, 0, 0x104},
		{ "opt-lf",						  no_argument, 0, 0x105},
		{ "opt-lfc",					  no_argument, 0, 0x106},
		{ "opt-ns",						  no_argument, 0, 0x107},
		{ "opt-ratio",					  no_argument, 0, 0x108},
		{ "count",						  no_argument, 0, 0x109 },
		{ "opt-cont-ns",				  no_argument, 0, 0x10A},
		{ 0, 0, 0, 0 } };
		/* getopt_long stores the option index here. */

		int option_index = 0;
		int c;

		//single-character=short option-name. x: means x takes a required argument, x:: means x takes an optional argument.
		c = getopt_long(argc, argv, "hk:", long_options, &option_index);

		if (c == -1)
			break;

		switch (c)
		{
			case 0:
				// If this option set a flag, do nothing
				break;
			case 0x100:
				linFormFileName = optarg;
				break;
			case 0x101:
				polynomialFileName = optarg;
				break;
			case 0x102:
				boxFileName = optarg;
				break;
			case 'h':
				printHelpMenu();
				exit(0);
				break;
			case 'k':
				userK  = atoi(optarg);
				break;
			case 0x103:
				epsilon = atof(optarg);
				break;
			case 0x104:
				cmd = "range";
				break;
			case 0x105:
				cmd = "opt-lf";
				break;
			case 0x106:
				cmd = "opt-lfc";
				break;
			case 0x107:
				cmd = "opt-ns";
				break;
			case 0x108:
				cmd = "opt-ratio";
				break;
			case 0x109:
				cmd = "count";
				break;
			case 0x10A:
				cmd = "opt-cont-ns";
				break;
			default:
				cerr << "Unknown command/option " << endl;
				THROW_LATTE(LattException::ue_BadCommandLineOption);
				exit(1);
		}
	}//while.


	if (boxFileName.length() == 0) {
		cerr << "box file is missing" << endl;
		THROW_LATTE(LattException::ue_FileNameMissing);
	}

	linFormSum originalLinearForm;
	if ( linFormFileName.length())
	{
		ifstream linFormFile(linFormFileName.c_str());
		string linFormStr;
		getline(linFormFile, linFormStr);
		cout << "lin form str: " << linFormStr.c_str() << endl;
		loadLinForms(originalLinearForm, linFormStr.c_str()); //incorrect. need to mult or maybe divide by M!
		linFormFile.close();
	}

	monomialSum originalPolynomial;
	if ( polynomialFileName.length())
	{
		ifstream polynomialFile(polynomialFileName.c_str());
		string polyStr;
		getline(polynomialFile, polyStr);
		cout << "poly str: " << polyStr.c_str() << endl;
		loadMonomials(originalPolynomial, polyStr.c_str());
		polynomialFile.close();
	}


	ifstream boxFile(boxFileName.c_str());
	int dim;
	boxFile >> dim;

	vec_RR lowerBound, upperBound;
	lowerBound.SetLength(dim);
	upperBound.SetLength(dim);
	for(int i = 0; i < dim; ++i)
		boxFile >> lowerBound[i] >> upperBound[i];

	cout << "lb: " << lowerBound << endl;
	cout << "ub: " << upperBound << endl;


	RR N;
	//cout << "epsilon=" << epsilon << endl;
	N = 1;
	for (int i = 0;  i < lowerBound.length(); ++i)
		N *= (upperBound[i] - lowerBound[i] + 1);
	N = ceil((1.0 + inv(epsilon))* log(N));
	k = INT_MAX;
	if ( N < INT_MAX)
		k = to_int(N);
	else
		cout << "Warning: k is larger than " << INT_MAX << endl;



	if (userK > 0)
		k = userK;

	cout << "k is " << k << endl;

	totalTime.start();
	if ( cmd == "count")
	{
		if ( polynomialFileName.length())
		{
			BTrieIterator<RationalNTL, int>* polynomialItr = new BTrieIterator<RationalNTL, int> ();
			originalLinearForm.termCount = 0;
			originalLinearForm.varCount = originalPolynomial.varCount;
			polynomialItr->setTrie(originalPolynomial.myMonomials, originalPolynomial.varCount);
			decompose(polynomialItr, originalLinearForm);
			destroyMonomials(originalPolynomial);
		}//decompose polynomial into power of linear forms

		mpq_class weightedCount = computeWeightedCountingBox(conv<vec_ZZ>(lowerBound), conv<vec_ZZ>(upperBound), originalLinearForm);
		cout << "Final count: " << weightedCount << endl;
		destroyLinForms(originalLinearForm);
	}


	if ( cmd == "opt-lf")
	{
		BoxOptimization bo;
		bo.setPolynomial(originalPolynomial);
		bo.setBounds(conv<vec_ZZ>(lowerBound), conv<vec_ZZ>(upperBound));

		if ( bo.isTrivial())
		{
			bo.enumerateProblem(conv<vec_ZZ>(lowerBound), conv<vec_ZZ>(upperBound), originalPolynomial);
			cout << "optimal: " << bo.L << " <= f(x) <= " << bo.U << endl;
			return 0;
		}

		Timer timeSummation("Total Time for lf decomp + counting");
		timeSummation.start();
		bo.setPower(k);
		bo.decomposePoly(BoxOptimization::lf);
		bo.findSPolynomial(BoxOptimization::lf, conv<vec_ZZ>(lowerBound), conv<vec_ZZ>(upperBound));
		timeSummation.stop();
		cout << timeSummation << endl;
		bo.printSpolynomial();

		cout << "old range " << bo.L << " <=f(x)<= " << bo.U << " old max bounds: " << bo.maximumLowerBound() << " <=max(f)<= " << bo.maximumUpperbound() << " gap %:" << (bo.maximumUpperbound() - bo.maximumLowerBound())/bo.maximumUpperbound() << endl;
		bo.findRange(10);
		cout << "new range " << bo.L << " <=f(x)<= " << bo.U << " new max bounds: " << bo.maximumLowerBound() << " <=max(f)<= " << bo.maximumUpperbound() << " gap %:" << (bo.maximumUpperbound() - bo.maximumLowerBound())/bo.maximumUpperbound() << endl;

	}


	if ( cmd == "opt-lfc")
	{
		BoxOptimization bo;
		bo.setPolynomial(originalPolynomial);
		bo.setBounds(conv<vec_ZZ>(lowerBound), conv<vec_ZZ>(upperBound));

		if ( bo.isTrivial())
		{
			bo.enumerateProblem(conv<vec_ZZ>(lowerBound), conv<vec_ZZ>(upperBound), originalPolynomial);
			cout << "optimal: " << bo.L << " <= f(x) <= " << bo.U << endl;
			return 0;
		}

		Timer timeSummation("Total Time for decomp + making the table ");
		Timer timeSpoly("Total Time for finding s-poly ");
		timeSummation.start();
		bo.setPower(k);
		bo.decomposePoly(BoxOptimization::lfCache);
		timeSummation.stop();


		for(int i = 0; i <1000; ++i)
		{

			for(int a = 0; a < 5; ++a)
			if ( rand() % 2)
			{
				int j = rand() % lowerBound.length();
				lowerBound[j] = to_int((lowerBound[j] + upperBound[j])/2);
			}
			else
			{
				int j = rand() % lowerBound.length();
				upperBound[j] = to_int((lowerBound[j] + upperBound[j])/2);
			}

			cout << "\nLow: " << lowerBound << " up: " << upperBound << endl;

			timeSpoly.start();
			bo.setBounds(conv<vec_ZZ>(lowerBound), conv<vec_ZZ>(upperBound));
			bo.findSPolynomial(BoxOptimization::lfCache,conv<vec_ZZ>(lowerBound), conv<vec_ZZ>(upperBound));
			timeSpoly.stop();
			bo.printSpolynomial();
			cout << "WARNING: lfc is broken for rational polynomials. Remove this code in production" << endl;
			THROW_LATTE(LattException::bug_Unknown);


			cout << "old range " << bo.L << " <=f(x)<= " << bo.U << " old max bounds: " << bo.maximumLowerBound() << " <=max(f)<= " << bo.maximumUpperbound() << " gap %:" << (bo.maximumUpperbound() - bo.maximumLowerBound())/bo.maximumUpperbound() << endl;
			bo.findRange(10);
			cout << "new range " << bo.L << " <=f(x)<= " << bo.U << " new max bounds: " << bo.maximumLowerBound() << " <=max(f)<= " << bo.maximumUpperbound() << " gap %:" << (bo.maximumUpperbound() - bo.maximumLowerBound())/bo.maximumUpperbound() << endl;
			cout << "i=" << i << endl;

			if ( bo.N <= 50)
				break;

		}
		cout << timeSummation << endl;
		cout << timeSpoly << endl;
	}


	if ( cmd == "opt-ns")
	{

		vec_ZZ zzLowerBound(conv<vec_ZZ>(lowerBound));
		vec_ZZ zzUpperBound(conv<vec_ZZ>(upperBound));

		BoxOptimization bo;
		bo.setPolynomial(originalPolynomial);
		bo.setBounds(zzLowerBound, zzUpperBound);

		/*
		if ( bo.isTrivial())
		{
			bo.enumerateProblem(lowerBound, upperBound, originalPolynomial);
			cout << "optimal: " << bo.L << " <= f(x) <= " << bo.U << endl;
			return 0;
		}
		*/



		bo.setPower(k);
		bo.decomposePoly(BoxOptimization::naturalSummation);
		bo.findSPolynomial(BoxOptimization::naturalSummation, conv<vec_ZZ>(lowerBound), conv<vec_ZZ>(upperBound));
		bo.printSpolynomial();

		cout << "old range " << bo.L << " <=f(x)<= " << bo.U << " old max bounds: " << bo.maximumLowerBound() << " <=max(f)<= " << bo.maximumUpperbound() << " gap %:" << (bo.maximumUpperbound() - bo.maximumLowerBound())/bo.maximumUpperbound() << endl;
		bo.findRange(10);
		cout << "new range " << bo.L << " <=f(x)<= " << bo.U << " new max bounds: " << bo.maximumLowerBound() << " <=max(f)<= " << bo.maximumUpperbound() << " gap %:" << (bo.maximumUpperbound() - bo.maximumLowerBound())/bo.maximumUpperbound() << endl;

	}


	if ( cmd == "opt-ratio")
	{
		BoxOptimization bo1, bo2;
		bo1.setPolynomial(originalPolynomial);
		bo1.setBounds(conv<vec_ZZ>(lowerBound), conv<vec_ZZ>(upperBound));

		bo2.setPolynomial(originalPolynomial);
		bo2.setBounds(conv<vec_ZZ>(lowerBound), conv<vec_ZZ>(upperBound));

		//todo: if trivial?

		bo1.setPower(k);
		bo1.decomposePoly(BoxOptimization::naturalSummation);
		bo1.findSPolynomial(BoxOptimization::naturalSummation, conv<vec_ZZ>(lowerBound), conv<vec_ZZ>(upperBound));


		bo2.setPower(k+1);
		bo2.decomposePoly(BoxOptimization::naturalSummation);
		bo2.findSPolynomial(BoxOptimization::naturalSummation, conv<vec_ZZ>(lowerBound), conv<vec_ZZ>(upperBound));




		bo1.printSpolynomial();
		bo2.printSpolynomial();



		cout << "old range " << bo1.L << " <=f(x)<= " << bo1.U << endl;;
		bo1.findRange(10);
		bo2.findRange(10);

		cout << "new range k   " << bo1.L << " <=f(x)<= " << bo1.U << " new max bounds: " << bo1.maximumLowerBound() << " <=max(f)<= " << bo1.maximumUpperbound() << " gap %:" << (bo1.maximumUpperbound() - bo1.maximumLowerBound())/bo1.maximumUpperbound() << endl;
		cout << "new range k+1 " << bo2.L << " <=f(x)<= " << bo2.U << " new max bounds: " << bo2.maximumLowerBound() << " <=max(f)<= " << bo2.maximumUpperbound() << " gap %:" << (bo2.maximumUpperbound() - bo2.maximumLowerBound())/bo2.maximumUpperbound() << endl;

		cout << bo2.L + bo2.currentMap.eval(-bo2.L) / bo1.currentMap.eval(-bo2.L) << " max(f)" << endl;

	}

	if (cmd == "opt-cont-ns")
	{

		ifstream bf(boxFileName.c_str());
		int dim;
		bf >> dim;

		vec_RR rrLowerBound, rrUpperBound;
		rrLowerBound.SetLength(dim);
		rrUpperBound.SetLength(dim);
		for(int i = 0; i < dim; ++i)
			bf >> rrLowerBound[i] >> rrUpperBound[i];
		bf.close();

		RR::SetPrecision(2000);
		cout << "current precision " << RR::precision() << endl;

		BoxOptimizationContinuous bo1;
		bo1.setPolynomial(originalPolynomial);



		//vec_RationalNTL rrLowerBound(lowerBound);
		//vec_RationalNTL rrUpperBound(upperBound);


		bo1.setPower(k);
		bo1.setBounds(rrLowerBound, rrUpperBound);
		bo1.findSPolynomial(BoxOptimizationContinuous::naturalSummation, rrLowerBound, rrUpperBound);
		bo1.printSpolynomial();

		BoxOptimizationContinuous bo2;
		bo2 = bo1;
		bo2.setPower(k+1);
		bo2.setBounds(rrLowerBound, rrUpperBound);
		bo2.findSPolynomial(BoxOptimizationContinuous::naturalSummation, rrLowerBound, rrUpperBound);
		bo2.printSpolynomial();
		//bo2.findRange(10);
		cout << "Range helped: " << bo2.findRange(10) << endl;



		cout << "results " << endl;
		cout << " range " << bo1.L << " <=f(x)<= " << bo1.U << endl;
		cout << " new max bounds: " << bo1.maximumLowerBound() << " <=max(f)<= " << bo1.maximumUpperbound() << endl; //" gap %:" << (bo1.maximumUpperbound() - bo1.maximumLowerBound())/bo1.maximumUpperbound() << endl;
        cout << " new ratio lower bound: " << bo2.L + bo2.evalSpoly(-bo2.L)/ bo1.evalSpoly(-bo2.L) << endl;
        cout << "                       bo2.L=" << bo2.L << endl;
        cout << "                       \\int f^{k+1}=" << bo2.evalSpoly(-bo2.L)<< endl;
        cout << "                       \\int f^k    =" << bo1.evalSpoly(-bo2.L)<< endl;


        bo2.printStats();
	}
	totalTime.stop();
	cout << totalTime << endl;



	return 0;
}

