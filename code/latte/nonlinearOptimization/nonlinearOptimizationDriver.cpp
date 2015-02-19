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
    "  --opt                                  Does optimization\n"
    "  --range                                Finds range of function\n"
	"  --count                                Finds the weighted average of f(x)\n"
    "Other options:\n"
    "  --k, -k INT                            Sets (f(x)+s)^k\n"
    "  --epi                                  If --opt, sets k\n"
    "  --help, -h                             Prints this help message\n"
    "\nExamples:\n"
    "  ./boxOpt --boxFile box20d.box --monomials deg5.poly --opt -k 5\n"
    << endl;

}


RationalNTL evaluate(monomialSum & poly, const vec_ZZ & point)
{
	RationalNTL ans;
	BTrieIterator<RationalNTL, int>* pItr =	new BTrieIterator<RationalNTL, int> ();
	pItr->setTrie(poly.myMonomials,	poly.varCount);
	pItr->begin();


	term<RationalNTL, int>* term;
	int * exp;
	exp = new int[poly.varCount];

	for (term = pItr->nextTerm(); term; term = pItr->nextTerm())
	{
		RationalNTL value;
		value = 1;
		for (int currentPower = 0; currentPower < poly.varCount; ++currentPower)
		{
			value *= power(point[currentPower], term->exps[currentPower]);
			value *= term->coef;
		}
		ans += value;
	}
	return ans;
	delete pItr;
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
		{ "opt",						  no_argument, 0, 0x105},
		{ "count",						  no_argument, 0, 0x106 },
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
				cmd = "opt";
				break;
			case 0x106:
				cmd = "count";
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

	vec_ZZ lowerBound, upperBound;
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
		N *= to_RR(upperBound[i] - lowerBound[i] + 1);
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

		mpq_class weightedCount = computeWeightedCountingBox(lowerBound, upperBound, originalLinearForm);
		cout << "Final count: " << weightedCount << endl;
		destroyLinForms(originalLinearForm);
	}
/*
	if ( cmd == "range" )
	{
		BoxOptimization bo;
		bo.setPolynomial(originalPolynomial);
		bo.setBounds(lowerBound, upperBound);
		if ( bo.isTrivial())
		{
			bo.enumerateProblem(lowerBound, upperBound, originalPolynomial);
			cout << "optimal: " << bo.L << " <= f(x) <= " << bo.U << endl;
			return 0;
		}

		bo.setPower(k);
		bo.findRange(10);

	}
*/
	if ( cmd == "opt")
	{
		BoxOptimization bo;
		bo.setPolynomial(originalPolynomial);
		bo.setBounds(lowerBound, upperBound);

		if ( bo.isTrivial())
		{
			bo.enumerateProblem(lowerBound, upperBound, originalPolynomial);
			cout << "optimal: " << bo.L << " <= f(x) <= " << bo.U << endl;
			return 0;
		}

		Timer timeSummation("Total Time for decomp + making the table ");
		Timer timeSpoly("Total Time for finding s-poly ");
		timeSummation.start();
		bo.setPower(k, false);
		timeSummation.stop();


		for(int i = 0; i <1000; ++i)
		{

			for(int a = 0; a < 5; ++a)
			if ( rand() % 2)
			{
				int j = rand() % lowerBound.length();
				lowerBound[j] = (lowerBound[j] + upperBound[j])/2;
			}
			else
			{
				int j = rand() % lowerBound.length();
				upperBound[j] = (lowerBound[j] + upperBound[j])/2;
			}

			/*
			for(int j = 0; j < lowerBound.length(); ++j)
			{
				lowerBound[j] = (rand() % 10) - 5;
				upperBound[j] = lowerBound[j] + (rand() % 5);
			}
			*/


			cout << "\nLow: " << lowerBound << " up: " << upperBound << endl;

			timeSpoly.start();
			bo.setBounds(lowerBound, upperBound);
			bo.findSPolynomial(lowerBound, upperBound);
			timeSpoly.stop();


			cout << "old range " << bo.L << " <=f(x)<= " << bo.U << " old max bounds: " << bo.maximumLowerBound() << " <=max(f)<= " << bo.maximumUpperbound() << " gap %:" << (bo.maximumUpperbound() - bo.maximumLowerBound())/bo.maximumUpperbound() << endl;
			bo.findRange(10);
			cout << "new range " << bo.L << " <=f(x)<= " << bo.U << " new max bounds: " << bo.maximumLowerBound() << " <=max(f)<= " << bo.maximumUpperbound() << " gap %:" << (bo.maximumUpperbound() - bo.maximumLowerBound())/bo.maximumUpperbound() << endl;
			cout << "i=" << i << endl;

			if ( bo.N <= 10)
				break;

		}
		cout << timeSummation << endl;
		cout << timeSpoly << endl;
	}
	totalTime.stop();
	cout << totalTime << endl;

/*
	if( cmd == "loop")
	{
		ZZ i;
		for(i = lowerBound[0]; i <= upperBound[0]; ++i)
		{
			vec_ZZ current;
			current = lowerBound;
			current[0] = i;
			RationalNTL value;
			value = evaluate(originalPolynomial, current);
		}
		totalTime.stop();
		i = 1;
		for(int j = 1; j < 20; ++j)
			i *= (upperBound[0] - lowerBound[0] + 1);
		cout << "Total time  = " << totalTime.get_seconds() << endl;
		cout << "Total time for dim 20 = " << to_RR(i)*totalTime.get_seconds() << endl;

		exit(1);
	}

	totalTime.stop();
	cout  << totalTime << endl;
	cout << "Total time for dim 20 = " << 20*totalTime.get_seconds() << endl;
*/

	return 0;
}

