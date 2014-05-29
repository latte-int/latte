/*
 * nonlinearOptimizationDriver.cpp
 *
 *  Created on: Apr 23, 2014
 *      Author: bedutra
 */

#include <iostream>
#include <fstream>
#include <string>
#include <climits>

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

using namespace std;

int main2(int argc, const char *argv[]) ;
int main1(int argc, const char *argv[]) ;

int main1(int argc, const char *argv[]) {
	if (argv[1][0] == '1')
		main1(argc - 1, argv + 1);
	else
		main2(argc -1, argv + 1);

	return 0;
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

}

int main(int argc, const char *argv[]) {

	string linFormFileName, polynomialFileName, boxFileName;
	string cmd;
	int userK = -1;
	RR epsilon;
	Timer totalTime("Total Time");

	epsilon = 0.10;

	latte_banner(cerr);

	cerr << "Invocation: ";
	for (int i = 0; i < argc; i++) {
		cerr << argv[i] << " ";
	}
	cerr << endl;

	for (int i = 1; i < argc; i++) {
		if (strncmp(argv[i], "--linear-forms=", 15) == 0){
			linFormFileName = string(argv[i] + 15);
		} else if (strncmp(argv[i], "--monomials=", 12) == 0){
			polynomialFileName = string(argv[i] + 12);
		} else if (strncmp(argv[i], "--boxFile=", 10) == 0) {
			boxFileName = string(argv[i] + 10);
		} else if (strncmp(argv[i], "--count", 10) == 0) {
			cmd = "count";
		} else if (strcmp(argv[i], "--opt") == 0){
			cmd = "opt";
		} else if (strcmp(argv[i], "--range") == 0) {
			cmd = "range";
		} else if (strncmp(argv[i], "--k=",4) == 0) {
			userK = atoi(argv[i]+4);
		} else if ( strncmp(argv[i], "--epi=", 6) == 0) {
			epsilon = to_RR(argv[i]+ 6);
		} else if ( strncmp(argv[i], "--loop", 6) == 0) {
			cmd = "loop";
		} else if (strcmp(argv[i], "--help") == 0)	{
			cout << "--linear-forms=FILE\n"
				 << "--monomials=FILE\n"
				 << "--boxFile=FILE\n"
				 << "--count\n"
				 << "--range\n"
				 << "--opt\n";
			exit(0);
		} else {
			cerr << "Unknown command/option " << argv[i] << endl;
			THROW_LATTE_MSG(LattException::ue_BadCommandLineOption, argv[i]);
		}
	} //for i.

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
		loadLinForms(originalLinearForm, linFormStr.c_str());
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

	if ( cmd == "range" )
	{
		BoxOptimization bo;
		bo.setPolynomial(lowerBound, upperBound, originalPolynomial);
		bo.setPower(5);
		bo.findRange(10);
	}

	if ( cmd == "opt")
	{
		BoxOptimization bo;
		bo.setPolynomial(lowerBound, upperBound, originalPolynomial);


		RR N;
		cout << "epsilon=" << epsilon << endl;
		N = 1;
		for (int i = 0;  i < lowerBound.length(); ++i)
			N *= to_RR(upperBound[i] - lowerBound[i] + 1);
		N = ceil((1.0 + inv(epsilon))* log(N));
		int k = INT_MAX;
		if ( N < INT_MAX)
			k = to_int(N);
		else
			cout << "Warning: k is larger than " << k << endl;


		cout << "starting k was " << k << endl;
		if (userK > 0)
			k = userK;
		k--;
		while ( bo.maximumUpperbound() - bo.maximumLowerBound() > 0.01)
		{
			k++;
			bo.setPower(k);
			bo.findRange(10);
			bo.findNewUpperbound();
			cout << "k: " << k << " " << bo.L << " <= f(x) <= " << bo.U << "\n";
			cout << "k: " << k << " " << bo.maximumLowerBound() << " <= max f(x) <= " << bo.maximumUpperbound() << "\n";
			cout << "gap: " << bo.maximumUpperbound() - bo.maximumLowerBound() << "\n";
			//if (k > 10)
				break;
		}
		cout << "k: " << k << " " << bo.L << " <= f(x) <= " << bo.U << "\n";
		cout << "k: " << k << " " << bo.maximumLowerBound() << " <= max f(x) <= " << bo.maximumUpperbound() << "\n";
		cout << "gap: " << bo.maximumUpperbound() - bo.maximumLowerBound() << "\n";
	}


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


	return 0;
}


int main2(int argc, const char *argv[]) {
	ReadPolyhedronData read_polyhedron_data;
	string linFormFileName;

	struct BarvinokParameters *params = new BarvinokParameters;

	latte_banner(cerr);

	cerr << "Invocation: ";
	for (int i = 0; i < argc; i++) {
		cerr << argv[i] << " ";
	}
	cerr << endl;

	params->substitution = BarvinokParameters::PolynomialSubstitution;
	params->decomposition = BarvinokParameters::DualDecomposition;
	params->max_determinant = 1;
	for (int i = 1; i < argc; i++) {
		if (read_polyhedron_data.parse_option(argv[i])) {
		} 
		else if (strncmp(argv[i], "--linFile=", 10) == 0){
			linFormFileName = string(argv[i] + 10);
		} else {
			cerr << "Unknown command/option " << argv[i] << endl;
			THROW_LATTE_MSG(LattException::ue_BadCommandLineOption, argv[i]);
		}
	} //for i.

	if (read_polyhedron_data.expect_filename) {
		cerr << "Filename missing" << endl;
		THROW_LATTE(LattException::ue_FileNameMissing);
	}

	const char *fileName = read_polyhedron_data.filename.c_str();

	Polyhedron *Poly = read_polyhedron_data.read_polyhedron(params);

	cout << "fix me**********************************************" << endl;
	params->Number_of_Variables = Poly->numOfVars;
	
	if (Poly->cones != NULL && Poly->cones->rays == NULL) {

			// Only facets computed, for instance by using the 4ti2
			// method of computing vertex cones.  So dualize twice to
			// compute the rays.
			cerr << "(First computing their rays... ";
			cerr.flush();
			dualizeCones(Poly->cones, Poly->numOfVars, params);
			dualizeCones(Poly->cones, Poly->numOfVars, params); // just swaps
			cerr << "done; rays are now computed) \n";
			cerr.flush();

	}
	//printListCone(Poly->cones, Poly->numOfVars);
	//cout << "*****************" << endl;
			
	//cout << "numrays=" << lengthListVector(Poly->cones->rays) << endl;
	//cout << "num vars" << params->Number_of_Variables << endl;
	//cout << Poly->numOfVars;
	
	assert(Poly->cones->rays != NULL);

	//Weighted_Exponential_Single_Cone_Parameters *exp_param =
	//		new Weighted_Exponential_Single_Cone_Parameters(*params);
	//delete params;
	//params = exp_param;

	if ( linFormFileName.length() == 0) {
		cerr << "lin fomr file missing: " << endl;
		THROW_LATTE(LattException::ue_FileNameMissing);
	}
	
	ifstream linFormFile(linFormFileName.c_str());
	string linFormStr;
	getline(linFormFile, linFormStr);
	cout << "lin form str: " << linFormStr.c_str() << endl;
	linFormSum originalLinearForm;
	loadLinForms(originalLinearForm, linFormStr.c_str());
		

	cout << "what should happen next in the general setting?" << endl;
	mpq_class weighted_count = computeWeightedExponentialResidue(Poly->cones, params, originalLinearForm);
	cout << "Final count: " << weighted_count << endl;
	destroyLinForms(originalLinearForm);

	return 0;
}
