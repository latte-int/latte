/*
 * testIntegrationHyperrectanglesDriver.cpp
 *
 *  Created on: Jul 21, 2010
 *      Author: bedutra
 *
 *  usage: exe correct-answer [p or l] [polynomial-file or linear-form-file] latte-file
 *  checks the result from the valuation class with the passed-in correct value.
 *  This exe is used by a maple script.
 */

#include <iostream>
#include <cstring>
#include "../valuation.h"
#include "rational.h"

using namespace std;


int checkPolynomial(const char *correctAns, const char * integrandFile, const char * latteFile)
{

	RationalNTL correctAnswer(correctAns);
	Valuation::ValuationContainer valuationContainerLaw, valuationContainerTri;

	char * options[4];

	options[0] = (char*)"./exe";
	options[1] = (char*)"--valuation-alg=poly-lf-cone";
	options[2] = new char[strlen(integrandFile) + 13];
	strcpy(options[2], "--monomials=");
	strcat(options[2], integrandFile);
	options[3] = (char*)latteFile;


	cout << "Running main valuation Driver with Lawrence" << endl;
	valuationContainerLaw = Valuation::mainValuationDriver((const char **)options, 4);

	options[1] = (char*)"--valuation-alg=poly-lf-triangulation";
	cout << "Running main valuation Driver with Triangulate" << endl;
	valuationContainerTri = Valuation::mainValuationDriver((const char **)options, 4);

	delete options[2];

	RationalNTL triAns, lawAns;
	for(int i = 0; i < valuationContainerLaw.answers.size(); ++i)
	{
		if (valuationContainerLaw.answers[i].valuationType == PolytopeValuation::integratePolynomialAsLinearFormCone)
			lawAns = valuationContainerLaw.answers[i].answer;
	}


	for(int i = 0; i < valuationContainerTri.answers.size(); ++i)
	{
		if (valuationContainerTri.answers[i].valuationType == PolytopeValuation::integratePolynomialAsLinearFormTriangulation)
			triAns = valuationContainerTri.answers[i].answer;
	}



	if ( triAns != correctAnswer || lawAns != correctAnswer)
	{
		cout << "******ERROR*******" << endl;
		cout << "correct answer    = " << correctAnswer << endl;
		cout << "Lawrence answer   = " << lawAns << endl;
		cout << "Triangulate anser = " << triAns << endl;
		cout << "Input options:" << endl;
		for(int i = 0; i < 4; ++i)
			cout << "     " << options[i] << endl;
			return 2;
	}


	return 0;
}

int checkLinearForm(const char *correctAns, const char * integrandFile, const char * latteFile)
{

	RationalNTL correctAnswer(correctAns);
	Valuation::ValuationContainer valuationContainerLaw, valuationContainerTri;

	char * options[4];

	options[0] = (char*)"./exe";
	options[1] = (char*)"--valuation-alg=lf-cone";
	options[2] = new char[strlen(integrandFile) + 16];
	strcpy(options[2], "--linear-forms=");
	strcat(options[2], integrandFile);
	options[3] = (char*)latteFile;


	cout << "Running main valuation Driver with Lawrence" << endl;
	valuationContainerLaw = Valuation::mainValuationDriver((const char **)options, 4);

	options[1] = (char*)"--valuation-alg=lf-triangulation";
	cout << "Running main valuation Driver with Triangulate" << endl;
	valuationContainerTri = Valuation::mainValuationDriver((const char **)options, 4);

	delete options[2];

	RationalNTL triAns, lawAns;
	for(int i = 0; i < valuationContainerLaw.answers.size(); ++i)
	{
		if (valuationContainerLaw.answers[i].valuationType == PolytopeValuation::integrateLinearFormCone)
			lawAns = valuationContainerLaw.answers[i].answer;
	}


	for(int i = 0; i < valuationContainerTri.answers.size(); ++i)
	{
		if (valuationContainerTri.answers[i].valuationType == PolytopeValuation::integrateLinearFormTriangulation)
			triAns = valuationContainerTri.answers[i].answer;
	}



	if ( triAns != correctAnswer || lawAns != correctAnswer)
	{
		cout << "******ERROR*******" << endl;
		cout << "correct answer    = " << correctAnswer << endl;
		cout << "Lawrence answer   = " << lawAns << endl;
		cout << "Triangulate anser = " << triAns << endl;
		cout << "Input options:" << endl;
		for(int i = 0; i < 4; ++i)
			cout << "     " << options[i] << endl;
			return 2;
	}
	return 0;
}

int checkProductLinearForm(const char *correctAns, const char * integrandFile, const char * latteFile)
{

	RationalNTL correctAnswer(correctAns);
	Valuation::ValuationContainer valuationContainerTri;

	char * options[4];

	options[0] = (char*)"./exe";
	options[1] = (char*)"--valuation-alg=plf-triangulation";
	options[2] = new char[strlen(integrandFile) + 25];
	strcpy(options[2], "--product-linear-forms=");
	strcat(options[2], integrandFile);
	options[3] = (char*)latteFile;

	cout << "Running main valuation Driver with Triangulate" << endl;
	valuationContainerTri = Valuation::mainValuationDriver((const char **)options, 4);

	delete options[2];

	RationalNTL triAns;

	for(int i = 0; i < valuationContainerTri.answers.size(); ++i)
	{
		if (valuationContainerTri.answers[i].valuationType == PolytopeValuation::integrateProductLinearFormsTriangulation)
			triAns = valuationContainerTri.answers[i].answer;
	}



	if ( triAns != correctAnswer)
	{
		cout << "******ERROR*******" << endl;
		cout << "correct answer    = " << correctAnswer << endl;
		cout << "Triangulate anser = " << triAns << endl;
		cout << "Input options:" << endl;
		for(int i = 0; i < 4; ++i)
			cout << "     " << options[i] << endl;
			return 2;
	}
	return 0;
}

int main(int argc, char *argv[])
{
	if ( argc != 5)
	{
		cout << "Usage: " << argv[0] << " correct-answer [p or l or d] [polynomial-file or linear-form-file or product-linear-form-file] latte-file\n";
		cout << "the correct-answer could be a fraction, but should not have spaces.\n\n";
		cout << "Description: Tests that integrating a polynomial or power of a linear form or a product of linear forms over the polytope gives the same answer as what is passed in." << endl;
		return 1;
	}


	if ( argv[2][0] == 'p')
		return checkPolynomial(argv[1], argv[3], argv[4]); //checkPolynomial(const char *correctAns, const char * integrandFile, const char * latteFile);
	else if (argv[2][0] == 'l')
		return checkLinearForm(argv[1], argv[3], argv[4]);
	else if ( argv[2][0] == 'd')
		return checkProductLinearForm(argv[1], argv[3], argv[4]);
	else
		THROW_LATTE(LattException::ue_BadCommandLineOption);

}
//main()
