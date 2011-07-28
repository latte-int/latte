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

int main(int argc, char *argv[])
{
	if ( argc != 5)
	{
		cout << "Usage: " << argv[0] << " correct-answer [p or l] [polynomial-file or linear-form-file] latte-file\n";
		cout << "the correct-answer could be a fraction, but should not have spaces.\n\n";
		cout << "Description: Tests that integrating a polynomial or power of a linear form over the polytope gives the same answer as what is passed in." << endl;
		return 1;
	}


	RationalNTL correctAnswer(argv[1]);
	Valuation::ValuationContainer valuationContainerLaw, valuationContainerTri;

	char * options[5];

	options[0] = "./exe";
	options[1] = "--valuation=integrate";
	options[2] = "--all";
	if ( strcmp(argv[2], "p") == 0)
	{
		options[3] = new char[strlen(argv[3]) + 13];
		strcpy(options[3], "--monomials=");
		strcat(options[3], argv[3]);
	}
	else if ( strcmp(argv[2], "l") == 0)
	{
		options[3] = new char[strlen(argv[3]) + 16];
		strcpy(options[3], "--linear-forms=");
		strcat(options[3], argv[3]);
	}
	else
	{
		cout << "Expected 'p' or 'l', but recieved " << argv[2] << endl;
		exit(1);
	}
	options[4] = argv[4];


	char *lawrenceOptions[5], *triOptions[5];
	for (int i = 0; i < 5; ++ i)
	{
		lawrenceOptions[i] = options[i];
		triOptions[i] = options[i];
	}


	//if we want to test both methods starting from the vertex-cone or the lifted-cone representation.
#if 1
	lawrenceOptions[2] = "--cone-decompose";
	triOptions[2] = "--triangulate";

	//run lawrence and triangulate on their own (not using --all) to make sure the tangent and lifted cone methods are working.
	cout << "Running main valuation Driver with Lawrence" << endl;
	valuationContainerLaw = Valuation::mainValuationDriver((const char **)lawrenceOptions, 5);

	cout << "Running main valuation Driver with Triangulate" << endl;
	valuationContainerTri = Valuation::mainValuationDriver((const char **)triOptions, 5);
#else //else, test both methods starting from vertex-cones.
	cout << "Running main valuation Driver with both" << endl;
	valuationContainerLaw = Valuation::mainValuationDriver((const char **)lawrenceOptions, 5);
	valuationContainerTri = valuationContainerLaw;
#endif

	delete options[3];

	RationalNTL triAns, lawAns;
	for(int i = 0; i < valuationContainerLaw.answers.size(); ++i)
	{
		if (valuationContainerLaw.answers[i].valuationType == Valuation::ValuationData::integrateLawrence)
			lawAns = valuationContainerLaw.answers[i].answer;
	}


	for(int i = 0; i < valuationContainerTri.answers.size(); ++i)
	{
		if (valuationContainerTri.answers[i].valuationType == Valuation::ValuationData::integrateTriangulation)
			triAns = valuationContainerTri.answers[i].answer;
	}

	if ( triAns != correctAnswer || lawAns != correctAnswer)
	{
		cout << "******ERROR*******" << endl;
		cout << "correct answer    = " << correctAnswer << endl;
		cout << "Lawrence answer   = " << lawAns << endl;
		cout << "Triangulate anser = " << triAns << endl;
		cout << "Input options:" << endl;
		for(int i = 0; i < 5; ++i)
			cout << "     " << options[i] << endl;

		return 2;
	}

	return 0; //no errors!
}
//main()
