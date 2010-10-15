/*
 * testIntegrationHyperrectangles.cpp
 *
 *  Created on: Jul 21, 2010
 *      Author: bedutra
 *
 *  usage: exe correct-answer polynomial-file latte-file
 */

#include <iostream>
#include <cstring>
#include "valuation.h"
#include "rational.h"

using namespace std;

int main(int argc, char *argv[])
{
	if ( argc != 4)
	{
		cout << "Usage: " << argv[0] << " correct-answer polynomial-file latte-file\n";
		cout << "the correct-answer could be a fraction, but should not have spaces.\n\n";
		cout << "Description: Tests that integrating the polynomial over the polytope gives the same answer as what is passed in." << endl;
		return 1;
	}


	RationalNTL correctAnswer(argv[1]);
	Valuation::ValuationContainer valuationContainer;

	char * options[5];

	options[0] = "./exe";
	options[1] = "--valuation=integrate";
	options[2] = "--triangulate";
	options[3] = new char[strlen(argv[2]) + 13];
	strcpy(options[3], "--monomials=");
	strcat(options[3], argv[2]);
	options[4] = argv[3];

	cout << "Running main valuation Driver" << endl;
	valuationContainer = Valuation::mainValuationDriver((const char **)options, 5);
	delete options[3];


	if ( valuationContainer.answers[0].answer != correctAnswer)
	{
		cout << "******ERROR*******" << endl;
		cout << "correct answer  = " << correctAnswer << endl;
		cout << "computed answer = " << valuationContainer.answers[0].answer << endl;

		return 2;
	}

	return 0; //no errors!
}
//main()
