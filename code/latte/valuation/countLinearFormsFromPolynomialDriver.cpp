/*
 * countLinearFormsFromPolynomialDrover.cpp
 *
 *  Created on: March 27, 2010
 *      Author: Brandon Dutra and Gregory Pinto
 *
 *  Given a list of polynomial files, counts the number of linear forms in the decomposition.
 */

#include <cstdlib>
#include <iostream>
#include <string>
#include <sstream>

#include <NTL/vec_ZZ.h>
#include <NTL/vec_RR.h>

/* Integration Headers */
#include "integration/PolyTrie.h"
#include "integration/newIntegration.h"


using namespace std;

int main(int argc, char * argv[])

{
	ZZ sumCount, minCount, maxCount;
	sumCount = 0;
	minCount = 0;
	maxCount = 0;
	
	
	if ( argc == 1)
	{
		cout << "usage: " << argv[0] << " polynomial files" << endl;
		exit(1);
	}
	
	for(int i = 1; i < argc; ++i)	
	{
		//read the polynomial from the file.
		ifstream inFile;
		string polynomialLine;
		polynomialLine = "";
		
		inFile.open(argv[i]);
		getline(inFile, polynomialLine, '\n');
		inFile.close();
	
		//convert to polynomial datastructure.
		monomialSum originalPolynomial;
		loadMonomials(originalPolynomial, polynomialLine);
		
		//convert to linear form
		linFormSum linearForms;
		BTrieIterator<RationalNTL, int>* polynomialIterator =
			new BTrieIterator<RationalNTL, int> ();
		
		linearForms.termCount = 0;
		linearForms.varCount = originalPolynomial.varCount;
	
			
	
		polynomialIterator->setTrie(originalPolynomial.myMonomials,
			originalPolynomial.varCount);
		decompose(polynomialIterator, linearForms);

		destroyMonomials(originalPolynomial);


		//finally, count the terms.
		sumCount += linearForms.termCount;
		if ( minCount > linearForms.termCount || minCount == 0)
			minCount = linearForms.termCount;
		if ( maxCount < linearForms.termCount || maxCount == 0)
			maxCount = linearForms.termCount;
		
		//delete the linear form.
		destroyLinForms(linearForms);
	}//for each file.

	cout << "Total count " << sumCount
		<< "\nAvg count " << to_RR(sumCount)/argc
		<< "\nMin Count " << minCount
		<< "\nMax Count " << maxCount;

	return 0;
}

