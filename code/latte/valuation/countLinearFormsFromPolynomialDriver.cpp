/*
 * countLinearFormsFromPolynomialDrover.cpp
 *
 *  Created on: March 27, 2010
 *      Author: Brandon Dutra and Gregory Pinto
 *
 *  Given a list of monomial files, creates a file containing the number of linear forms from that file.
 */

#include <cstdlib>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <cassert>

#include <NTL/vec_ZZ.h>
#include <NTL/vec_RR.h>

/* Integration Headers */
#include "integration/PolyTrie.h"
#include "integration/newIntegration.h"


using namespace std;

int main(int argc, char * argv[])

{
	if ( argc == 1)
	{
		cout << "usage: " << argv[0] << " polynomial files" << endl;
		exit(1);
	}
	
	ofstream file;
	file.open((string(argv[0]) + ".log").c_str(), ios::app);

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

		//find the degree of the monomial.
		int degree;
		term<RationalNTL, int>* temp;
		polynomialIterator->begin();
		temp = polynomialIterator->nextTerm();
		degree = 0;
		for(int k = 0; k < originalPolynomial.varCount; ++k)
			degree += temp->exps[k];
		temp = polynomialIterator->nextTerm();
		assert(temp == NULL);


		decompose(polynomialIterator, linearForms);

		destroyMonomials(originalPolynomial);

		file << linearForms.varCount << " " << degree << " " << linearForms.termCount << " " << argv[i] << "\n";
		//delete the linear form.
		destroyLinForms(linearForms);
	}//for each file.
	file.close();

	cout << "Results printed to file " << argv[0] << ".log" << endl;
	return 0;
}

