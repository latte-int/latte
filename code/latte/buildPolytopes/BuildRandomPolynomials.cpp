/*
 * BuildRandomPolynomials.cpp
 *
 *  Created on: Nov 30, 2010
 *      Author: bedutra
 */

#include <vector>
#include <sstream>
#include "BuildRandomPolynomials.h"
#include <cstdlib>

/**
 * makes a random monomial in dim variables of a set degree.
 * The monomial is monic.
 * If you want only 1 monomial in your polynomial, you must call makeRandomPolynomial or enclose the string in the final '[]' brackets.
 */
string makeRandomMonomial(const int dim, int totalDegree)
{
	//We treat totalDegree as the amount left over...how much powers we still have to add.
	const int stepSize = 2; //add no more than 1 powers to a term at a time.
	vector<int> powers;
	int i; //i = the x_i who will be getting more powers.
	int newPower;
	stringstream poly;

	powers.resize(dim);


	while (totalDegree > 0 )
	{
		i = rand()%dim; //find the x_i who will get more powers.
		newPower = rand()%stepSize; //find the additional power.

		powers[i] += newPower;
		totalDegree -= newPower;
	}//add more powers to the polynomial

	if (totalDegree < 0)
	{
		powers[i] += totalDegree; //totalDegree is neg!
	}//if we added too much, subtract from the last term we added to.

	//now make a string monomial.
	poly << "[1,[";
	for(size_t j = 0; j < powers.size(); ++j)
	{
		poly << powers[j];
		if (j != powers.size()-1)
			poly << ',';
	}
	poly << "]]";
	return poly.str();
}


/* makes many random monic monomials of the same degree.
 * If you want 1 monomial, you must use this function and just set numMonomials to 1.
 *
 */
string makeRandomPolynomial(const int dim, const int totalDegree, const int numMonomials)
{
	stringstream poly;

	poly << "[";
	for(int i = 0; i < numMonomials; ++i)
	{
		poly << makeRandomMonomial(dim, totalDegree);
		if ( i < numMonomials-1)
			poly << ',';
	}
	poly << "]";
	return poly.str();
}

