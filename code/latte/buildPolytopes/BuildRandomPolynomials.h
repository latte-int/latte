/*
 * BuildRandomPolynomials.h
 *
 *  Created on: Nov 30, 2010
 *      Author: bedutra
 */

#include <string>
using namespace std;

#ifndef BUILDRANDOMPOLYNOMIALS_H_
#define BUILDRANDOMPOLYNOMIALS_H_

//makes a random polynomial in dim variables of a set degree.
string makeRandomMonomial(const int dim, int totalDegree);

//makes many random monomials
string makeRandomPolynomial(const int dim, const int total, const int numMonomials);


#endif /* BUILDRANDOMPOLYNOMIALS_H_ */
