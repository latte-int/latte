/* hello
 * polynomial.h
 *
 *  Created on: May 13, 2010
 *      Author: bedutra
 */

#ifndef POLYNOMIAL_H_
#define POLYNOMIAL_H_

#include <iomanip>
#include <vector>
#include "gmp.h"
#include <gmpxx.h>

using namespace std;

class PolynomialInt
{
private:
    int colSize;
    int rowSize;
    mpq_class **matrix;
    int initCounter; //counts the number of rows added during initialization.
    
	
	//mpf_class floating point?
	//mpz_t this is GMP integer
public:
	PolynomialInt(unsigned int degree);
	
	void addPoint(mpq_class x, mpq_class f);
	void addMultRows(mpq_class value, int fromRow, int toRow);
	void timesRow(int row, mpq_class value);
	void swap(int a, int b);
	void GE();
	void solve() {};
	bool isSingular();
	
	vector<mpq_class> getSolutionVector();
	
	void printMatrix();
	
	PolynomialInt & operator=(const PolynomialInt & rhs);
	

};//PolynomialInt



#endif /* POLYNOMIAL_H_ */

