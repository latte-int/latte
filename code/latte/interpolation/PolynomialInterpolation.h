/*
 *
 *
 *  Created on: May 13, 2010
 *      Author: bedutra
 */

#ifndef POLYNOMIAL_INTERPOLATION_H_
#define POLYNOMIAL_INTERPOLATION_H_

#include <iomanip>
#include <vector>
#include "gmp.h"
#include <gmpxx.h>

using namespace std;

class PolynomialInterpolation
{
private:
    int colSize; 		//number of columns
    int rowSize;		//number of rows.
    mpq_class **matrix; // a rowSize X (colSize +1) matrix. matrix = [A | b]
    int initCounter;	//counts the number of rows added during initialization.
    
public:
    PolynomialInterpolation(unsigned int degree);
	
	//A B C D E F G H I J K L M N O P Q R S T U V W X Y Z

	void addMultRows(mpq_class &value, int fromRow,  int toRow);
	void addPoint(mpq_class x, mpq_class f);
	mpq_class evaluatePoly(mpq_class const & xValue, vector<mpq_class>  const & poly);
	void GE();
	vector<mpq_class> getSolutionVector();
	void printMatrix() const;
	bool isSingular();
	vector<mpq_class> solve();
	void swap(int a, int b);
	void timesRow(int row, mpq_class value);
	
	
	
	PolynomialInterpolation & operator=(const PolynomialInterpolation & rhs);
};//PolynomialInt



#endif /* POLYNOMIAL_INTERPOLATION_H_ */

