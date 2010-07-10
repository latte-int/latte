/*
 * PolynomialInterpolation.cpp
 *
 *  Created on: May 13, 2010
 *      Author: bedutra
 */
#include "PolynomialInterpolation.h"
#include <iostream>
#include <cstdlib>

using namespace std;


/** Find the coefficients of any degree d polynomial over Q. 
 * degree: the degree of the unknown poly.
 */
PolynomialInterpolation::PolynomialInterpolation(unsigned int degree):
	colSize(degree + 1),
	rowSize(degree + 1),
	matrix(NULL),
	initCounter(0)
{
	
	matrix = new mpq_class*[rowSize];
	for(int row = 0; row < rowSize; ++row)
	{
		matrix[row] = new mpq_class[colSize + 1]; // [A | b}
		for(int col = 0; col <= colSize; ++col)
		{
			//mpq_init(matrix[row][col]);
			matrix[row][col] = 0;
		}
	}
}//constructor


//A B C D E F G H I J K L M N O P Q R S T U V W X Y Z	

//sets matrix(toRow,:) = matirx(toRow,:) - value * matrix(fromRow, :) (matlab syntax)
void PolynomialInterpolation::addMultRows(mpq_class  &value, int fromRow, int toRow)
{
	value.canonicalize();
	for(int i  = 0; i <= colSize; ++i)
	{
		matrix[toRow][i] = matrix[toRow][i] - value * matrix[fromRow][i];
	}//for i
}//addMultRows



//adds the point (x, f(x)) to the matrix.
//  adds the row [1, x, x^2, ..., f(x)] to the matrix.
void PolynomialInterpolation::addPoint(mpq_class x, mpq_class f)
{
	mpq_class powerOfx(1);
	f.canonicalize();
	x.canonicalize();

	if ( initCounter >= rowSize)
	{
		cerr << "Ops, you have too many points" << endl;
		exit(1);
	}

	for(int col = 0; col < colSize; ++col)
	{
		matrix[initCounter][col] = powerOfx;
		powerOfx = powerOfx * x;
	}

	matrix[initCounter][colSize] = f;
	initCounter++;
}//addPoint


/*	Returns f(xValue) where f is a polynomial. poly[i] is the coefficient of x^i.
 *
 */
mpq_class PolynomialInterpolation::evaluatePoly(mpq_class const & xValue, vector<mpq_class>  const & poly)
{
	mpq_class power, sum;
	sum = 0;
	power = 1;

	sum.canonicalize();
	power.canonicalize();

	for( int i = poly.size() - 1; i > 0; --i)
	{
		sum = (sum + poly[i])* xValue;
		//sum.canonicalize();
		//if ( ! isReduced(sum) )
		//{
		//	cout << "ops, not in reduced form: " << sum << endl;
		//	exit(1);
		//}
	}//for i
	mpq_class finalSum = sum + poly[0];
	finalSum.canonicalize();
	return finalSum;
}//evaluates.


void PolynomialInterpolation::GE()
{
	//PolynomialInterpolation copy(rowSize - 1);
	//copy = *this;
	
	int perfectRow = 0; //the row with only zeros on the left.
	int currentColumn = 0; 
    while( currentColumn < colSize && perfectRow < rowSize)
    {

        //find row with non-zero base
        int currentRow = perfectRow;
        while( currentRow < rowSize && matrix[currentRow, currentColumn] == 0)
        {
            currentRow = currentRow + 1;
        }
        
        if(  currentRow > rowSize)
		{
			currentColumn++;
            continue; //found column of zeros
        }
        else
        {
        	swap(currentRow, perfectRow);
        }//swap rows so that matrix[perfectRow][currentColum] is non-zero

        if(matrix[perfectRow][currentColumn] == 0 )
        {
            cerr << "GE:assert matirx[perfectRow][currentColumn] != 0" << endl;
            //cout << "perfectRow = " << perfectRow << ", curCol=" << currentColumn << endl;
            //cout << "row, col size=" << rowSize << ", " << colSize << endl;
            //printMatrix();
            //cout << "origional matirx\n";
            //copy.printMatrix();
            exit(1);
        }
        
        timesRow(perfectRow, mpq_class(matrix[perfectRow][currentColumn].get_den(), matrix[perfectRow][currentColumn].get_num()));

        if (  matrix[perfectRow][currentColumn] != 1)
        {
        	cerr << "GE::assert matrix[perfectRow][currentColumn] == 1 failed" << endl;
        	//cerr << "pr = " << perfectRow << ", curCol=" << currentColumn << endl;
        	//printMatrix();
        	exit(1);
        }	

        //use row 'perfectRow' to eliminate every other row.
        for( int i = 0; i < rowSize; ++i)
        {
            if( matrix[i][currentColumn] != 0 && i != perfectRow)
            {  
            	mpq_class value(matrix[i][currentColumn].get_num(), matrix[i][currentColumn].get_den());

            	addMultRows(value, perfectRow, i);
                if ( matrix[i][currentColumn] != 0)
                {
                	cerr << "GE::matrix[i][currentColumn] != 0 failed." << endl;
                	exit(1);
                }
            }//do row additions on every row other than the current row.
        }//for i.
		perfectRow++;
		currentColumn++;
    }//end while
}//GE()


/** return b if the matrix is in the form [Id | b]
  * else calls exit if the matrix is singular.
  */
vector<mpq_class> PolynomialInterpolation::getSolutionVector()
{
	vector<mpq_class> answer(rowSize);
	if ( isSingular() )
	{
		cerr << "ops, there is a nullspace!" << endl;
		printMatrix();
		exit(1);
	}
	
	for(int i = 0; i < rowSize; ++i)
	{
		answer[i] = matrix[i][colSize];
	}
	return answer;
}//getSolutionVector


void PolynomialInterpolation::printMatrix() const
{
	cout << "PRINT MATRIX\n";
	for(int row = 0; row < rowSize; ++row)
	{
		for(int col = 0; col <= colSize; ++col)
		{
			cout << setw(10) << matrix[row][col] << ", ";
		}//for col
		cout << endl;
	}//for row.
}//printMatrix

/** Returns true if A != Id, where matrix = [A | b]
 *
 */
bool PolynomialInterpolation::isSingular()
{
	for(int i = 0; i < rowSize; ++i)
		for(int k = 0; k < colSize; ++k)
			if ( matrix[i][k] != 0 && i != k)
				return true;
			else if ( matrix[i][k] != 1 && i == k) //matrix[i][i] could be 0
				return true;
	return false;
}//isSingular


/** Does GE and returns the coeff. vector.
 *
 */
vector<mpq_class> PolynomialInterpolation::solve()
{
	GE();
	return getSolutionVector();
}

//swap two rows.
void PolynomialInterpolation::swap(int a, int b)
{
	mpq_class *temp;
	temp = matrix[a];
	matrix[a] = matrix[b];
	matrix[b] = temp;
}

/* Times row "row" (including the b-column) by "value"
 *
 */
void PolynomialInterpolation::timesRow(int row, mpq_class value)
{
	value.canonicalize();
	for(int i = 0; i <= colSize; ++i)
	{
		matrix[row][i] = matrix[row][i] * value;
	}//for i. Also times the constant term.

}//timesRow


PolynomialInterpolation & PolynomialInterpolation::operator=(const PolynomialInterpolation & rhs)
{

	if (this == &rhs)      // Same object?
      return *this;        // Yes, so skip assignment, and just return *this.
	
	colSize = rhs.colSize;
	rowSize = rhs.rowSize;
	initCounter = rhs.initCounter;
	
	matrix = new mpq_class*[rowSize];
	for(int row = 0; row < rowSize; ++row)
	{
		matrix[row] = new mpq_class[colSize + 1]; // [A | b}
		for(int col = 0; col <= colSize; ++col)
		{
			matrix[row][col] = rhs.matrix[row][col];
		}
	}//for row
}//operator=
	
