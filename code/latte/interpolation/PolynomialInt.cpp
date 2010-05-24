/*
 * polynomial.cpp
 *
 *  Created on: May 13, 2010
 *      Author: bedutra
 */
#include "PolynomialInt.h"
#include <iostream>

using namespace std;


/** Find the coefficients of any degree d polynomial over Q. 
 * degree: the degree of the unknown poly.
 */
PolynomialInt::PolynomialInt(unsigned int degree):
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

//adds the point (x, f(x)) to the matrix.
void PolynomialInt::addPoint(mpq_class x, mpq_class f)
{
	mpq_class powerOfx(1);
	//mpq_init(powerOfx);
	//mpq_set_ui(powerOfx, 1, 1);
	f.canonicalize();
	x.canonicalize();
	//mpq_canonicalize(x);
	//mpq_canonicalize(f);
	//mpq_canonicalize(powerOfx);
	
	if ( initCounter >= rowSize)
	{
		cerr << "Ops, you have too many points" << endl;
		//exit(1);
	}

	for(int col = 0; col < colSize; ++col)
	{
		//mpq_set(matrix[initCounter][col], powerOfx);
		//mpq_mul(powerOfx, powerOfx, x);
		matrix[initCounter][col] = powerOfx;
		powerOfx = powerOfx * x;
	}
	//mpq_set(matrix[initCounter][colSize], f);
	matrix[initCounter][colSize] = f;
	initCounter++;
}//addPoint



void PolynomialInt::swap(int a, int b)
{
	mpq_class *temp;
	temp = matrix[a];
	matrix[a] = matrix[b];
	matrix[b] = temp;
}

void PolynomialInt::timesRow(int row, mpq_class value)
{
	value.canonicalize();
	for(int i = 0; i <= colSize; ++i)
	{
		matrix[row][i] = matrix[row][i] * value; 
	}//for i. Also times the constant term.

}

void PolynomialInt::GE()
{

	int perfectRow = 0; //the row with only zeros on the left.
	int currentColumn = 0; 
    while( currentColumn < colSize && perfectRow < rowSize)
    {
    
    	//cout << "loop interation: curCol=" << currentColumn << ", perRow=" << perfectRow << endl;
    	//printMatrix();

        //find row with non-zero base
        int currentRow = perfectRow;
        while( currentRow < rowSize && matrix[currentRow, currentColumn] == 0)
        {
            currentRow = currentRow + 1;
        }
        
        if(  currentRow > rowSize)
		{
			currentColumn++;
            continue; //colum of zeros
        }
        else
        {
        	swap(currentRow, perfectRow);
            //temp = M(i,:);
            //M(i,:) = M(perfectRow,:);
            //M(perfectRow,:) = temp;
        }//swap rows so that matrix[k][k'] is non-zero

        if(matrix[perfectRow][currentColumn] == 0 )
        {
            cerr << "GE:assert matirx[perfectRow][currentColumn= != 0" << endl;
            exit(1);
        }
        
                
        timesRow(perfectRow, mpq_class(matrix[perfectRow][currentColumn].get_den(), matrix[perfectRow][currentColumn].get_num()));
        //cout << "after clearing perfRow=" << perfectRow << endl;
        //printMatrix();
        
        
        if (  matrix[perfectRow][currentColumn] != 1)
        {
        	cerr << "GE::assert matrix[perfectRow][currentColumn] == 1 failed" << endl;
        	cerr << "pr = " << perfectRow << ", curCol=" << currentColumn << endl;
        	printMatrix();
        	exit(1);
        }	

        for( int i = 0; i < rowSize; ++i)
        {
            if( matrix[i][currentColumn] != 0 && i != perfectRow)
            {  
            	//cout << "no seg fault yet" << endl;
            	mpq_class value(matrix[i][currentColumn].get_num(), matrix[i][currentColumn].get_den());
            	//value.canonicalize();
            	//cout << "making values seg fault?:: value= " << value << ", perRow=" << perfectRow << ", i= " << i << endl;
            	addMultRows(value, perfectRow, i);
            	//cout << "maybe addMultRow seg faults?" << endl;
                //M(i,:) = xor(M(i,:), M(perfectRow,:));  
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
    //B = M;
    //nullspace = getNullBasis(B);

}//GE()


/** return b if the matrix is in the form [Id | b]
  * else calles exit
  */
vector<mpq_class> PolynomialInt::getSolutionVector()
{
	vector<mpq_class> answer(rowSize);
	if (! isSingular() )
	{
		cerr << "ops, there is a nullspace!" << endl;
		exit(1);
	}
	
	for(int i = 0; i < rowSize; ++i)
		answer[i] = matrix[i][colSize+1];
		
	return answer;
}//getSolutionVector


//sets matrix(toRow,:) = matirx(toRow,:) - value * matrix(fromRow, :) (matlab syntax)
void PolynomialInt::addMultRows(mpq_class value, int fromRow, int toRow)
{
	value.canonicalize();
	for(int i  = 0; i <= colSize; ++i)
	{
		matrix[toRow][i] = matrix[toRow][i] - value * matrix[fromRow][i];
	}//for i
}//addMultRows



void PolynomialInt::printMatrix()
{
	for(int row = 0; row < rowSize; ++row)
	{
		for(int col = 0; col <= colSize; ++col)
		{
			//mpq_out_str(NULL, 10, matrix[row][col]);
			cout << setw(10) << matrix[row][col];
			cout << ", ";
		}
		cout << endl;
	}//for row.


}//printMatrix


/** Returns true if A != Id, where matrix = [A | b]
 *
 */
bool PolynomialInt::isSingular()
{
	for(int i = 0; i < rowSize; ++i)
		for(int k = 0; k < colSize; ++k)
			if ( matrix[i][k] != 0 && i != k)
				return true;
			else if ( matrix[i][k] != 1 && i == k) //matrix[i][i] could be 0
				return true;
	return false;
}
