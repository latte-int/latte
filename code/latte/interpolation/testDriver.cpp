/*
 * testDriver.cpp
 *
 *  Created on: May 13, 2010
 *      Author: bedutra
 */

#include <iostream>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include "PolynomialInterpolation.h"
#include <gmp.h>
#include <gmpxx.h>

#define MAXDEG 40 		//max degree of a poly.
#define PNEG .5 		//prop of a neg coeff.
#define MAXCOEF 200000 	//max coeff



using namespace std;


mpq_class evaluatePoly(mpq_class const & xValue, vector<mpq_class>  const & poly);
bool isReduced(mpq_class r);
void test1Poly();

int main()
{
	srand(44);
/*
	//c style 
    int polyTopC[] = {4, 5, -20, 4, -8, 9};
    int polyBotC[] = {1, 2,   1, 3,  5, 10};
    //int polyBotC[] = {1, 1,   1, 1,  1, 1};
    int const POLYSIZE = 6;

    int xValueTopC[] = {1, 2, -3, 5, -7, 10};
    int xValueBotC[] = {1, 3, 1 , 2,  4, 3};
    

    
    
    vector<mpq_class> xValue;
    vector<mpq_class> poly;
    
    xValue.reserve(POLYSIZE);
    poly.  reserve(POLYSIZE);
    
    mpq_class temp;
    
    for(int i = 0; i < POLYSIZE; ++i)
	{
		//xValue[i] = mpq_class(xValueTopC[i], xValueBotC[i]);
		//poly[i] = mpq_class(polyTopC[i], polyBotC[i]);
    	xValue.push_back(mpq_class(xValueTopC[i], xValueBotC[i]));
    	poly.push_back(mpq_class(polyTopC[i], polyBotC[i]));
		xValue[i].canonicalize();
		poly[i].canonicalize();
	}

    cout << "the Poly nominal = ";
    //for(int i = 0; i < POLYSIZE; ++i)
    for(int i = 0; i < (int) poly.size(); ++i)
    {
    	cout << poly[i];
    	if (i == 0)
    		cout << " ";
    	else
    		cout << "x^" << i << " ";
    }
    cout << endl;

cout << "hello world" << endl;
	
	PolynomialInterpolation p(POLYSIZE -1);
	for(int i = 0; i < POLYSIZE; ++i)
	{

		temp = evaluatePoly(xValue[i], poly);
		temp.canonicalize();
		cout << "( " << xValue[i] << ", " << temp << " )" << endl;
		//mpq_set(fValue[i], temp);
		p.addPoint(xValue[i], temp);
	}	

   	p.printMatrix();
   	p.GE();
   	p.printMatrix();
	
	cout << "testing one poly" << endl;
*/
	for(int k = 0; k < 100000; ++k)
	{
		cout << "k=" << k << "\n";
		test1Poly();
	}//test many different random poly.
	
	return 0;
}//main


/**
 * 	returns f(xValue), where f is the polynomial given by poly. poly[i] is the coeff of x^i.
 */
mpq_class evaluatePoly(mpq_class const & xValue, vector<mpq_class>  const & poly)
{
	mpq_class power, sum;
	sum = 0;
	power = 1;

	sum.canonicalize();
	power.canonicalize();

	for( int i = poly.size() - 1; i > 0; --i)
	{
		sum = (sum + poly[i])* xValue;
		sum.canonicalize();
		if ( ! isReduced(sum))
		{
			cout << "ops, not in reduced form: " << sum << endl;
			exit(1);
		}
	}//for i
	mpq_class finalSum = sum + poly[0];
	finalSum.canonicalize();
	return finalSum;
}//evaluatePoly



//checks to see gcd(r.num, r.den) = 1.
bool isReduced(mpq_class r)
{
	mpz_class a(r.get_num()), b(r.get_den());
	mpz_class c;

    while(1)
    {
		c = a%b;
		if(c==0)
		{
			return b == 1 || b == -1;
		}

		a = b;
		b = c;
    }//while
}//isReduced


/*	Randomly generate polynomials, and create the coefficient matirx. Then solves the matrix and
 *  	makes sure the returnd answer is the same as the polynomial we randomly generated.
 *
 * 	The polynomials have a max degree of MAXDEG and each coefficient is negative with prob. PNEG..
 * 		and coefficients are limited by MAXCOEF.
 */
void test1Poly()
{

	int degree = rand()% MAXDEG + 1;
	cout << "degree=" << degree << "\n";
	vector<mpq_class> poly;

	//make a poly of finite degree.
	for(int i = 0; i <= degree; ++i)
	{
		int num = rand() % MAXCOEF;
		if ( rand() < RAND_MAX * PNEG)
			num = -1 * num;
		poly.push_back(mpq_class(num, rand() % MAXCOEF +1));
		poly[i].canonicalize();
	}//for i
	;

	PolynomialInterpolation p(degree);
	PolynomialInterpolation pCopy(degree); //pCopy is printed as the origional matrix if there is an error.
	vector<mpq_class> allX(degree); //keep track of points added so far.

	//insert degree+1 many unique points (x, f(x)).
	for(int i = 0; i <= degree; ++i)
	{
		mpq_class x;
		bool notFound = true;

		//make sure we did not already picked/inserted this x.
		while ( notFound)
		{
			x = mpq_class(rand() % MAXCOEF, rand() % MAXCOEF+1);
			x.canonicalize();
			for(int j = 0; j < (int) allX.size(); ++j)
				if (x == allX[j])
					notFound = false;
			if ( notFound == false)
				notFound = true; //try again.
			else
			{
				allX.push_back(x);
				break; //this x works.
			}//else
		}//while

		//make the coefficient negative.
		if ( rand() < RAND_MAX * PNEG)
			x = x * (-1);

		mpq_class temp = evaluatePoly(x, poly);
		temp.canonicalize();
		p.addPoint(x, temp);
	}//add some points.

	pCopy = p;
	p.GE();
	vector<mpq_class> ans = p.getSolutionVector();
	
	//check answer has the correct size/degree
	if (poly.size() != ans.size() )
	{
		cerr << "answer != poly" << endl;
		exit(1);
	}

	//check answer has the correct coefficients.
	for(int i = 0; i < (int) poly.size(); ++i)
		if ( poly[i] != ans[i])
		{
			cerr << "i=" << i << ", " << poly[i] << " != " << ans[i] << endl;
			for(int k = 0; k < (int) poly.size(); ++k)
				cout << " " << poly[k] << " == " << ans[k] << "\n";
			p.printMatrix();
			cout << "\n" << endl;
			pCopy.printMatrix();
			exit(1);
		}
}//test1Poly()
