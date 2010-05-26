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
#include "PolynomialInt.h"
#include <gmp.h>
#include <gmpxx.h>

#define MAXDEG 35 //max degree of a poly.
#define PNEG .5 //prop of a neg coeff.
#define MAXCOEF 20000 //max coeff



using namespace std;

mpq_class evaluatePoly(mpq_class const & xValue, vector<mpq_class>  const & poly);
//void evaluatePoly(mpq_t &answer, mpq_t xValue, vector<mpq_class> const &poly);


bool isReduced(mpq_class r)
{
	mpz_class a(r.get_num()), b(r.get_den());

/*	if (false && mpz_cmp(a.get_mpz_t(), b.get_mpz_t()) )
	{
		mpz_class temp(a);
		a = b;
		b = temp;
	}
*/

	//cout << "find gcd of " << a << ", " << b;

	mpz_class c;

    while(1)
    {
		c = a%b;
		if(c==0)
		{
			//cout << " = " << b << endl;
			return b == 1 || b == -1;
		}

		a = b;
		b = c;
    }




}

int gcd ( int a, int b )
{
  int c;
  while ( a != 0 ) {
     c = a; a = b%a;  b = c;
  }
  return b;
}


void test1Poly()
{

	int degree = rand()% MAXDEG + 1;
	cout << "degree=" << degree << "\n";
	vector<mpq_class> poly;
	for(int i = 0; i <= degree; ++i)
	{
		int num = rand() % MAXCOEF;
		if ( rand() < RAND_MAX * PNEG)
			num = -1 * num;
		poly.push_back(mpq_class(num, rand() % MAXCOEF +1));
		poly[i].canonicalize();
	}//make a poly.
	
	//for(int i = 0; i < (int) poly.size(); ++i)
	//	cout << poly[i] << "x^" << i << " ";
	//cout << endl;
	
	PolynomialInt p(degree);
	PolynomialInt pCopy(degree);
	vector<mpq_class> allX(degree); //keep track of points added so far.
	for(int i = 0; i <= degree; ++i)
	{
	
		mpq_class x; 
		bool notFound = true;
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
		if ( rand() < RAND_MAX * PNEG)
			x = x * (-1);
		mpq_class temp = evaluatePoly(x, poly);
		temp.canonicalize();
		p.addPoint(x, temp);
		//cout << "  (" << x << ", " << temp << ") \n";
		
	}//add some points.	

	pCopy = p;
	p.GE();
	vector<mpq_class> ans = p.getSolutionVector();
	if (poly.size() != ans.size() )
	{
		cerr << "answer != poly" << endl;
		exit(1);
	}
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



int main()
{
	srand(44);

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
	
	PolynomialInt p(POLYSIZE -1);
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
	
	for(int k = 0; k < 100000; ++k)
	{
		test1Poly();
		cout << "k=" << k << "\n";
	}
	
	
	return 0;

}//main


mpq_class evaluatePoly(mpq_class const & xValue, vector<mpq_class>  const & poly)
{
	mpq_class power, sum;
	sum = 0;
	power = 1;

	sum.canonicalize();
	power.canonicalize();

	/*
	cout << "xValue = " << xValue << endl;
	cout << "evalutatePoly::poly = ";
	for (int i = 0; i < poly.size(); ++i)
		cout << poly[i] << " ";
	cout << endl;
	*/


	for( int i = poly.size() - 1; i > 0; --i)
	{
		sum = (sum + poly[i])* xValue;
		//sum = (sum + poly[i]);
		//sum.canonicalize();
		//sum = sum * xValue;
		sum.canonicalize();
		if ( ! isReduced(sum))
		{
			cout << "ops, not in reduced form: " << sum << endl;
			exit(1);
		}

	}
	mpq_class finalSum = sum + poly[0];
	finalSum.canonicalize();
	return finalSum;
	//return sum + poly[0];
}

/*
void evaluatePoly(mpq_t &answer, mpq_t xValue, vector<mpq_t> const &poly)
{
	mpq_t power, sum;
	mpq_init(power);
	mpq_init(sum);
	mpq_set_ui(power, 1, 1);
	
	mpq_t temp;
	mpq_init(temp);
	for(unsigned  int i = 0; i < poly.size(); ++i)
	{
		mpq_mul(temp, poly[i], power);
		mpq_add(sum, sum, temp);
		mpq_mul(power, power, xValue);
	}
	mpq_set(answer, sum);
}

*/
