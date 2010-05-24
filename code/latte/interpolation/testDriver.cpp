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



int main()
{

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
