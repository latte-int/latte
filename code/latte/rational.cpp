/* rational.cpp -- Functions for handling rational vectors.

 Copyright 2002-2005 Raymond Hemmecke, Ruriko Yoshida
 Copyright 2006 Matthias Koeppe

 This file is part of LattE.

 LattE is free software; you can redistribute it and/or modify it
 under the terms of the version 2 of the GNU General Public License
 as published by the Free Software Foundation.

 LattE is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with LattE; if not, write to the Free Software Foundation,
 Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
 */

#include <stdlib.h>
#include "print.h"
#include "ramon.h"
#include <cassert>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <string>
#include <sstream>
#include "rational.h"
/* ----------------------------------------------------------------- */
rationalVector::rationalVector(int dimension)
{
	enumerator.SetLength(dimension);
	denominator.SetLength(dimension);
	int i;
	for (i = 0; i < dimension; i++)
	{
		enumerator[i] = 0;
		denominator[i] = 1;
	}
	computed_integer_scale = false;
}

rationalVector::rationalVector(const vec_ZZ &numer, const ZZ &denom)
{
	int dimension = numer.length();
	enumerator = numer;
	denominator.SetLength(dimension);
	int i;
	for (i = 0; i < dimension; i++)
		denominator[i] = denom;
	integer_scale = numer;
	integer_scale_factor = denom;
	computed_integer_scale = true;
}

rationalVector::rationalVector(const vec_ZZ &numer, const vec_ZZ & denom)
{
	assert(numer.length() == denom.length());
	enumerator = numer;
	denominator = denom;
	computed_integer_scale = false;
}
// Construct a rational vector from 1  rational vectors.
rationalVector::rationalVector(const vector<RationalNTL> & rational)
{
	enumerator.SetLength(rational.size());
	denominator.SetLength(rational.size());

	for(size_t i = 0; i < rational.size(); ++i)
	{
		enumerator[i] = rational[i].getNumerator();
		denominator[i] = rational[i].getDenominator();
	}
	computed_integer_scale = false;
}

/* ----------------------------------------------------------------- */
/**
 * Computes new_rationalVector = old_rationalVector * numer / denom.
 * We also canonicalize the vector and find the integer scale.
 */
void rationalVector::scalarMultiplication(const ZZ &numer, const ZZ & denom)
{
	for (int i = 0; i < denominator.length(); ++i)
	{
		enumerator[i] *= numer;
		denominator[i] *= denom;
	}
	computed_integer_scale = false;
	canonicalizeRationalVector(this, enumerator.length()); //this function will also compute the integer scale.
}//scalarMultiplication


/* ----------------------------------------------------------------- */
rationalVector* createRationalVector(int numOfVars)
{
	return new rationalVector(numOfVars);
}

/* ----------------------------------------------------------------- */
rationalVector** createArrayRationalVector(int numOfVectors)
{
	rationalVector** w;

	w = new rationalVector*[numOfVectors + 1];

	if (w == 0)
		abort();
	return (w);
}

/* ----------------------------------------------------------------- */
rationalVector* normalizeRationalVector(rationalVector *z, int numOfVars)
{
	int i, j;
	ZZ d, g;

	for (i = 0; i < numOfVars; i++)
	{
		d = z->denominators()[i];
		if (d > 1)
		{
			for (j = 0; j < numOfVars; j++)
			{
				g = GCD(d, z->denominators()[j]);
				z->set_denominator(j, z->denominators()[j] / g);
				g = d / g;
				z->set_numerator(j, z->numerators()[j] * g);
			}
		}
	}
	return (z);
}

/* ------------------------------------------------------- */
RationalNTL::RationalNTL()
{
	denominator = 1; //num = 0
}

RationalNTL::RationalNTL(const ZZ &num, const ZZ& denom) :
	numerator(num), denominator(denom)
{
	canonicalize();
}

RationalNTL::RationalNTL(const ZZ &num, const int denom) :
	numerator(num)
{
	denominator = denom;
	canonicalize();
}

RationalNTL::RationalNTL(const int num, const ZZ& denom):denominator(denom)
{
	numerator = num;
	canonicalize();
}

RationalNTL::RationalNTL(const int num, const int denom)
{
	numerator = num;
	denominator = denom;
	canonicalize();
}

RationalNTL::RationalNTL(const string &num, const string &denom)
{
	numerator = to_ZZ(num.c_str());
	denominator = to_ZZ(denom.c_str());
	canonicalize();
}

RationalNTL::RationalNTL(const string &number)
{
	for (size_t i = 0; i < number.length(); ++i)
		if (number[i] == '/')
		{
			numerator = to_ZZ(number.substr(0, i).c_str());
			denominator = to_ZZ(
					number.substr(i + 1, number.length() - i - 1).c_str());
			canonicalize();
			return;
		}
	numerator = to_ZZ(number.c_str());
	denominator = 1;
}

/**
 * factors out common terms and makes the denominator positive.
 */
void RationalNTL::canonicalize()
{

	if (denominator < 0)
	{
		denominator *= -1;
		numerator *= -1;
	}//make denominator positive.

	ZZ gcd = GCD(numerator, denominator); //GCD(numerator, denominator);

	if (gcd != 1)
	{
		//cout << numerator << "/" <<  denominator << endl;
		//cout << "n:" << numerator/gcd << endl;
		numerator /= gcd;
		//cout << "*" << flush;
		//cout << "d:" << denominator/gcd << flush;
		denominator /= gcd;
		//cout << " saved." << endl;
	}//if can divide.

}//canonicalize


RationalNTL & RationalNTL::add(const ZZ &num, const ZZ& denom)
{

	numerator = numerator * denom + num * denominator;
	denominator *= denom;
	assert(denom != 0);
	canonicalize();
	return *this;
}//add


RationalNTL & RationalNTL::add(const RationalNTL & rationalNTL)
{
	return add(rationalNTL.numerator, rationalNTL.denominator);
}//add


RationalNTL RationalNTL::operator+(const RationalNTL & rhs) const
{
	RationalNTL answer(*this);
	return answer.add(rhs.numerator, rhs.denominator);
}
RationalNTL RationalNTL::operator-(const RationalNTL & rhs) const
{
	RationalNTL answer(*this);
	return answer.add(rhs.numerator * -1, rhs.denominator);
}

RationalNTL & RationalNTL::operator+=(const RationalNTL &rhs)
{
	return add(rhs.numerator, rhs.denominator);
}
RationalNTL & RationalNTL::operator-=(const RationalNTL &rhs)
{
	return add(rhs.numerator * -1, rhs.denominator);
}

RationalNTL & RationalNTL::div(const ZZ & rhs)
{
	denominator *= rhs;
	canonicalize();
	return *this;
}

RationalNTL & RationalNTL::div(const RationalNTL & rhs)
{
	numerator *= rhs.denominator;
	denominator *= rhs.numerator;
	canonicalize();
	return *this;
}

RationalNTL RationalNTL::operator/(const RationalNTL & rhs) const
{
	RationalNTL answer(*this);
	return answer.div(rhs);
}



RationalNTL RationalNTL::operator/(const ZZ & rhs) const
{
	RationalNTL answer(*this);
	return answer.div(rhs);
}

const ZZ & RationalNTL::getNumerator() const
{
	return numerator;
}

const ZZ & RationalNTL::getDenominator() const
{
	return denominator;
}

RationalNTL & RationalNTL::mult(const ZZ &num, const ZZ &denum)
{
	numerator *= num;
	denominator *= denum;
	canonicalize();
	return *this;
}

RationalNTL & RationalNTL::mult(const RationalNTL & rationalNTL)
{
	return mult(rationalNTL.numerator, rationalNTL.denominator);
}

RationalNTL & RationalNTL::mult(const ZZ rhs)
{
	numerator *= rhs;
	canonicalize();
	return *this;
}

/**
 * computes (a/b)^e for positive, negative, or zero e.
 */
RationalNTL & RationalNTL::power(const long e)
{
	if (e > 0)
	{
		numerator = NTL::power(numerator, e);
		denominator = NTL::power(denominator, e);
	}// return (a/b) ^ e
	else if (e < 0)
	{
		assert(numerator != 0);
		ZZ oldNum;
		oldNum = numerator;
		numerator = NTL::power(denominator, e * -1);
		denominator = NTL::power(oldNum, e * -1);
	} // return (b/a)^|e|, where |.| is abs. value.
	else if (e == 0)
	{
		numerator = 1;
		denominator = 1;
	} // return (a/b)^0

	canonicalize();

	return *this;
}//power

//static method. (base)^e
RationalNTL RationalNTL::power(const RationalNTL & base, long e)
{
	RationalNTL answer(base);
	return answer.power(e);
}//power


RationalNTL RationalNTL::operator*(const RationalNTL & rhs) const
{
	RationalNTL answer(*this);
	return answer.mult(rhs.numerator, rhs.denominator);
}

RationalNTL RationalNTL::operator*(const ZZ & rhs) const
{
	RationalNTL answer(*this);
	return answer.mult(rhs);
}


RationalNTL RationalNTL::operator*(const int & rhs) const
{
	RationalNTL answer(*this);
	return answer.mult(to_ZZ(rhs));
}

RationalNTL & RationalNTL::operator*=(const RationalNTL & rhs)
{
	return mult(rhs.numerator, rhs.denominator);
}

RationalNTL & RationalNTL::operator*=(const ZZ & rhs)
{
	return mult(rhs);
}

RationalNTL & RationalNTL::operator*=(const int & rhs)
{
	return mult(to_ZZ(rhs));
}


RR RationalNTL::to_RR() const
{
	return NTL::to_RR(numerator) / NTL::to_RR(denominator);
}

//converts the fraction to a string.
string RationalNTL::str() const
{
	stringstream s;
	s << *this;
	return s.str();
}

//times by -1. does not try to reduce the fraction.
void RationalNTL::changeSign()
{
	numerator *= -1;
}


bool RationalNTL::operator==(const RationalNTL & rhs) const
{
	return numerator == rhs.numerator && denominator == rhs.denominator;
}

bool RationalNTL::operator==(const long rhs) const
{
	if (denominator == 1 && numerator == rhs)
		return true;
	else
		return false;
}

bool RationalNTL::operator==(const ZZ & rhs) const
{
	if (denominator == 1 && numerator == rhs)
		return true;
	else
		return false;
}

bool RationalNTL::operator!=(const RationalNTL & rhs) const
{
	return !(*this == rhs);
}

bool RationalNTL::operator!=(const long rhs) const
{
	return !(*this == to_ZZ(rhs));
}

RationalNTL & RationalNTL::operator=(const long rhs)
{
	numerator = rhs;
	denominator = 1;
	return *this;
}

RationalNTL & RationalNTL::operator=(const ZZ & rhs)
{
	numerator = rhs;
	denominator = 1;
	return *this;
}

RationalNTL & RationalNTL::operator=(const RationalNTL & rhs)
{
	if (this == &rhs)
		return *this;
	numerator = rhs.numerator;
	denominator = rhs.denominator;
	canonicalize();
	return *this;
}

ostream& operator <<(ostream &out, const RationalNTL & rationalNTL)
{
	out << rationalNTL.numerator;

	if (rationalNTL.denominator != 1)
		out << "/" << rationalNTL.denominator;
	return out;
}

//hack. to delete.
ZZ RationalNTL::myGCD(ZZ a, ZZ b) const
{

	a = a*sign(a);
	b = b*sign(b);
	long int numTimes2;
	numTimes2 = 0;

	if ( b == 1 || a == b) return b;

	while(true)
	{
		if (IsZero(a)) return b*power2_ZZ(numTimes2);
		if (IsZero(b)) return a*power2_ZZ(numTimes2);

		if ( a % 2 == 0)
		{
			if ( b % 2 == 0)
			{
				++numTimes2;
				a /= 2;
				b /= 2;
			}
			else
				a /= 2;
		}
		else if ( b % 2 == 0)
			b /= 2;
		else
		{
			if ( a < b)
				b = (b - a)/2;
			else
				a = (a - b)/2;
		}
	}
}

/**]
 * Friend function, looks at the stream for numbers in the form
 * +-a
 * +-a/b
 *
 * Example: if the stream is "   -12/-23]blah"
 * We will set rationalNTL to 12/23 and leave the stream at "]blah"
 *
 * Example: if the stream is "a12/4" this should break the function.
 *
 * Example: if the stream is "12 / 4hello" or "12 / 4   " this should read and leave
 * the stream with "hello" or "    ".
 */
istream& operator >>(istream &in, RationalNTL & rationalNTL)
{
	//goal break the string "   - 12/-23" in  to "-12"  and "-23"

	ZZ num, denum;

	num = rationalNTL.readNumber(in);


	while ( isspace(in.peek()) )
		in.get();

	if (in.peek() == '/')
	{
		in.get();
		denum = rationalNTL.readNumber(in);

	}
	else
	{
		denum = 1;
	}

	rationalNTL = RationalNTL(num,denum);

	return in;
} //operator >>.

/**
 * Reads in integer numbers and removes white space from the stream
 *
 * Example: if the steam is "   -2b  ", then we return -2 and the steam is now "b  "
 */

ZZ RationalNTL::readNumber(istream &in)
{
	//static int count = 0;
	//++count;
	stringstream s;
	char currentChar;

	while (isspace(in.peek()))
	{
		in.get();
		//cout << "space:" << in.get() << "." << endl;
	}

	currentChar = in.get();

	assert('+' == currentChar || '-' == currentChar || isdigit(currentChar));
	assert(in.eof() == false);

	//cout << "currentChar=" << currentChar << '.' << endl;
	s << currentChar;

	while( isdigit(in.peek()) )
	{
		//currentChar = in.get();
		//cout << "digit=" << currentChar << '.' << endl;
		//s << currentChar;
		s << (char) in.get(); //must cast it to a char, otherwise you will get it as an int.
	}
	return to_ZZ(s.str().c_str());
}//readNumber

/* ----------------------------------------------------------------- */


vec_RationalNTL::vec_RationalNTL()
{
}

vec_RationalNTL::vec_RationalNTL(const vec_RationalNTL& a)
{
	vec = a.vec;
}// copy constructor;


// assignment...performs an element-wise assignment
vec_RationalNTL& vec_RationalNTL::operator=(const vec_RationalNTL& a)
{
	vec = a.vec;
	return *this;
}

vec_RationalNTL::~vec_RationalNTL()
{
	vec.clear();
}

//Static method. returns <v1, v2>, where <.,.> is the std. inner product.
RationalNTL vec_RationalNTL::innerProduct(const vec_RationalNTL & v1,
		const vec_RationalNTL & v2)
{
	RationalNTL sum;
	assert(v1.length() == v2.length());
	for (int k = 0; k < v1.length(); ++k)
		sum.add(v1.vec[k] * v2.vec[k]);
	return sum;
}//innerProduct

// release space and set to length 0
void vec_RationalNTL::kill()
{
	vec.clear();
}

// current length
long vec_RationalNTL::length() const
{
	return vec.size();
}

// set current length to n, growing vector if necessary
void vec_RationalNTL::SetLength(long n)
{
	vec.resize(n);
}

// indexing operation, applied to non-const vec_RationalNTL::, and returns a non-const reference to a RationalNTL
RationalNTL& vec_RationalNTL::operator[](long i)
{
	return vec[i];
}

// applied to a const vec_T and returns a const reference to a T.
const RationalNTL& vec_RationalNTL::operator[](long i) const
{
	return vec[i];
}


/* ----------------------------------------------------------------- */
static ZZ lcm(const ZZ& a, const ZZ& b)
{
	return a * (b / GCD(a, b));
}

/* ----------------------------------------------------------------- */

vec_ZZ constructRay(rationalVector* v, rationalVector* w, int numOfVars)
{
	ZZ v_scale, w_scale;
	const vec_ZZ &v_int = scaleRationalVectorToInteger(v, numOfVars, v_scale);
	const vec_ZZ &w_int = scaleRationalVectorToInteger(w, numOfVars, w_scale);
	vec_ZZ result;
	result.SetLength(numOfVars);
	ZZ common_scale = lcm(v_scale, w_scale);
	// result = (common_scale / w_scale) * w_int - (common_scale / v_scale) * v_int;
	ZZ w_factor, v_factor;
	div(w_factor, common_scale, w_scale);
	div(v_factor, common_scale, v_scale);
	int i;
	ZZ tw, tv;
	for (i = 0; i < numOfVars; i++)
	{
		mul(tw, w_factor, w_int[i]);
		mul(tv, v_factor, v_int[i]);
		sub(result[i], tw, tv);
	}
	/* Removing common factors */
	ZZ g = result[0];
	for (i = 1; i < numOfVars; i++)
		GCD(g, g, result[i]);
	abs(g, g);
	if (g != 1)
	{
		for (i = 0; i < numOfVars; i++)
			result[i] /= g;
	}
	return result;
}

/* ----------------------------------------------------------------- */
int ReadCDD(ifstream & in, ZZ & numerator, ZZ & denominator)
{
	/*
	 Author: Ruriko Yoshida
	 Date: December 3rd, 2002
	 Update: December 7th, 2002
	 This program reads big rationals and returns ZZs for the
	 numerator and the denominator.

	 Log:
	 December 3rd:  Start writing this code.
	 December 4th:  Debug copying the string for numerator.
	 It did not copy right.  I needed to add new memory
	 everytime, it is called.
	 March 4th, 2005: Change tmpString[200] to tmpString[2000].
	 */

	/*
	 * Jan 18, 2011: Brandon: removed everything, and used the RationalNTL class.
	 *
	 *
	int i, len;
	char* tmpString = new char[2000];
	in >> tmpString;
	//  cerr << endl;
	//  cerr << tmpString << endl;
	len = strlen(tmpString);
	int flag = 0, index = 0;
	int sign = 0;

	if (tmpString[0] == '-')
		sign = 1;

	char* t2 = new char[strlen(tmpString) + 1];
	char* s2 = new char[strlen(tmpString) + 1];
	for (i = 0; i < len + 1; i++)
	{
		s2[i] = 0;
		t2[i] = 0;
	}
	//  cerr << "s2 = " << t2 << endl;
	//  cerr << "t2 = " << t2 << endl;
	conv(denominator, 1);
	for (i = 0; i < len; i++)
		if (tmpString[i] == '/')
		{
			index = i;
			flag = 1;
		}

	//  cerr << "flag = " << flag << ", index = " << index << endl;
	//  cerr << "t2 = " << t2 << endl;

	if (flag == 1)
		strncat(t2, tmpString, index);
	else
		strcpy(t2, tmpString);

	//  cerr << "t2 = " << t2 << endl;

	HugInt x(t2); //cerr << t2 << endl;
	numerator = x.BigInt;

	//  if (abs(numerator)>10) exit(1);

	if (sign == 1)
		numerator = -numerator;
	if (flag == 1)
	{
		for (i = 0; i < (len - index); i++)
		{
			s2[i] = tmpString[index + i + 1];
		}
	}

	HugInt y(s2);
	if (flag == 1)
		denominator = y.BigInt;
	//  return 1;
	delete[] s2;
	delete[] t2;
	delete[] tmpString;
	return 1;
	*/
	string rationalString;
	in >> rationalString;
	RationalNTL rational(rationalString);
	numerator   = rational.getNumerator();
	denominator = rational.getDenominator();
	return 1;
}
/* ----------------------------------------------------------------- */

rationalVector* copyRationalVector(const rationalVector *v)
{
	return new rationalVector(*v);
}

/* ----------------------------------------------------------------- */

void rationalVector::compute_integer_scale()
{
	integer_scale_factor = 1;
	int i;
	int numOfVars = numerators().length();
	integer_scale.SetLength(numOfVars);
	for (i = 0; i < numOfVars; i++)
		integer_scale_factor = lcm(integer_scale_factor, denominators()[i]);
	for (i = 0; i < numOfVars; i++)
		integer_scale[i] = numerators()[i] * (integer_scale_factor
				/ denominators()[i]);
	computed_integer_scale = true;
}

const vec_ZZ &scaleRationalVectorToInteger(rationalVector *vec, int numOfVars,
		ZZ &scale_factor)
{
	assert(numOfVars == vec->denominators().length()
			&& numOfVars == vec->numerators().length());
	if (!vec->computed_integer_scale)
	{
		vec->compute_integer_scale();
	}
	scale_factor = vec->integer_scale_factor;
	return vec->integer_scale;
}

/* ----------------------------------------------------------------- */
void canonicalizeRationalVector(rationalVector *vec, int numOfVars)
{
	int i;
	assert(numOfVars == vec->denominators().length()
			&& numOfVars == vec->numerators().length());
	for (i = 0; i < numOfVars; i++)
	{
		ZZ g = GCD(vec->numerators()[i], vec->denominators()[i]);
		if (g != 1)
		{
			vec->enumerator[i] /= g;
			vec->denominator[i] /= g;
			vec->computed_integer_scale = false;
		}
	}
	if (!vec->computed_integer_scale)
		vec->compute_integer_scale();
}
