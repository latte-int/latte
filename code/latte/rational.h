// This is a -*- C++ -*- header file.

/* rational.cpp -- Functions for handling rational vectors.

 Copyright 2002 Ruriko Yoshida
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

#ifndef RATIONAL__H
#define RATIONAL__H

#include "latte_ntl.h"

#include <fstream>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

using namespace std;


//classes in this file.
class HugInt;
class rationalVector;
class RationalNTL;
class vec_RationalNTL;

//There are some function headers at the end of the file too.

/**
 * Brandon: observation:
 * This class takes a string (of max length  200) and converts it into a ZZ
 * The string must be positive and integer.
 * Note that this can be done in RationalNTL
 */
class HugInt
{
public:
	int* integer;
	ZZ BigInt;
	HugInt(const char* tmpString)
	{
		integer = new int[200];
		int i, j;
		for (i = 0; i < 200; i++)
			integer[i] = 0;

		int Len = 0;
		int tmpLen = strlen(tmpString);
		for (i = 0; i < tmpLen; i++)
			if (tmpString[i] != '\0')
				Len = i;
		for (i = Len; i >= 0; i--)
			if (isdigit(tmpString[Len - i]))
				integer[i] = tmpString[Len - i] - '0';

		ZZ tmp;
		for (i = 0; i <= Len; i++)
		{
			tmp = 0;
			if (integer[i] != 0)
			{
				conv(tmp, integer[i]);
				for (j = Len - 1; j > Len - 1 - i; j--)
					tmp *= 10;
				BigInt += tmp;
			}
		}
	}
	~HugInt()
	{
		BigInt.kill();
		delete[] integer;
	}
};
int ReadCDD(ifstream & in, ZZ & numerator, ZZ & denominator);

class rationalVector
{
private:
	vec_ZZ enumerator;
	vec_ZZ denominator;
	bool computed_integer_scale;
	vec_ZZ integer_scale;
	ZZ integer_scale_factor;
	void compute_integer_scale();
public:
	// Construct a zero vector.
	rationalVector(int dimension = 0);
	// Construct a rational vector from an integer vector and a scalar denominator.
	rationalVector(const vec_ZZ & numer, const ZZ & denom);
	// Construct a rational vector from two integer vectors.
	rationalVector(const vec_ZZ & numer, const vec_ZZ & denom);
	// Construct a rational vector from 1  rational vectors.
	rationalVector(const vector<RationalNTL> & rational);
	const vec_ZZ &numerators() const
	{
		return enumerator;
	}
	const vec_ZZ &denominators() const
	{
		return denominator;
	}

	void scalarMultiplication(const ZZ &numer, const ZZ & denom);
	void set_entry(int i, const ZZ &numer, const ZZ &denom, bool lazy = false)
	{
		enumerator[i] = numer;
		denominator[i] = denom;
		if (lazy)
			computed_integer_scale = false;
		else
			compute_integer_scale();
	}

	void getEntry(const int i, ZZ &numer, ZZ &denom) const
	{
		if (i < enumerator.length())
		{
			numer = enumerator[i];
			denom = denominator[i];
		}//get data
	}//getEntry

	void set_numerator(int i, const ZZ &numer)
	{
		enumerator[i] = numer;
		compute_integer_scale();
	}
	void set_denominator(int i, const ZZ &denom)
	{
		denominator[i] = denom;
		compute_integer_scale();
	}
	friend rationalVector* normalizeRationalVector(rationalVector*, int);
	/* Compute an integer vector RESULT and a SCALE_FACTOR such that
	 VEC = RESULT/SCALE_FACTOR.  Return RESULT.
	 */
	friend const vec_ZZ &scaleRationalVectorToInteger(rationalVector *vec,
			int numOfVars, ZZ &scale_factor);
	/* Bring each coordinate in VEC to canonicalized form, i.e.,
	 gcd(numerator, denominator) = 1. */
	friend void canonicalizeRationalVector(rationalVector *vec, int numOfVars);
};

/**
 * RationalNTL is a wrapper for ZZ fractions.
 * I could not find a similar class in the NTL or Latte lib (I did not want to use a rationalVector...it might be confusing for others to
 * read because the vector would only be holding one element.)
 * --Brandon July 30, 2010.
 */
class RationalNTL
{
private:
	ZZ numerator, denominator;
public:
	//a b c d e f g h i j k l m n o p q r s t u v w x y z

	// CONSTRUCTORS
	RationalNTL(); //initialize to 0
	RationalNTL(const ZZ &num, const ZZ& denom);
	RationalNTL(const ZZ &num, const int denom);
	RationalNTL(const int num, const int denom);
	RationalNTL(const string &num, const string &denom);
	RationalNTL(const string &number);


	void canonicalize(); // reduces the fraction to lowest terms, and makes
						//	the denominator positive if it can be.

	//ADDITION
	RationalNTL & add(const ZZ &num, const ZZ& denom); // adds fractions and then reduces them.
	RationalNTL & add(const RationalNTL & rationalNTL);
	RationalNTL operator+(const RationalNTL & rhs) const; //rhs stands for "right hand size: object + rhs.
	RationalNTL operator-(const RationalNTL & rhs) const;
	RationalNTL & operator+=(const RationalNTL &rhs);
	RationalNTL & operator-=(const RationalNTL &rhs);

	//DIVISION
	RationalNTL & div(const ZZ & rhs); // divides by rhs.
	RationalNTL & div(const RationalNTL & rhs);
	RationalNTL operator/(const RationalNTL & rhs) const;

	RationalNTL operator/(const ZZ & rhs) const;

	//GET FUNCTIONS
	const ZZ & getNumerator() const;
	const ZZ & getDenominator() const;

	//MULTIPLICATION
	RationalNTL & mult(const ZZ &num, const ZZ &denum); // mult. two fractions and then reduces them.
	RationalNTL & mult(const RationalNTL & rationalNTL);
	RationalNTL & mult(const ZZ rhs);
	RationalNTL & power(const long e); //(a/b)^e for positive, negative, or zero e.
	static RationalNTL power(const RationalNTL & base, long e); //(base)^e
	RationalNTL operator*(const RationalNTL & rhs) const;
	RationalNTL operator*(const ZZ & rhs) const;
	RationalNTL & operator*=(const RationalNTL & rhs);
	RationalNTL & operator*=(const ZZ & rhs);

	//OTHER
	RR to_RR() const; // converts the fraction to a float.
	string str() const; //converts the fraction to a string.
	void changeSign(); //times by -1. does not try to reduce the fraction.

	// I/O
	friend ostream& operator <<(ostream &out, const RationalNTL & rationalNTL);
	friend istream& operator >>(istream &in, RationalNTL & rationalNTL);
	static ZZ readNumber(istream &in);

	//COMPARE
	bool operator==(const RationalNTL & rhs) const;
	bool operator==(const long rhs) const;
	bool operator==(const ZZ & rhs) const;
	bool operator!=(const RationalNTL & rhs) const;
	bool operator!=(const long rhs) const;

	//ASSIGNMENT
	RationalNTL & operator=(const long rhs);
	RationalNTL & operator=(const ZZ & rhs);
	RationalNTL & operator=(const RationalNTL & rhs);

};

//this class follows the vec_T interface for NTL datatypes, but not all the functions have been implemented.
class vec_RationalNTL
{
private:
	vector<RationalNTL> vec;
public:
	//a b c d e f g h i j k l m n o p q r s t u v w x y z

	vec_RationalNTL(); // initially length 0
	~vec_RationalNTL();
	vec_RationalNTL(const vec_RationalNTL& a);	// copy constructor;

	vec_RationalNTL& operator=(const vec_RationalNTL& a);
	// assignment...performs an element-wise assignment

	static RationalNTL innerProduct(const vec_RationalNTL & v1, const vec_RationalNTL & v2); //returns <v1, v2>, where <.,.> is the std. inner product.
	void kill(); 				// release space and set to length 0
	long length() const; 		// current length
	void SetLength(long n);		// set current length to n, growing vector if necessary



	RationalNTL& operator[](long i);
	const RationalNTL& operator[](long i) const;
	// indexing operation, starting from 0.
	// The first version is applied to non-const vec_T,
	// and returns a non-const reference to a T, while the second version
	// is applied to a const vec_T and returns a const reference to a T.




	//long position(const T& a) const;
	// returns position of a in the vector, or -1 if it is not there.
	// The search is conducted from position 0 to MaxAlloc()-1 of the vector,
	// and an error is raised if the object is found at position MaxLength()
	// or higher (in which case a references an uninitialized object).
	// Note that if NTL_CLEAN_PTR flag is set, this routine takes
	// linear time, and otherwise, it takes constant time.

	//long position1(const T& a) const;
	// returns position of a in the vector, or -1 if it is not there.
	// The search is conducted from position 0 to length()-1 of the vector.
	// Note that if NTL_CLEAN_PTR flag is set, this routine takes
	// linear time, and otherwise, it takes constant time.
}; //vec_RationalNTL

rationalVector* createRationalVector(int);
rationalVector** createArrayRationalVector(int);
vec_ZZ constructRay(rationalVector*, rationalVector*, int);
rationalVector* copyRationalVector(const rationalVector *);

#endif
