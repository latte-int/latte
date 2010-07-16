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
#include "rational.h"
/* ----------------------------------------------------------------- */
rationalVector::rationalVector(int dimension)
{
  enumerator.SetLength(dimension);
  denominator.SetLength(dimension);
  int i;
  for (i=0; i<dimension; i++) {
    enumerator[i]=0;
    denominator[i]=1;
  }
  computed_integer_scale = false;
}

rationalVector::rationalVector(const vec_ZZ &numer, const ZZ &denom)
{
  int dimension = numer.length();
  enumerator = numer;
  denominator.SetLength(dimension);
  int i;
  for (i = 0; i<dimension; i++)
    denominator[i] = denom;
  integer_scale = numer;
  integer_scale_factor = denom;
  computed_integer_scale = true;
}


/* ----------------------------------------------------------------- */
/**
 * Computes new_rationalVector = old_rationalVector * numer / denom.
 * We also canonicalize the vector and find the integer scale.
 */
void rationalVector::scalarMultiplication(const ZZ &numer, const ZZ & denom)
{
	for(int i = 0; i < denominator.length(); ++i)
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
rationalVector** createArrayRationalVector(int numOfVectors) {
  rationalVector** w;

  w = new rationalVector*[numOfVectors+1];

  if (w==0) abort();
  return (w);
}



/* ----------------------------------------------------------------- */
rationalVector* normalizeRationalVector(rationalVector *z, int numOfVars) {
  int i,j; 
  ZZ d,g;

  for(i=0; i<numOfVars; i++) {
    d=z->denominators()[i];
    if (d>1) {
      for(j=0; j<numOfVars; j++) {
	g=GCD(d,z->denominators()[j]);
	z->set_denominator(j, z->denominators()[j]/g);
	g=d/g;
	z->set_numerator(j, z->numerators()[j]*g);
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

RationalNTL::RationalNTL(const ZZ &num, const ZZ& denom): numerator(num), denominator(denom)
{
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

/**
 * factors out common terms and makes the denominator positive.
 */
void RationalNTL::canonicalize()
{
	if ( denominator < 0)
	{
		denominator *= -1;
		numerator *= -1;
	}//make denominator positive.

	ZZ gcd = GCD(numerator, denominator);
	if ( gcd != 1)
	{
		numerator /= gcd;
		denominator/= gcd;
	}//if can divide.

}//canonicalize


RationalNTL & RationalNTL::add(const ZZ &num, const ZZ& denom)
{

	numerator = numerator * denom + num * denominator;
	denominator *= denom;
	canonicalize();
	return *this;
}//add


RationalNTL & RationalNTL::add(const RationalNTL & rationalNTL)
{
	return add(rationalNTL.numerator, rationalNTL.denominator);
}//add


RationalNTL RationalNTL::operator+(const RationalNTL & rhs)
{
	RationalNTL answer(*this);
	return answer.add(rhs.numerator, rhs.denominator);
}
RationalNTL RationalNTL::operator-(const RationalNTL & rhs)
{
	RationalNTL answer(*this);
	return answer.add(rhs.numerator * -1, rhs.denominator);
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

RationalNTL  RationalNTL::operator*(const RationalNTL & rhs)
{
	RationalNTL answer(*this);
	return answer.mult(rhs.numerator, rhs.denominator);
}

RationalNTL  RationalNTL::operator*(const ZZ & rhs)
{
	RationalNTL answer(*this);
	return answer.mult(rhs, to_ZZ(1));
}

RR RationalNTL::to_RR() const
{
	return NTL::to_RR(numerator) / NTL::to_RR(denominator);
}

bool RationalNTL::operator==(const RationalNTL & rhs) const
{
	return numerator == rhs.numerator && denominator == rhs.denominator;
}

bool RationalNTL::operator!=(const RationalNTL & rhs) const
{
	return !(*this == rhs);
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
	canonicalize(); //should not be needed.
	return *this;
}

ostream& operator <<(ostream &out, const RationalNTL & rationalNTL)
{
	out << rationalNTL.numerator;

	if ( rationalNTL.denominator != 1)
		cout << "/" << rationalNTL.denominator;
	return out;
}

/* ----------------------------------------------------------------- */
static ZZ
lcm(const ZZ& a, const ZZ& b)
{
  return a * ( b / GCD(a, b));
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
  for (i = 0; i<numOfVars; i++) {
    mul(tw, w_factor, w_int[i]);
    mul(tv, v_factor, v_int[i]);
    sub(result[i], tw, tv);
  }
  /* Removing common factors */
  ZZ g = result[0];
  for (i=1; i<numOfVars; i++)
    GCD(g, g, result[i]);
  abs(g, g);
  if (g!=1) {
    for (i = 0; i<numOfVars; i++)
      result[i] /= g;
  }
  return result;
}

/* ----------------------------------------------------------------- */
int ReadCDD(ifstream & in, ZZ & numerator, ZZ & denominator) {
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
  int i, len;
  char* tmpString = new char[2000];
  in >> tmpString;
//  cerr << endl;
//  cerr << tmpString << endl; 
  len=strlen(tmpString);
  int flag = 0, index = 0;
  int sign = 0;

  if(tmpString[0] == '-') sign = 1;

  char* t2 = new char[strlen(tmpString) + 1];
  char* s2 = new char[strlen(tmpString) + 1];
  for (i=0;i<len+1; i++) {s2[i]=0; t2[i]=0;}
//  cerr << "s2 = " << t2 << endl;
//  cerr << "t2 = " << t2 << endl;
  conv(denominator, 1);
  for (i = 0; i<len; i++)
    if (tmpString[i] == '/') {
      index = i; 
      flag = 1;
    }

//  cerr << "flag = " << flag << ", index = " << index << endl;
//  cerr << "t2 = " << t2 << endl;

  if(flag == 1)
    strncat(t2, tmpString, index);
  else
    strcpy(t2, tmpString);

//  cerr << "t2 = " << t2 << endl;

  HugInt x(t2);   //cerr << t2 << endl;
  numerator = x.BigInt;

//  if (abs(numerator)>10) exit(1);

  if(sign == 1) numerator = - numerator;
  if(flag == 1) {
    for(i=0; i<(len-index); i++) {
      s2[i]=tmpString[index+i+1];
    }
  }
  
  HugInt y(s2);
  if(flag == 1)
    denominator = y.BigInt;
//  return 1;
  delete [] s2;
  delete [] t2;
  delete [] tmpString;
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
    for (i = 0; i<numOfVars; i++)
      integer_scale_factor = lcm(integer_scale_factor, denominators()[i]);
    for (i = 0; i<numOfVars; i++)
      integer_scale[i] = numerators()[i] * (integer_scale_factor / denominators()[i]);
    computed_integer_scale = true;
}

const vec_ZZ &scaleRationalVectorToInteger(rationalVector *vec,
					   int numOfVars,
					   ZZ &scale_factor)
{
  assert(numOfVars == vec->denominators().length()
	 && numOfVars == vec->numerators().length());
  if (!vec->computed_integer_scale) {
    vec->compute_integer_scale();
  }
  scale_factor = vec->integer_scale_factor;
  return vec->integer_scale;
}

/* ----------------------------------------------------------------- */
void canonicalizeRationalVector(rationalVector *vec,
				int numOfVars)
{
  int i;
  assert(numOfVars == vec->denominators().length()
	 && numOfVars == vec->numerators().length());
  for (i = 0; i<numOfVars; i++) {
    ZZ g = GCD(vec->numerators()[i], vec->denominators()[i]);
    if (g != 1) {
      vec->enumerator[i] /= g;
      vec->denominator[i] /= g;
      vec->computed_integer_scale = false;
    }
  }
  if (!vec->computed_integer_scale)
    vec->compute_integer_scale();
}
