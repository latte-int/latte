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
using namespace std;

class HugInt{
public:
  int* integer;
  ZZ BigInt;
  HugInt(const char* tmpString) {
    integer = new int[200];
    int i,j;
    for(i=0; i<200; i++) integer[i]=0;
    
    int Len = 0;
    int tmpLen=strlen(tmpString);
    for(i=0; i<tmpLen; i++)
      if(tmpString[i] != '\0') Len = i;
    for(i = Len; i >= 0; i--)
      if(isdigit(tmpString[Len-i])) integer[i] = tmpString[Len-i] - '0';

    ZZ tmp;           
    for(i = 0; i <= Len; i++) {
      tmp = 0;
      if(integer[i] != 0) { 
	conv(tmp, integer[i]);
	for(j=Len-1; j>Len-1-i; j--) tmp *= 10;
	BigInt += tmp;
      }
    }
  }
  ~HugInt(){
   BigInt.kill();
   delete [] integer;
   }
};
int ReadCDD(ifstream & in, ZZ & numerator, ZZ & denominator);

class rationalVector {
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
  rationalVector(const vec_ZZ &numer, const ZZ &denom); 
  const vec_ZZ &numerators() const { return enumerator; }
  const vec_ZZ &denominators() const { return denominator; }
  void set_entry(int i, const ZZ &numer, const ZZ &denom, bool lazy=false) {
    enumerator[i] = numer;
    denominator[i] = denom;
    if (lazy)
      computed_integer_scale = false;
    else
      compute_integer_scale();
  }

  void getEntry(const int i , ZZ &numer, ZZ &denom)
  {
	if ( i < enumerator.length())
	{
		numer = enumerator[i];
		denom = denominator[i];
	}//get data
  }//getEntry

  void set_numerator(int i, const ZZ &numer) {
    enumerator[i] = numer;
    compute_integer_scale();
  }
  void set_denominator(int i, const ZZ &denom) {
    denominator[i] = denom;
    compute_integer_scale();
  }
  friend rationalVector* normalizeRationalVector(rationalVector*, int);
  /* Compute an integer vector RESULT and a SCALE_FACTOR such that
     VEC = RESULT/SCALE_FACTOR.  Return RESULT.
  */
  friend const vec_ZZ &scaleRationalVectorToInteger(rationalVector *vec,
						    int numOfVars,
						    ZZ &scale_factor);
  /* Bring each coordinate in VEC to canonicalized form, i.e.,
     gcd(numerator, denominator) = 1. */
  friend void canonicalizeRationalVector(rationalVector *vec,
					 int numOfVars);
};

rationalVector* createRationalVector(int);
rationalVector** createArrayRationalVector(int);
vec_ZZ constructRay(rationalVector*, rationalVector*, int);
rationalVector* copyRationalVector(const rationalVector *);

#endif
