/* This is a -*- C++ -*- header file.
   
  Author: Ruriko Yoshida
  Date: December 3rd, 2002
  Update: December 4th, 2002
  This program reads big rationals and returns ZZs for the
  numerator and the denominator.

  Log:
     December 3rd:  Start writing this code.
     December 4th:  Debug copying the string for numerator.
                    It did not copy right.  I needed to add new memory
                    everytime, it is called.
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
  
public:
  rationalVector(int dimension = 0);
  const vec_ZZ &numerators() { return enumerator; }
  const vec_ZZ &denominators() { return denominator; }
  void set_entry(int i, const ZZ &numer, const ZZ &denom) {
    enumerator[i] = numer;
    denominator[i] = denom;
    computed_integer_scale = false;
  }
  void set_numerator(int i, const ZZ &numer) {
    enumerator[i] = numer;
    computed_integer_scale = false;
  }
  void set_denominator(int i, const ZZ &denom) {
    denominator[i] = denom;
    computed_integer_scale = false;
  }
  friend rationalVector* normalizeRationalVector(rationalVector*, int);
};

rationalVector* createRationalVector(int);
rationalVector** createArrayRationalVector(int);
rationalVector* addRationalVectorsWithUpperBoundOne(rationalVector*, 
						    rationalVector*, int);
rationalVector* subRationalVector(rationalVector*, rationalVector*, int);
rationalVector* addRationalVector(rationalVector*, rationalVector*, int);
vec_ZZ constructRay(rationalVector*, rationalVector*, int);
rationalVector* copyRationalVector(const rationalVector *);

/* Compute an integer vector RESULT and a SCALE_FACTOR such that
   VEC = RESULT/SCALE_FACTOR.  Return RESULT.
*/
vec_ZZ scaleRationalVectorToInteger(const rationalVector *vec,
				    int numOfVars,
				    ZZ &scale_factor);

/* Bring each coordinate in VEC to canonicalized form, i.e.,
   gcd(numerator, denominator) = 1. */
void canonicalizeRationalVector(rationalVector *vec,
				int numOfVars);

#endif
