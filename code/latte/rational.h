/*
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
#endif

rationalVector* createRationalVector(int);
rationalVector** createArrayRationalVector(int);
rationalVector* normalizeRationalVector(rationalVector*, int);
rationalVector* addRationalVectorsWithUpperBoundOne(rationalVector*, 
						    rationalVector*, int);
rationalVector* subRationalVector(rationalVector*, rationalVector*, int);
vec_ZZ constructRay(rationalVector*, rationalVector*, int);
vec_ZZ* subtractRowFromRow(vec_ZZ*, int, int, int, vec_ZZ*, int);
rationalVector* solveLinearSystem(vec_ZZ*, vec_ZZ, int, int);
