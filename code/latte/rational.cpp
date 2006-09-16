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
#include "myheader.h"
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

/* ----------------------------------------------------------------- */
static ZZ
lcm(const ZZ& a, const ZZ& b)
{
  return a * ( b / GCD(a, b));
}

/* ----------------------------------------------------------------- */

vec_ZZ constructRay(rationalVector* v, rationalVector* w, int numOfVars) {
  int i;
  ZZ d,g,factorV,factorW;
  rationalVector* z;

  z=createRationalVector(numOfVars);

  /* Constructing z=w-v. */

  for(i=0; i<numOfVars; i++) {
    g=lcm(v->denominators()[i],w->denominators()[i]);
    factorV=g/(v->denominators()[i]);
    factorW=g/(w->denominators()[i]);
    z->set_entry(i, (w->numerators()[i])*factorW-(v->numerators()[i])*factorV, g);

    d=GCD(z->numerators()[i],z->denominators()[i]);
    if (d!=1) {
      z->set_entry(i, (z->numerators()[i])/d, (z->denominators()[i])/d);
    }
    if (z->denominators()[i]<0) {
      z->set_entry(i, -(z->numerators()[i]), -(z->denominators()[i]));
    }
  }

  /* Removing common factors from enumerators of z. */

  g=z->numerators()[0];
  for(i=1; i<numOfVars; i++) g=GCD(g,z->numerators()[i]);
  g=abs(g);

  if (g!=1) {
    for (i=0; i<numOfVars; i++) 
      z->set_numerator(i, (z->numerators()[i])/g);
  }

  /* Normalizing rational vector z to integer vector. */

  z=normalizeRationalVector(z,numOfVars);

  vec_ZZ result = z->numerators();
  delete z;
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
//  cout << endl;
//  cout << tmpString << endl; 
  len=strlen(tmpString);
  int flag = 0, index = 0;
  int sign = 0;

  if(tmpString[0] == '-') sign = 1;

  char* t2 = new char[strlen(tmpString) + 1];
  char* s2 = new char[strlen(tmpString) + 1];
  for (i=0;i<len+1; i++) {s2[i]=0; t2[i]=0;}
//  cout << "s2 = " << t2 << endl;
//  cout << "t2 = " << t2 << endl;
  conv(denominator, 1);
  for (i = 0; i<len; i++)
    if (tmpString[i] == '/') {
      index = i; 
      flag = 1;
    }

//  cout << "flag = " << flag << ", index = " << index << endl;
//  cout << "t2 = " << t2 << endl;

  if(flag == 1)
    strncat(t2, tmpString, index);
  else
    strcpy(t2, tmpString);

//  cout << "t2 = " << t2 << endl;

  HugInt x(t2);   //cout << t2 << endl;
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
      vec->set_entry(i, vec->numerators()[i] / g, vec->denominators()[i] / g);
    }
  }
}
