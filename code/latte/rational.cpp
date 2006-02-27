/* ----------------------------------------------------------------- */
/*                                                                   */
/* LattE (Lattice Point Enumeration)                                 */
/*                                                                   */
/* Functions for handling rational vectors.                          */
/*                                                                   */
/* Author     : Raymond Hemmecke, modified by R. Yoshida             */
/*                                                                   */
/* Created    : 13-AUG-02                                            */
/* Last Update: 25-JUL-05                                            */
/*                                                                   */
/* ----------------------------------------------------------------- */
#include <stdlib.h>
#include "myheader.h"
#include "print.h"
#include "ramon.h"
#include "barvinok/Cone.h"

#include <fstream>
#include <cstdlib>
#include <cstring>
#include <string>
#include "rational.h"
/* ----------------------------------------------------------------- */
rationalVector* createRationalVector(int numOfVars) {
  int i;
  vec_ZZ x,y;
  rationalVector* z;

  x = createVector(numOfVars+1);
//    if (x==0) exit(0);
  y = createVector(numOfVars+1);
//    if (y==0) exit(0);
  z = new rationalVector;
//    z = (rationalVector*)malloc(sizeof(rationalVector));
//    if (z==0) exit(0);

  for (i=0; i<numOfVars+1; i++) {
    x[i]=0;
    y[i]=1;
  }

  z->enumerator=x;
  z->denominator=y;

  return (z);
}
/* ----------------------------------------------------------------- */
rationalVector** createArrayRationalVector(int numOfVectors) {
  rationalVector** w;

  w = new rationalVector*[numOfVectors+1];

//    w = (rationalVector**)calloc(sizeof(rationalVector*),(numOfVectors+1));
  if (w==0) exit(0);
  return (w);
}
/* ----------------------------------------------------------------- */
rationalVector* normalizeRationalVector(rationalVector *z, int numOfVars) {
  int i,j; 
  ZZ d,g;

  for(i=0; i<numOfVars; i++) {
    d=z->denominator[i];
    if (d>1) {
      for(j=0; j<numOfVars; j++) {
	g=GCD(d,z->denominator[j]);
	z->denominator[j]=z->denominator[j]/g;
	g=d/g;
	z->enumerator[j]=z->enumerator[j]*g;
      } 
    }
  }
  return (z);
}
/* ----------------------------------------------------------------- */
rationalVector* addRationalVectorsWithUpperBoundOne(rationalVector *x, 
						    rationalVector *y, 
						    int numOfVars) {
  int i;
  ZZ a,b,c,d,g,m,n,s,t;
  rationalVector *z;

/* Returns x+y if x+y<1 in every component. Otherwise 0 is returned. */

/*  printf("x and y are:\n");
  printRationalVector(x,numOfVars);
  printRationalVector(y,numOfVars); */

  z=createRationalVector(numOfVars);

  for (i=0; i<numOfVars;i++) {
    a=x->enumerator[i];
    b=x->denominator[i];
    c=y->enumerator[i];
    d=y->denominator[i];
/* Compute m/n = a/b + c/d. */
    g=lcm(b,d);
    s=g/b;
    t=g/d;

    m=s*a+t*c;
    n=g;

    if (m>=n) {
      free(z);
      return (0);
    } 

    g=GCD(m,n);
    if (g!=1) {
      m=m/g;
      n=n/g;
    }
    z->enumerator[i]=m;
    z->denominator[i]=n;
  }

  return (z);
}
/* ----------------------------------------------------------------- */
rationalVector* subRationalVector(rationalVector *x, rationalVector *y, 
				  int numOfVars) {
  int i;
  ZZ a,b,c,d,g,m,n,s,t;
  rationalVector *z;

/* Returns x-y. */

  z=createRationalVector(numOfVars);

  if (x==0) {
    for (i=0; i<numOfVars; i++) {
      z->enumerator[i]=-y->enumerator[i];
      z->denominator[i]=y->denominator[i];
    }
    return (z);
  }

  if (y==0) {
    for (i=0; i<numOfVars; i++) {
      z->enumerator[i]=x->enumerator[i];
      z->denominator[i]=x->denominator[i];
    }
    return (z);
  }

  for (i=0; i<numOfVars;i++) {
    a=x->enumerator[i];
    b=x->denominator[i];
    c=y->enumerator[i];
    d=y->denominator[i];
/* Compute m/n = a/b - c/d. */
    g=abs(lcm(b,d));
    s=g/b;
    t=g/d;

    m=s*a-t*c;
    n=g;

    g=abs(GCD(m,n));
    if (g!=1) {
      m=m/g;
      n=n/g;
    }
    z->enumerator[i]=m;
    z->denominator[i]=n;
  }

  return (z);
}
/* ----------------------------------------------------------------- */
vec_ZZ constructRay(rationalVector* v, rationalVector* w, int numOfVars) {
  int i;
  ZZ d,g,factorV,factorW;
  rationalVector* z;

  z=createRationalVector(numOfVars);

  /* Constructing z=w-v. */

  for(i=0; i<numOfVars; i++) {
    g=lcm(v->denominator[i],w->denominator[i]);
    factorV=g/(v->denominator[i]);
    factorW=g/(w->denominator[i]);
    z->enumerator[i]=(w->enumerator[i])*factorW-(v->enumerator[i])*factorV;
    z->denominator[i]=g;

    d=GCD(z->enumerator[i],z->denominator[i]);
    if (d!=1) {
      z->enumerator[i]=(z->enumerator[i])/d;
      z->denominator[i]=(z->denominator[i])/d;
    }
    if (z->denominator[i]<0) {
      z->enumerator[i]=-(z->enumerator[i]);
      z->denominator[i]=-(z->denominator[i]);
    }
  }

  /* Removing common factors from enumerators of z. */

  g=z->enumerator[0];
  for(i=1; i<numOfVars; i++) g=GCD(g,z->enumerator[i]);
  g=abs(g);

  if (g!=1) {
    for (i=0; i<numOfVars; i++) 
      z->enumerator[i]=(z->enumerator[i])/g;
  }

  /* Normalizing rational vector z to integer vector. */

  z=normalizeRationalVector(z,numOfVars);
  z->denominator.kill();
  return (z->enumerator);
}
/* ----------------------------------------------------------------- */
vec_ZZ* subtractRowFromRow(vec_ZZ* M, int x, int y, int k, vec_ZZ* w, 
			   int numOfVars) {
  int i;
  ZZ g,factorX,factorY;

/* Subtract row x from row y using column k. */

  g=lcm(M[x][k],M[y][k]);
  factorX=g/(M[x][k]);
  factorY=g/(M[y][k]);

  for (i=0; i<numOfVars; i++) {
    M[y][i]=factorY*M[y][i]-factorX*M[x][i];
  }
  (*w)[y]=factorY*((*w)[y])-factorX*((*w)[x]);

  g=(*w)[y];
  for (i=0; i<numOfVars; i++) {
    g=GCD(g,M[y][i]);
  }
  g=abs(g);

  if (g>1) {
    for (i=0; i<numOfVars; i++)
      M[y][i]=M[y][i]/g;
    (*w)[y]=(*w)[y]/g;
  }

  return (M);
}
/* ----------------------------------------------------------------- */
rationalVector* solveLinearSystem(vec_ZZ* A, vec_ZZ rhs, int numOfEquations, 
				  int numOfVars) {
  int i,j,m;
  ZZ k;
  vec_ZZ w,z;
  vec_ZZ *matrix;
  rationalVector *lambda;

/* This function assumes that the system of linear equations has
   a unique solution! */

/*    printf("Matrix: (%d x %d)\n",numOfEquations,numOfVars); */
/*    for (i=0; i<numOfEquations; i++) printVector(A[i],numOfVars); */

/*    printf("rhs: "); */
/*    printVector(rhs,numOfEquations); */

  matrix=createArrayVector(numOfEquations);
  for (i=0; i<numOfEquations; i++) {
    matrix[i]=copyVector(A[i],numOfVars);
  }
  w=copyVector(rhs,numOfEquations);

  for (i=0; i<numOfVars; i++) {
    m=i;
    while ((m<numOfEquations) && (matrix[m][i]==0)) m++;
    if (m!=i) {
      /* Swap rows. */
      z=matrix[i];
      matrix[i]=matrix[m];
      matrix[m]=z;
      k=w[i];
      w[i]=w[m];
      w[m]=k;
    }

    for (j=0; j<numOfEquations; j++) {
      if ((j!=i) && ((matrix[j][i])!=0)) {
        matrix=subtractRowFromRow(matrix,i,j,i,&w,numOfVars);
      }
    }
  }

  for (i=0; i<numOfEquations; i++) {
    if (matrix[i][i]<0) {
      for (j=0; j<numOfEquations; j++) {
	matrix[i][j]=-matrix[i][j];
      }
      w[i]=-w[i];
    }
  }

  lambda=createRationalVector(numOfVars);
  for (i=0; i<numOfEquations; i++) {
    if ((matrix[i][i])>=0) {
      lambda->enumerator[i]=w[i];
      lambda->denominator[i]=matrix[i][i];
    } else {
      lambda->enumerator[i]=-w[i];
      lambda->denominator[i]=-matrix[i][i];
    }
    if (w[i]==0) lambda->denominator[i]=1;
  }
  w.kill();
  z.kill();
  delete [] matrix;
  return (lambda);
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
