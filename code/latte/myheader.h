/* This is a -*- C++ -*- header file. */
#ifndef MYHEADER__H
#define MYHEADER__H
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>

#include "latte_ntl.h"

using namespace std;

#ifdef SUN
  #define GIVETIME times(&tms_buf);z=tms_buf.tms_utime;z=z/CLK_TCK;
  #include <sys/param.h>
  #include <sys/types.h>
  #include <time.h>
  #include <sys/times.h>
#endif

#ifdef LINUX
  #define GIVETIME z=getcputime();
#endif

typedef ZZ Integer;

struct listVector {
  vec_ZZ first;
  struct listVector *rest;
};

typedef struct rationalVector {
  vec_ZZ enumerator;
  vec_ZZ denominator;
} rationalVector;

typedef struct listRationalVector {
  rationalVector* first;
  struct listRationalVector *rest;
} listRationalVector;

typedef struct listCone {
  int coefficient;
  ZZ determinant;
  rationalVector* vertex;
  listVector *rays;
  listVector *facets;
  listVector *latticePoints;
  struct listCone *rest;
} listCone;

typedef struct inputFile {
  ifstream handle;
  char char_buffer;
  char token[256];
} inputFile;

#endif

