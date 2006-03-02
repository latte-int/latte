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

typedef ZZ Integer;

struct listVector {
  vec_ZZ first;
  struct listVector *rest;
};

typedef struct rationalVector {
  vec_ZZ enumerator;
  vec_ZZ denominator;
} rationalVector;

typedef struct listCone {
  int coefficient;
  ZZ determinant;		// with sign
  rationalVector* vertex;
  listVector *rays;
  listVector *facets;
  listVector *latticePoints;
  struct listCone *rest;
} listCone;

#endif

