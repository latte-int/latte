#ifndef MYHEADER__H
#define MYHEADER__H
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <NTL/vec_ZZ.h>
#include <NTL/ZZ.h>

using namespace NTL_NAMESPACE;
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

typedef vec_ZZ vector;

typedef struct listVector {
  vector first;
  struct listVector *rest;
}
listVector;

typedef struct PtrCone {
  bool sign;
  listVector *Generator;
}
PtrCone;

typedef struct rationalVector {
  vector enumerator;
  vector denominator;
} rationalVector;

typedef struct listRationalVector {
  rationalVector* first;
  struct listRationalVector *rest;
} listRationalVector;

typedef struct listCone {
  int coefficient;
  rationalVector* vertex;
  listVector *rays;
  listVector *facets;
  listVector *latticePoints;
  int numOfLatticePoints;
  struct listCone *rest;
} listCone;

typedef struct inputFile {
  ifstream handle;
  char char_buffer;
  char token[256];
} inputFile;

#endif

