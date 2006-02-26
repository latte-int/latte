/***********************************************************************
  Author: Ruriko Yoshida
  July 24th, 2002
  Update: September 11th, 2002
  This is a program for Barvinok's decomposition of cones.
  This is a class file.

************************************************************************/
#ifndef BINARYSEARCHIP__H
#define BINARYSEARCHIP__H

#include <fstream>
#include <cstdlib>
#include <cstring>
#include <math.h>
#include <algorithm>
#include <time.h>

#include "myheader.h"
#include "cone.h"
#include "print.h"
#include "ramon.h"
#include "rational.h"
#include "timing.h"
#include "flags.h"
#include "vertices/cdd.h"
#include <NTL/mat_ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_RR.h>
#include <NTL/vec_RR.h>
#include <NTL/config.h>
#include <NTL/LLL.h>
#include <NTL/HNF.h>
#include <NTL/ZZ.h>

  ZZ binarySearch(listVector* matrix, listVector* ineq, vec_ZZ cost, int numOfVars, char* min);

#endif

