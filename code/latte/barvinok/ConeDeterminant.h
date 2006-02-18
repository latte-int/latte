/***********************************************************************
  Author: Ruriko Yoshida
  Octorber 9th, 2002
  Update: October 9th, 2002
  This is a program computes the determinant of cones.

************************************************************************/
#ifndef CONEDETERMINANT__H
#define CONEDETERMINANT__H

#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <time.h>
#include <NTL/mat_ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_RR.h>
#include <NTL/vec_RR.h>
#include <NTL/config.h>
#include <NTL/LLL.h>
#include <NTL/HNF.h>
#include <NTL/ZZ.h>
using namespace NTL_NAMESPACE;
using namespace std;

int ComputeDet(int**, int, int);

#endif
