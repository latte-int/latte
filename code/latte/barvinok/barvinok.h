/*********************************************************** -*- C++ -*-
  Author: Ruriko Yoshida
  July 24th, 2002
  Update: Febrary 3rd, 2003
  This is a program for Barvinok's decomposition of cones.
  This is a class file.

************************************************************************/
#ifndef BARVINOK__H
#define BARVINOK__H
#include <list>

#include <fstream>
#include <cstdlib>
#include <cstring>
#include <math.h>
#include <algorithm>
#include <time.h>
#include "../myheader.h"
#include "../ramon.h"
#include "Cone.h"
#include "../RudyResNTL.h"
#include "../flags.h"
using namespace std;

/* Do a signed decomposition, modulo lower-dimensional cones, of the
   cone spanned by B, until the determinants of all cones are at most
   MAX_DETERMINANT.  Return the cones in UNI. */ 
int barvinok(mat_ZZ &B, list< PtrCone > &Uni, int &numOfUniCones,
	     int max_determinant = 1);

int barvinok_Single (mat_ZZ &, int &, Single_Cone_Parameters *, rationalVector *vertex);

int barvinok_DFS(Cone *, Single_Cone_Parameters *);


#endif
