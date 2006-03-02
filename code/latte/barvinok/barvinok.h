/*********************************************************** -*- C++ -*-
  Author: Ruriko Yoshida
  July 24th, 2002
  Update: Febrary 3rd, 2003
  This is a program for Barvinok's decomposition of cones.
  This is a class file.

************************************************************************/
#ifndef BARVINOK__H
#define BARVINOK__H

#include "myheader.h"
#include "Cone.h"
#include "flags.h"

/* Do a signed decomposition, modulo lower-dimensional cones, of the
   SIMPLICIAL cone spanned by the ROW VECTORS of B with apex at
   VERTEX, until the determinants of all cones are at most
   PARAMETERS->max_determinant.

   Call PARAMETERS->ConsumeCone() for each of the small cones.
*/ 
int
barvinok_Single(mat_ZZ B, Single_Cone_Parameters *Parameters,
		rationalVector *vertex);

/* Likewise, but the cone is given by an instance of class Cone. */
int
barvinok_DFS(Cone *cone, Single_Cone_Parameters *Parameters,
	     rationalVector *vertex);


#endif
