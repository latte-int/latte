/*********************************************************** -*- C++ -*-
  Author: Ruriko Yoshida
  July 24th, 2002
  Update: Febrary 3rd, 2003
  This is a program for Barvinok's decomposition of cones.
  This is a class file.

************************************************************************/
#ifndef CONEDECOM__H
#define CONEDECOM__H

#include "myheader.h"
#include "flags.h"

/* Do a signed decomposition, modulo lower-dimensional cones, of the
   polyhedral cone spanned by the ROW VECTORS of B with apex at
   CONE_VERTEX, until the determinants of all cones are at most
   PARAMETERS->max_determinant.  (When the cone is not simplicial, we
   triangulate it first.)

   Call PARAMETERS->ConsumeCone() for each of the small cones.
*/ 
int
barvinokDecomposition_Single (const mat_ZZ &B, rationalVector *Cone_Vertex,
			      Single_Cone_Parameters *Parameters);

/* Likewise, but the cone is given as a listCone (whose "rest" we
   ignore). */ 
int
barvinokDecomposition_Single (listCone *cone,
			      Single_Cone_Parameters *Parameters);

#endif






