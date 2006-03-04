/*********************************************************** -*- c++ -*-
  Author: Ruriko Yoshida
  July 24th, 2002
  Update: October 4th, 2002
  This is a program for Barvinok's decomposition of cones.
  This is a class file.

************************************************************************/
#ifndef CONE__H
#define CONE__H

#include "latte_ntl.h"

vec_ZZ ComputeOmega( const mat_ZZ & B, const mat_ZZ &Dual,
		     long m, int x, int y);
vec_ZZ ComputeOmega_2(mat_ZZ &B, long m);

#endif
