// This is a -*- C++ -*- header file.

/* Cone.cpp -- Barvinok's decomposition of a cone.

   Copyright 2002 Ruriko Yoshida
   Copyright 2006 Matthias Koeppe

   This file is part of LattE.
   
   LattE is free software; you can redistribute it and/or modify it
   under the terms of the version 2 of the GNU General Public License
   as published by the Free Software Foundation.

   LattE is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with LattE; if not, write to the Free Software Foundation,
   Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
*/

#ifndef CONE__H
#define CONE__H

#include "latte_ntl.h"

vec_ZZ ComputeOmega( const mat_ZZ & B, const mat_ZZ &Dual,
		     long m, int x, int y);
vec_ZZ ComputeOmega_2(mat_ZZ &B, long m);

#endif
