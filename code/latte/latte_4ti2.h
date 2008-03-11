// This is a -*- C++ -*- header file.

/* latte_4ti2.h -- Interface to 4ti2
	       
   Copyright 2007 Matthias Koeppe

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

#ifndef LATTE_4TI2__H
#define LATTE_4TI2__H

#include "groebner/VectorArray.h"
#include "cone.h"

/* Create a matrix whose rows are the RAYS.
   NUM_HOMOGENIZATION_VARS extra coordinates, each set to zero,
   are introduced in front of the ray data. */
_4ti2_::VectorArray *
rays_to_4ti2_VectorArray(listVector *rays, int numOfVars,
			 int num_homogenization_vars = 0,
			 int num_extra_rows = 0);

#endif
