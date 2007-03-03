// This is a -*- C++ -*- header file.

/* RegularTriangulation.h -- Support for regular triangulations
	       
   Copyright 2006, 2007 Matthias Koeppe

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

#ifndef REGULARTRIANGULATION_H
#define REGULARTRIANGULATION_H

/* Helper functions. */
#include <vector>
#include "latte_gmp.h"

std::vector<listVector *>
ray_array(listCone *cone);

/* Height functions for regular triangulations. */
typedef void
height_function_type(mpq_t height, const vec_ZZ &ray, void *data);

void
random_height(mpq_t height, const vec_ZZ &ray, void *data);

void
delone_height(mpq_t height, const vec_ZZ &ray, void *data);

/* Data points to an integer PERCENTAGE.
   Give this PERCENTAGE of the rays a height of 2,
   the rest a height of 1. */
void
biased_random_height(mpq_t height, const vec_ZZ &ray, void *data);

#endif
