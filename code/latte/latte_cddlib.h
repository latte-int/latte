// This is a -*- C++ -*- header file.
   
/* latte_cddlib.h -- Interface to cddlib

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

#ifndef LATTE_CDDLIB__H
#define LATTE_CDDLIB__H

#include <gmp.h>
#include "config.h"

#ifdef HAVE_CDDLIB_C_LINKAGE
extern "C" {
#endif
#define GMPRATIONAL
#include <setoper.h>
#include <cdd.h>
#ifdef HAVE_CDDLIB_C_LINKAGE
}
#endif

#include "cone.h"

/* Ensure initialization of the library by having a global constructor
   called. */
class init_cddlib_class {
public:
  init_cddlib_class(); 
};
extern init_cddlib_class init_cddlib;

/* Check for and handle cddlib errors. */
void check_cddlib_error(dd_ErrorType error, const char *proc);

/** Conversion from LattE data structures to cddlib data structures. **/

/* Create a matrix whose rows are the RAYS.
   NUM_HOMOGENIZATION_VARS extra coordinates, each set to zero,
   are introduced in front of the ray data. */
dd_MatrixPtr
rays_to_cddlib_matrix(listVector *rays, int numOfVars,
		      int num_homogenization_vars = 1,
		      int num_extra_rows = 0);

dd_PolyhedraPtr
cone_to_cddlib_polyhedron(listCone *cone, int numOfVars);

/* MATRIX should be a CDDLIB representation of a cone.
   Return a corresponding listCone. */
listCone *
cddlib_matrix_to_cone(dd_MatrixPtr matrix);

#endif

