// This is a -*- C++ -*- header file.

/* NTL_to_LiDIA.h -- Convert between NTL and LiDIA integers.

   Copyright 2006 Susan Margulies, Matthias Koeppe

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

#ifndef NTL_TO_LIDIA_H 
#define NTL_TO_LIDIA_H 

#include "latte_ntl.h"
#include <LiDIA/bigint.h>
#include <LiDIA/bigint_matrix.h>
#include "../ramon.h"

using namespace LiDIA;

/* NTL to LiDIA conversions */
bigint_matrix
convert_mat_ZZ_to_bigint_matrix(const mat_ZZ &);
bigint_matrix
convert_listVector_to_bigint_matrix(listVector *);
bigint*
convert_vec_ZZ_to_bigint_array(const vec_ZZ &);
void
copy_vec_ZZ_to_bigint_array(bigint *, const vec_ZZ &);

/* LiDIA to NTL conversions */
mat_ZZ
convert_bigint_matrix_to_mat_ZZ(const bigint_matrix &);

/* debug/diagnostic functions */
void
print_debug_matrix(const bigint_matrix &); 
void
print_debug_vector(const bigint *, int);

#endif

