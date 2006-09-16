// This is a -*- C++ -*- header file.

/* latte_gmp.h -- Interface between GMP and NTL integers

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

#ifndef LATTE_GMP_H
#define LATTE_GMP_H

#include <gmpxx.h>
#include <vector>
#include "latte_ntl.h"

typedef std::vector<mpq_class> mpq_vector;
typedef std::vector<mpz_class> mpz_vector;

/* Even though NTL is based on GMP, there are no functions that allow
   to convert to GMP.  Here are our own functions.
*/

mpz_class
convert_ZZ_to_mpz(const ZZ &zz);

ZZ
convert_mpz_to_ZZ(const mpz_class &mpz);

#endif

