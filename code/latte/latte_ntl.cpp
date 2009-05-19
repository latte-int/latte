/* latte_ntl.cpp -- Interface to NTL

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

#include <cassert>
#include <iostream>
#include <climits>
#include "latte_ntl.h"
#include "latte_gmp.h"

void
InnerProductModulo(ZZ &result, const vec_ZZ &a, const vec_ZZ &b, const ZZ &module)
{
#if 1
  InnerProduct(result, a, b);
  rem(result, result, module);
#else
  /* very slow! */
  result = 0;
  assert(a.length() == b.length());
  int dimension = a.length();
  int i;
  ZZ ai, bi, p;
  for (i = 0; i<dimension; i++) {
    rem(ai, a[i], module);
    rem(bi, b[i], module);
    MulMod(p, ai, bi, module);
    AddMod(result, result, p, module);
  }
#endif
}

int
convert_ZZ_to_int(const ZZ &zz)
{
  mpz_class z = convert_ZZ_to_mpz(zz);
  if (abs(z) > INT_MAX) {
    std::cerr << "Numbers too large for conversion to machine integer" << std::endl;
    abort();
  }
  return mpz_get_si(z.get_mpz_t());
}
