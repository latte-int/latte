/* latte_gmp.cpp -- Interface between GMP and NTL integers

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

#include "latte_gmp.h"
#include <cassert>

mpz_class
convert_ZZ_to_mpz(const ZZ &zz)
{
  mpz_class mpz;
  convert_ZZ_to_mpz(zz, mpz);
  return mpz;
}

void
convert_ZZ_to_mpz(const ZZ &zz, mpz_class &mpz)
{
  long size = NumBytes(zz);
  unsigned char *data = new unsigned char[size];
  int sig = sign(zz);
  BytesFromZZ(data, zz, size);

  mpz_import(mpz.get_mpz_t(), size, -1, 1, 1, 0, data);
  if (sig == -1)
    mpz = -mpz;
  delete[] data;
}

ZZ
convert_mpz_to_ZZ(const mpz_class &mpz)
{
  size_t count;
  int sig = sgn(mpz);
  int size = 1;
  int nail = 0;
  int numb = 8*size - nail;
  count = (mpz_sizeinbase (mpz.get_mpz_t(), 2) + numb-1) / numb;
  unsigned char *data = new unsigned char[count * size];
  mpz_export(data, &count,
	     /*order:*/ -1, size, /*endian:*/ 1, nail,
	     mpz.get_mpz_t());
  ZZ result = ZZFromBytes(data, count);
  if (sig == -1)
    result = -result;
  delete[] data;
  return result;
}

mpq_class
convert_ZZ_to_mpq(const ZZ &zz)
{
  return mpq_class(convert_ZZ_to_mpz(zz));
}

ZZ
convert_mpq_to_ZZ(mpq_t mpq)
{
  mpq_class elt(mpq);
  assert(elt.get_den() == 1);
  return convert_mpz_to_ZZ(elt.get_num());
}

ZZ
convert_mpq_to_ZZ(mpq_class elt)
{
  assert(elt.get_den() == 1);
  return convert_mpz_to_ZZ(elt.get_num());
}

mpz_vector
convert_vec_ZZ_to_mpz_vector(const vec_ZZ &vec)
{
  mpz_vector result(vec.length());
  convert_vec_ZZ_to_mpz_vector(vec, result);
  return result;   
}

void
convert_vec_ZZ_to_mpz_vector(const vec_ZZ &vec, mpz_vector &result)
{
  int j;
  for (j = 0; j<vec.length(); j++) {
    result[j] = convert_ZZ_to_mpz(vec[j]);
  }
}
