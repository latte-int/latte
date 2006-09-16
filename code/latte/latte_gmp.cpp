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

mpz_class
convert_ZZ_to_mpz(const ZZ &zz)
{
  long size = NumBytes(zz);
  unsigned char *data = new unsigned char[size];
  int sig = sign(zz);
  BytesFromZZ(data, zz, size);

  mpz_class mpz;
  mpz_import(mpz.get_mpz_t(), size, -1, 1, 1, 0, data);
  if (sig == -1)
    mpz = -mpz;
  delete[] data;
  return mpz;
}

ZZ
convert_mpz_to_ZZ(const mpz_class &mpz)
{
  size_t size;
  int sig = sgn(mpz);
  void *data = mpz_export(NULL, &size,
			  -1, 1, 1, 0,
			  mpz.get_mpz_t());
  ZZ result = ZZFromBytes((unsigned char *)data, size);
  if (sig == -1)
    result = -result;
  free(data);
  return result;
}

