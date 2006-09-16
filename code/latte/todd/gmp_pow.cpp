/* gmp_pow.cpp -- Exponentiation of GMP rationals

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

#include "todd/gmp_pow.h"

mpq_class 
pow(mpq_class x, unsigned int n)
{
  mpz_class result_num;
  mpz_class result_den;

  mpz_pow_ui(result_num.get_mpz_t(), x.get_num_mpz_t(), n);
  mpz_pow_ui(result_den.get_mpz_t(), x.get_den_mpz_t(), n);

  return mpq_class(result_num,
		   result_den);
}

mpz_class 
pow(mpz_class x, unsigned int n)
{
  mpz_class result;
  mpz_pow_ui(result.get_mpz_t(), x.get_mpz_t(), n);
  return result;
}
