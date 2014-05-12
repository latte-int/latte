// This is a -*- C++ -*- header file.

/* todd-expansion.h -- Evaluation of the Todd polynomial

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

#ifndef TODD_EXPANSION_H
#define TODD_EXPANSION_H

#include <vector>
#include "latte_gmp.h"

using namespace std;

/* Return the Taylor series for exp(t) at t=0. */
mpq_vector
taylor_exponential(int order);

/* Let the Taylor series of a function f(t) be given up to some order.
   Compute the coefficients of the Taylor series 1/f(t) up to the same
   order.
*/
mpq_vector
taylor_reciprocal(const mpq_vector &taylor);

/* Compute the Taylor series for t/(1-exp(t)) at t=0. */
mpq_vector
taylor_for_todd(int order);

/* Let the Taylor series of several functions be given up to some
   (common) order.  Compute the coefficients of the Taylor series of
   their product up to the same order. */
mpq_vector
taylor_product(const vector<mpq_vector> &taylors);

/* Compute the Taylor series for
   \prod_{i=1}^d (x_i t)/(1 - exp(x_i t))
   at t=0.
   The Taylor coefficients are the evaluations of the Todd polynomials
   at x. */
mpq_vector
evaluate_todd(const mpz_vector &x);

/* Likewise, but compute the Taylor series up to t=order */
mpq_vector
evaluate_todd(const mpz_vector &x, int order);

#endif
