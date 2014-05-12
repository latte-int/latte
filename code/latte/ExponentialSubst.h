// This is a -*- C++ -*- header file.

/* ExponentialSubst.h -- the exponential substitution

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

#ifndef EXPONENTIALSUBST_H
#define EXPONENTIALSUBST_H

#include "latte_gmp.h"
#include "barvinok/dec.h"

// The exponential "Memory save" mode: A class of
// Single_Cone_Parameters that immediately performs exponential
// residue calculations at each subdivided cone in the tree, don't
// store the cones.
class Exponential_Single_Cone_Parameters
  : public Generic_Vector_Single_Cone_Parameters {
public:
  mpq_class result;
  Exponential_Single_Cone_Parameters() :
    result(0) {};
  Exponential_Single_Cone_Parameters(const BarvinokParameters &params) :
    Generic_Vector_Single_Cone_Parameters(params),
    result(0) {};
  virtual void InitializeComputation();
  virtual int ConsumeCone(listCone *cone);
};

/* Compute the weights w_k for the contribution of CONE
   in the counting formula
   
           \sum_{k=0}^d w_k \sum_{x\in P} <c, x>^k

   where c is the GENERIC_VECTOR, P is the set of the
   integer points in the fundament parallelepiped of CONE,
   and d is the dimension of CONE (<= NUMOFVARS).

   Also return the prod_ray_scalar_products, which might
   be useful for scaling purposes.
*/

mpq_vector /* FIXME: This version can probably go away */
computeExponentialResidueWeights(const vec_ZZ &generic_vector,
				 const listCone *cone, int numOfVars)
  throw(NotGenericException);

mpq_vector /* FIXME: This version can probably go away */
computeExponentialResidueWeights(const vec_ZZ &generic_vector,
				 mpz_class &prod_ray_scalar_products,
				 const listCone *cone, int numOfVars)
  throw(NotGenericException);

mpq_vector
computeExponentialResidueWeights(const vec_ZZ &generic_vector, const listCone *cone, int numOfVars, const vec_ZZ &linForm, int M)
throw(NotGenericException);

ZZ
scalar_power(const vec_ZZ &generic_vector,
	     const vec_ZZ &point,
	     int exponent);

/* Compute \sum_{x\in P} <c, x>^k for k = 0,...,d. */
vec_ZZ
compute_sums_of_scalar_powers(listCone *cone,
			      int numOfVars,
			      const vec_ZZ &generic_vector,
			      BarvinokParameters *params);

/* likewise, but for k=0,...,order*/
vec_ZZ
compute_sums_of_scalar_powers(listCone *cone,
			      int numOfVars,
			      const vec_ZZ &generic_vector,
			      BarvinokParameters *params, int order);

/* likewise for k=0,...,d but returns mpz_vector */
mpz_vector
compute_sums_of_scalar_powers_mpz(listCone *cone,
				  int numOfVars,
				  const vec_ZZ &generic_vector,
				  BarvinokParameters *params);

  /* Compute the limit of the generating function for the integer points
   in CONE, whose "latticePoints" (in the fundamental parallelepiped)
   are already computed, for z -> 1, by making an exponential
   substitution z |-> exp(t lambda) where lambda is a generic vector,
   as described in [Barvinok--Pommersheim].
   The function takes the "coefficient" of the cone into consideration.
*/
mpq_class
computeExponentialResidue_Single(const vec_ZZ &lambda,
				 listCone *cone, int numOfVars, BarvinokParameters *params);

/* likewise, but for a weighted counting by <linForm, x>^M
 * Fixme: generic_vector is not being used, we are not processing the case if linFrom is orthogonal to some vertex-ray */
mpq_class
computeExponentialResidue_Single(const vec_ZZ &generic_vector, listCone *cone, int numOfVars, BarvinokParameters *params, const vec_ZZ & linForm, int M);

/* Likewise, but for the whole list of CONES, summing up the
   results. */
ZZ
computeExponentialResidue(listCone *cones, int numOfVars, BarvinokParameters *params);

/* Likewise, but do simplicial triangulation and Barvinok
   decomposition, then perform it on the summands. */
ZZ
decomposeAndComputeExponentialResidue(listCone *cones, 
				      Exponential_Single_Cone_Parameters &param);


#endif
