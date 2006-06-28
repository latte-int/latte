// This is a -*- C++ -*- header file.

#ifndef EXPONENTIALSUBST_H
#define EXPONENTIALSUBST_H

#include "myheader.h"
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
				 mpz_class &prod_ray_scalar_products,
				 const listCone *cone, int numOfVars)
  throw(NotGenericException);

mpq_vector
computeExponentialResidueWeights(const vec_ZZ &generic_vector,
				 const listCone *cone, int numOfVars)
  throw(NotGenericException);

Integer
scalar_power(const vec_ZZ &generic_vector,
	     const vec_ZZ &point,
	     int exponent);

  /* Compute the limit of the generating function for the integer points
   in CONE, whose "latticePoints" (in the fundamental parallelepiped)
   are already computed, for z -> 1, by making an exponential
   substitution z |-> exp(t lambda) where lambda is a generic vector,
   as described in [Barvinok--Pommersheim].
   The function takes the "coefficient" of the cone into consideration.
*/
mpq_class
computeExponentialResidue_Single(const vec_ZZ &lambda,
				 listCone *cone, int numOfVars);

/* Likewise, but for the whole list of CONES, summing up the
   results. */
Integer
computeExponentialResidue(listCone *cones, int numOfVars);

/* Likewise, but do simplicial triangulation and Barvinok
   decomposition, then perform it on the summands. */
Integer
decomposeAndComputeExponentialResidue(listCone *cones, 
				      Exponential_Single_Cone_Parameters &param);


#endif
