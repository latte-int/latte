// This is a -*- C++ -*- header file.

#ifndef EXPONENTIALSUBST_H
#define EXPONENTIALSUBST_H

#include "myheader.h"
#include "latte_gmp.h"
#include "flags.h" // for Single_Cone_Parameters

// The exponential "Memory save" mode: A class of
// Single_Cone_Parameters that immediately performs exponential
// residue calculations at each subdivided cone in the tree, don't
// store the cones.
class Exponential_Single_Cone_Parameters : public Single_Cone_Parameters {
public:
  vec_ZZ generic_vector;
  mpq_class result;
  virtual int ConsumeCone(listCone *cone);
};

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
