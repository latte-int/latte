// This is a -*- C++ -*- header file.

#ifndef EXPONENTIALSUBST_H
#define EXPONENTIALSUBST_H

#include "myheader.h"

/* Compute the limit of the generating function for the integer points
   in CONE, whose "latticePoints" (in the fundamental parallelepiped)
   are already computed, for z -> 1, by making an exponential
   substitution z |-> exp(t lambda) where lambda is a generic vector,
   as described in [Barvinok--Pommersheim].
   The function takes the "coefficient" of the cone into consideration.
*/
Integer
computeExponentialResidue_Single(listCone *cone, int numOfVars);

/* Likewise, but for the whole list of CONES, summing up the
   results. */
Integer
computeExponentialResidue(listCone *cones, int numOfVars);

#endif
