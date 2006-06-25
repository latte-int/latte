// This is a -*- C++ -*- header file.

#ifndef IRRATIONAL_H
#define IRRATIONAL_H

#include "cone.h"

/* Move the vertex of the simplicial CONE without changing the set of
   integer points in CONE, such that all these points are in the
   relative interior of CONE.  More strongly, every CONE_2 whose apex
   is the same as that of CONE with index(CONE_2) <= index(CONE) has
   integer points only in the relative interior, not on proper faces.
*/
void
irrationalizeCone(listCone *cone, int numOfVars);

/* Likewise, for the whole list of cones. */
void
irrationalizeCones(listCone *cones, int numOfVars);

/* Return whether CONE does not contain any integer points on (the
   affine hulls of) its proper faces.
*/
bool
isConeIrrational(listCone *cone, int numOfVars);

/* Functions can throw this exception when they discover the passed
   cone was not irrational. */
struct NotIrrationalException {};

/* Check that CONE is irrational; otherwise throw a
   NotIrrationalException exception. */
void
checkConeIrrational(listCone *cone, int numOfVars);

/* Check that CONE1 and CONE1 with NEW_VERTEX contain the same integer points. */
void
assertConesIntegerEquivalent(listCone *cone1, rationalVector *new_vertex,
			     int numOfVars, const char *message);

#endif

