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

#endif

