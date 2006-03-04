// This is a -*- C++ -*- header file.

#ifndef DUAL_H
#define DUAL_H

#include "cone.h"

/* Dualize the polyhedral cones in the list. */
listCone* dualizeCones(listCone*, int);

/* Fill the slots `determinant' and `facets' of CONE.  The facet
   vectors are *not* made primitive.  The facet matrix is
   an nonnegative integer multiple of -transpose(inverse(rays)).
*/
void computeDetAndFacetsOfSimplicialCone(listCone *cone, int numOfVars);

/* Destructively dualize the simplicial cones in the list.
   The ray vectors of the resulting cones are made primitive.
*/
listCone* dualizeBackCones(listCone*, int);

#endif
