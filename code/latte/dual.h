// This is a -*- C++ -*- header file.

#ifndef DUAL_H
#define DUAL_H

#include "cone.h"

/* Dualize the polyhedral cones in the list. */
listCone* dualizeCones(listCone*, int);

/* Fill the slots `determinant', `facets', and `facet_divisors' of
   CONE.  The facet vectors are made primitive.  
*/
void computeDetAndFacetsOfSimplicialCone(listCone *cone, int numOfVars);

/* Destructively dualize the simplicial cones in the list.
   The ray vectors of the resulting cones are made primitive.
*/
listCone* dualizeBackCones(listCone*, int);

/* Fill the slot `facet' of all CONES by determining which of the
   INEQUALITIES are tight at the respective vertex. */
void computeTightInequalitiesOfCones(listCone *cones,
				     listVector *inequalities,
				     int numOfVars);

#endif
