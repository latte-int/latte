#ifndef GENFUNCTION_PIPED_H
#define GENFUNCTION_PIPED_H

vec_ZZ movePoint(vec_ZZ, rationalVector*, rationalVector*, vec_ZZ*, int, int);
listVector* pointsInParallelepiped(rationalVector*, listVector*, listVector*,
				   int);
listVector* pointsInParallelepipedOfUnimodularCone(rationalVector*, 
						   listVector*, int);

/* For all cones in the linked list CONES, compute their latticePoints
   slot. */
void computePointsInParallelepipeds(listCone *cones, int numOfVars);

#endif
