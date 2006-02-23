#ifndef GENFUNCTION_PIPED_H
#define GENFUNCTION_PIPED_H

vector movePoint(vector, rationalVector*, rationalVector*, vector*, int, int);
listVector* pointsInParallelepiped(rationalVector*, listVector*, listVector*,
				   int);
listVector* pointsInParallelepipedOfUnimodularCone(rationalVector*, 
						   listVector*, int);

/* For all cones in the linked list CONES, compute their latticePoints
   slot. */
void computePointsInParallelepipeds(listCone *cones, int numOfVars);

#endif
