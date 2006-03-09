// This is a -*- C++ -*- header file.

// Conversions between data types; no interesting computations here.

#ifndef CONVERT_H
#define CONVERT_H

#include "cone.h"

listVector *
transformArrayBigVectorToListVector(const mat_ZZ &A,
				    int numOfVectors,
				    int numOfVars);

/* Create a matrix whose ROWS are the ray vectors of CONE. */
mat_ZZ
createConeDecMatrix(const listCone *cone, int numOfRays, int numOfVars);

/* Create a matrix whose ROWS are the facet vectors of CONE. */
mat_ZZ
createFacetMatrix(const listCone *cone, int numOfFacets, int numOfVars);

#endif
