// This is a -*- C++ -*- header file.

// Conversions between data types; no interesting computations here.

#ifndef CONVERT_H
#define CONVERT_H

#include "myheader.h"

listVector *
transformArrayBigVectorToListVector(const mat_ZZ &A,
				    int numOfVectors,
				    int numOfVars);

mat_ZZ
createConeDecMatrix(const listCone *cone, int numOfRays, int numOfVars);

#endif
