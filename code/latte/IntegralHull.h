#ifndef INTEGRALHULL__H
#define INTEGRALHULL__H

#include <stdarg.h>
#include "myheader.h"
#include "ramon.h"
#include <list>

#include <fstream>
#include <cstdlib>
#include <cstring>
#include <math.h>
#include <algorithm>
#include <time.h>

#include <NTL/mat_ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_RR.h>
#include <NTL/vec_RR.h>
#include <NTL/config.h>
#include <NTL/LLL.h>
#include <NTL/HNF.h>
#include <NTL/ZZ.h>

using namespace std;

listCone* FindRationalFunction(listCone* cones, vector a, vector cost, int numOfVars);
listVector* Push_Vector(listVector* head, listVector* tail, int numOfVars);
vector SolveIP(listCone* cones, listVector* matrix,  vector cost, int numOfVars, int SINGLE_CONE);
int CheckVertices(listVector* vertices, listVector* newVertices);
listVector* GetVertices(listCone* cones,  listVector* matrix,  listVector* hyperplane, int numOfVars, int flag);
listVector* GetHRepresentation(listVector* vertices, int numOfVars);
listVector* IntegralHull(listCone* cones, listVector* matrix, int numOfVars);
ZZ Calculate_Polytope_Width (listCone *cones,listVector *matrix,int numOfVars);
#endif


