/* IntegralHull.h -- 

   Copyright 2004 Jesus A. De Loera, David Haws, Raymond Hemmecke,
      Peter Huggins, Jeremy Tauzer, Ruriko Yoshida

   This file is part of LattE.
   
   LattE is free software; you can redistribute it and/or modify it
   under the terms of the version 2 of the GNU General Public License
   as published by the Free Software Foundation.

   LattE is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with LattE; if not, write to the Free Software Foundation,
   Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
*/

#ifndef INTEGRALHULL__H
#define INTEGRALHULL__H

#include <stdarg.h>
#include "ramon.h"
#include <list>

#include <fstream>
#include <cstdlib>
#include <cstring>
#include <math.h>
#include <algorithm>
#include <time.h>

listCone* FindRationalFunction(listCone* cones, vec_ZZ a, vec_ZZ cost, int numOfVars);
listVector* Push_Vector(listVector* head, listVector* tail, int numOfVars);
vec_ZZ SolveIP(listCone* cones, listVector* matrix,  vec_ZZ cost, int numOfVars, int SINGLE_CONE);
int CheckVertices(listVector* vertices, listVector* newVertices);
listVector* GetVertices(listCone* cones,  listVector* matrix,  listVector* hyperplane, int numOfVars, int flag);
listVector* GetHRepresentation(listVector* vertices, int numOfVars);
listVector* IntegralHull(listCone* cones, listVector* matrix, int numOfVars);
ZZ Calculate_Polytope_Width (listCone *cones,listVector *matrix,int numOfVars);
#endif


