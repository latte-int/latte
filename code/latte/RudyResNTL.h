#ifndef RUDYRESNTL__H
#define RUDYRESNTL__H

#include "myheader.h"
#include "cone.h"
#include "barvinok/dec.h" // for Standard_Single_Cone_Parameters

void ResidueFunction(listCone* cones, int numOfVars, int print_flag, int degree, int output_cone);

// Returns -1 if a Dot Product is zero in the denominator, otherwise 1 if ok
int
ResidueFunction_Single_Cone (listCone *cones,
			     Standard_Single_Cone_Parameters *Residue_Parameters); 

#endif

