#ifndef GROBNER__H
#define GROBNER__H

#include "myheader.h"
#include "cone.h"
#include "print.h"
#include "ramon.h"
#include "ReadingFile.h"
#include <NTL/ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include <NTL/LLL.h>

listVector* Grobner(listVector *equations, 
		    listVector *inequalities, vector **generators,
		    int *numOfVars, listVector **newVec, int & oldnumber,
		    int bignum);
void SolveGrobner(char * filename, char * nonneg, char * dualApproach,
		  char * grobner, char * equationsPresent, char* cdd);
void CheckGrobner(char* filename, char * cdd);
#endif

