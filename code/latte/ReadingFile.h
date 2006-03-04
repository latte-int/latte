#ifndef READINGFILE__H
#define READINGFILE__H

#include "myheader.h"
#include "latte_ntl.h"

void CheckRed(char* Filename, char *equ, char* max, char* nonneg, char* interior, char* dil, int dilation);
void IntVector(listVector* basis, int numOfVars);
void IntCone(listCone* cones, int numOfVars) ;
listCone* IntCone2(listCone* cones, int numOfVars);
void CheckInputFileCDDRep1(char *InputFile);
listCone* ProjectUp(listCone* cone, int & oldNumOfVars, int & newNumOfVars, 
		    listVector *equations);
listCone* ProjectUp2(listCone* cone, int & oldNumOfVars, int & newNumOfVars, 
		     mat_ZZ AA, vec_ZZ b);
void CheckInputFileCDDRep(char *InputFile);
void CheckInputFileCDDRep3(char *InputFile);
void CheckInputFileCDDRep4(char *InputFile);
void CheckInputFile(char *InputFile);
void CheckInputFileVrep(char *InputFile);
void CheckLength(char * filename, char * equ);
void CheckLength2(char * filename, char* equ);
void readLatteProblem(char *fileName, listVector **equations,
		      listVector **inequalities, 
		      char *equationsPresent,
                      int *numOfVars, char *nonneg, char* dual,
		      char* grobner, char* max, vec_ZZ & cost, char* Vrep);

ZZ FindBigElt(listVector* equation, int numOfVars);

int CDDstylereadLatteProblem(char *fileName, listVector **equations,
		      listVector **inequalities, 
		      char *equationsPresent,
                      int *numOfVars, char *nonneg, char* dual,
                      char* taylor, int & degree, 
                      char* rat, int & cone_output, int & flags,
		      char* Memory_Save, char* uni, char* inthull,
		      char* grobner);


#endif
