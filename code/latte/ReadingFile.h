#ifndef READINGFILE__H
#define READINGFILE__H
#include <list>
#include <NTL/mat_ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_RR.h>
#include <NTL/vec_RR.h>
#include <NTL/config.h>
#include <NTL/LLL.h>
#include <NTL/HNF.h>
#include <NTL/ZZ.h>

#include "myheader.h"
#include "barvinok/dec.h"
#include "barvinok/barvinok.h"
#include "barvinok/cone.h"
#include "barvinok/ConeDecom.h"
#include "barvinok/Triangulation.h"
#include "vertices/cdd.h"
#include "genFunction/maple.h"
#include "genFunction/piped.h"
#include "cone.h"
#include "ConeDeterminant.h"
#include "dual.h"
#include "RudyResNTL.h"
//  #include "jesus.h"
#include "preprocess.h"
#include "print.h"
#include "ramon.h"
#include "rational.h"
#include "timing.h"
#include "flags.h"
//#include "testing.h"
#include "IntegralHull.h"

using namespace std;

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
		      char* grobner, char* max, vector & cost, char* Vrep);

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
