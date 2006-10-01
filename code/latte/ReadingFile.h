/* ReadingFile.h -- Reading file and check the input is correct.

   Copyright 2002, 2003 Raymond Hemmecke, Ruriko Yoshida

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

#ifndef READINGFILE__H
#define READINGFILE__H

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
