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

#include <string>

#include "cone.h"
#include "latte_ntl.h"

// Changes Filename!
void CheckRed(std::string &Filename, char *equ, char* max, char* nonneg, char* interior, char* dil, int dilation);
void CheckRed(char *Filename, char *equ, char* max, char* nonneg, char* interior, char* dil, int dilation);

void CheckInputFileCDDRep1(const char *InputFile);
listCone* ProjectUp(listCone* cone, int & oldNumOfVars, int & newNumOfVars, 
		    listVector *equations);
listCone* ProjectUp2(listCone* cone, int & oldNumOfVars, int & newNumOfVars, 
		     mat_ZZ AA, vec_ZZ b);
void CheckInputFileCDDRep(const char *InputFile);
void CheckInputFileCDDRep3(const char *InputFile);
void CheckInputFileCDDRep4(const char *InputFile);
void CheckInputFile(const char *InputFile);
void CheckInputFileVrep(const char *InputFile);
void CheckLength(const char * filename, char * equ);
void CheckLength2(const char * filename, char* equ);
void readLatteProblem(const char *fileName, listVector **equations,
		      listVector **inequalities, 
		      char *equationsPresent,
                      int *numOfVars, char *nonneg, char* dual,
		      char* grobner, char* Vrep);

ZZ FindBigElt(listVector* equation, int numOfVars);

int CDDstylereadLatteProblem(const char *fileName, listVector **equations,
		      listVector **inequalities, 
		      char *equationsPresent,
                      int *numOfVars, char *nonneg, char* dual,
                      char* taylor, int & degree, 
                      char* rat, int & cone_output, 
		      char* Memory_Save, char* uni, char* inthull,
		      char* grobner);


#endif
