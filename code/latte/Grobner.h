/* Grobner.h -- computing a Gr\"obner basis of a toric ideal.

   Copyright 2003 Raymond Hemmecke
   Copyright 2003 Ruriko Yoshida

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

#ifndef GROBNER__H
#define GROBNER__H

#include "cone.h"
#include "print.h"
#include "ramon.h"
#include "ReadingFile.h"

listVector* Grobner(listVector *equations, 
		    listVector *inequalities, vec_ZZ **generators,
		    int *numOfVars, listVector **newVec, int & oldnumber,
		    int bignum);
ZZ SolveGrobner(const char * filename, char * nonneg, char * dualApproach,
		  char * grobner, char * equationsPresent, char* cdd);
void CheckGrobner(const char* filename, char * cdd);
#endif

