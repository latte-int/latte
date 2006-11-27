/* ramon.h -- Helper functions

   Copyright 2002-2004 Raymond Hemmecke
   Copyright 2006 Matthias Koeppe

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

#ifndef RAMON__H
#define RAMON__H

#include "cone.h"

vec_ZZ createVector(int);
vec_ZZ* createArrayVector(int);
listVector* createListVector(vec_ZZ);
vec_ZZ copyVector(vec_ZZ, int);
vec_ZZ addVector(vec_ZZ, vec_ZZ, int);
vec_ZZ subVector(vec_ZZ, vec_ZZ, int);

/* Remove the successor of P from the list. */
void removeListVector( listVector* p );

/* Free the whole list of vectors. */
void freeListVector(listVector *p);

vec_ZZ negativeVector(vec_ZZ, int);
listVector* updateBasis(listVector*, listVector*);
int isVectorEqualToVector(vec_ZZ, vec_ZZ, int);
int isEqual(listVector *, listVector*);
int isVectorInListVector(vec_ZZ, listVector*, int);
int isVectorInListVector(const vec_ZZ &, listVector *);
listVector* readListVector(char*,int*);

#endif
