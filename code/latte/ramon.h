#ifndef RAMON__H
#define RAMON__H

#include "cone.h"

listVector* appendVectorToListVector(const vec_ZZ &, listVector*);
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
int lengthListVector(listVector*);
listVector* updateBasis(listVector*, listVector*);
int isVectorEqualToVector(vec_ZZ, vec_ZZ, int);
int isEqual(listVector *, listVector*);
int isVectorInListVector(vec_ZZ, listVector*, int);
int isVectorInListVector(const vec_ZZ &, listVector *);
listVector* readListVector(char*,int*);

#endif
