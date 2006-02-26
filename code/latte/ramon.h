#ifndef RAMON__H
#define RAMON__H
#include "myheader.h"

listVector* appendVectorToListVector(vec_ZZ, listVector*);
vec_ZZ createVector(int);
vec_ZZ* createArrayVector(int);
listVector* createListVector(vec_ZZ);
vec_ZZ copyVector(vec_ZZ, int);
vec_ZZ addVector(vec_ZZ, vec_ZZ, int);
vec_ZZ subVector(vec_ZZ, vec_ZZ, int);
void removeListVector( listVector* p );
vec_ZZ negativeVector(vec_ZZ, int);
int lengthListVector(listVector*);
listVector* updateBasis(listVector*, listVector*);
vec_ZZ* transformListVectorToArrayVector(listVector*, vec_ZZ*);
listVector* transformArrayVectorToListVector(vec_ZZ*, int);
int isVectorEqualToVector(vec_ZZ, vec_ZZ, int);
int isVectorInListVector(vec_ZZ, listVector*, int);
listVector* readListVector(char*,int*);
listVector* readListVectorMLP(char*,int*);

#endif

