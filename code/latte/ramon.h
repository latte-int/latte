#ifndef RAMON__H
#define RAMON__H
#include "myheader.h"

listVector* appendVectorToListVector(vector, listVector*);
vector createVector(int);
vector* createArrayVector(int);
listVector* createListVector(vector);
vector copyVector(vector, int);
vector addVector(vector, vector, int);
vector subVector(vector, vector, int);
void removeListVector( listVector* p );
vector negativeVector(vector, int);
int lengthListVector(listVector*);
listVector* updateBasis(listVector*, listVector*);
vector* transformListVectorToArrayVector(listVector*, vector*);
listVector* transformArrayVectorToListVector(vector*, int);
int isVectorEqualToVector(vector, vector, int);
int isVectorInListVector(vector, listVector*, int);
listVector* readListVector(char*,int*);
listVector* readListVectorMLP(char*,int*);

#endif

