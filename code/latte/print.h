// This is a -*- C++ -*- header file.

#ifndef PRINT_H
#define PRINT_H

#include "cone.h"

void printVector(vec_ZZ, int);
void printListVector(listVector*, int);
listVector* Changeproj(listVector* basis, listVector* equ, int numOfVars);
void printRationalVector(rationalVector*, int);
void printCone(listCone*, int);
void printListCone(listCone*, int);
void printVectorToFile(ostream &, vec_ZZ, int);
void printListVectorToFile(ostream &, listVector*, int);
void printVectorToFileWithoutBrackets(ostream &, vec_ZZ, int);
void printListVectorToFileWithoutBrackets(ostream &, listVector*, int);
void printRationalVectorToFile(ostream &, rationalVector*, int);
void printRationalVectorToFileWithoutBrackets(ostream &, rationalVector*, 
					      int);
void printConeToFile(ostream &out, listCone* cones, int numOfVars);
void printListConeToFile(char*, listCone*, int); 
void printResidueFile(char*, listCone*, int);

void
print_debug_matrix(const mat_ZZ &);
void
print_debug_vector(const vec_ZZ &);

#endif
