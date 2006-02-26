#include "myheader.h"
#include "barvinok/Cone.h"
#include <fstream>
#include <NTL/ZZ.h>
#include <NTL/vec_ZZ.h>
/* ----------------------------------------------------------------- */
listVector* appendVectorToListVector(vec_ZZ v, listVector *REST) {
  // listVector *LIST;
  listVector * LIST = new listVector;  
  // LIST=(listVector *)malloc(sizeof(listVector));
  LIST->first = v;
  LIST->rest = REST;
  return (LIST);
}
/* ----------------------------------------------------------------- */
vec_ZZ createVector(int numOfVars) {
  vec_ZZ w;

  w.SetLength(numOfVars);
  return (w);
}
/* ----------------------------------------------------------------- */
vec_ZZ* createArrayVector(int numOfVectors) {
  vec_ZZ* w;

  w = new vec_ZZ[numOfVectors+1];

//    w = (vec_ZZ*)calloc(sizeof(vec_ZZ*),(numOfVectors+1));
  if (w==0) exit(0);
  return (w);
}
/* ----------------------------------------------------------------- */
listVector* createListVector(vec_ZZ v) {
  return (appendVectorToListVector(v,0));
}
/* ----------------------------------------------------------------- */
 void removeListVector( listVector* p )
        {

            if( p->rest != NULL )
            {
                listVector *oldNode = p->rest;
                p->rest = p->rest->rest;  // Bypass deleted node
                delete oldNode;
            }
        }
/* ----------------------------------------------------------------- */
vec_ZZ copyVector(vec_ZZ v, int numOfVars) {
  int i;
  vec_ZZ w;

  w = createVector(numOfVars);
  for (i=0; i<numOfVars; i++) w[i] = v[i];
  return (w);
}
/* ----------------------------------------------------------------- */
vec_ZZ addVector(vec_ZZ v, vec_ZZ w, int numOfVars) {
  int i;

  for (i=0; i<numOfVars; i++) v[i]=v[i]+w[i];
  return (v);
}
/* ----------------------------------------------------------------- */
vec_ZZ subVector(vec_ZZ v, vec_ZZ w, int numOfVars) {
  int i;

  for (i=0; i<numOfVars; i++) v[i]=v[i]-w[i];
  return (v);
}
/* ----------------------------------------------------------------- */
vec_ZZ negativeVector(vec_ZZ v, int numOfVars) {
  int i;

  for (i=0; i<numOfVars; i++) v[i]=-v[i];
  return (v);
}
/* ----------------------------------------------------------------- */
int lengthListVector(listVector* LIST) {
  int len=0;

  while (LIST) {len++; LIST = LIST->rest;}
  return (len);
}
/* ----------------------------------------------------------------- */
listVector* updateBasis(listVector *v, listVector *endBasis) {
  endBasis->rest = v;
  v->rest=0;
  endBasis = endBasis->rest;
  return (endBasis);
}
/* ----------------------------------------------------------------- */
vec_ZZ* transformListVectorToArrayVector(listVector *L, vec_ZZ* A) {
  int i;
  listVector *tmp;

  i=0;

  tmp=L;
  while (tmp) {
    A[i]=tmp->first;
    tmp=tmp->rest;
    i++;
  }

  return (A);
}
/* ----------------------------------------------------------------- */
listVector* transformArrayVectorToListVector(vec_ZZ *A, int numOfVectors) {
  int i;
  listVector *L, *endL;

  L=createListVector(createVector(1));
  endL=L;

  for (i=0; i<numOfVectors; i++) {
    endL->rest = createListVector(A[i]);
    endL = endL->rest;
  }

  return (L->rest);
}
/* ----------------------------------------------------------------- */
int isVectorEqualToVector(vec_ZZ v, vec_ZZ w, int numOfVars) {
  int i;

//    if ((v==0) || (w==0)) return (0);
  for (i=0; i<numOfVars; i++) if (!(v[i]==w[i])) return (0);
  return (1);
}
/* ----------------------------------------------------------------- */
int isVectorInListVector(vec_ZZ v, listVector* LIST, int numOfVars) {
  vec_ZZ w;

  while (LIST) {
    w = LIST->first;
    LIST = LIST->rest;
    if (isVectorEqualToVector(v,w,numOfVars)==1) return (1);
  }
  return (0);
}
/* ----------------------------------------------------------------- */
listVector* readListVector(char *fileName, int* numOfVars) {
  int i,j,numOfVectors;
  listVector *basis, *endBasis;
  vec_ZZ b;

  /* Reads numOfVars, numOfVectors, list of vectors. */

  setbuf(stdout,0);
  ifstream in(fileName);
  if(!in){
    cerr << "Cannot open input file in readListVector." << endl;
    exit(1);
  }

  in >> numOfVectors;
  in >> (*numOfVars);

  basis = createListVector(createVector(*numOfVars));
  endBasis = basis;

  for (i=0; i<numOfVectors; i++) {
    b=createVector(*numOfVars);
    for (j=0; j<(*numOfVars); j++) in >> b[j];
    endBasis->rest = createListVector(b);
    endBasis=endBasis->rest;
  }

/*  printf("List of vectors:\n");
  printf("================\n");
  printListVector(basis->rest,numOfVars); */

  return(basis->rest);
}
/* ----------------------------------------------------------------- */
listVector* readListVectorMLP(char *fileName, int* numOfVars) {
  int i,j,numOfVectors;
  listVector *basis, *endBasis;
  vec_ZZ b;

  /* Reads numOfVars, numOfVectors, list of vectors. */

  setbuf(stdout,0);

  ifstream in(fileName);
  if(!in){
    cerr << "Cannot open input file in readListVectorMLP." << endl;
    exit(1);
  }

  in >> (*numOfVars);
  in >> numOfVectors;

  basis = createListVector(createVector(*numOfVars));
  endBasis = basis;

  for (i=0; i<numOfVectors; i++) {
    b=createVector(*numOfVars);
    for (j=0; j<(*numOfVars); j++) in >> b[j];
    b=negativeVector(b,*numOfVars);
    endBasis->rest = createListVector(b);
    endBasis=endBasis->rest;
  }

/*  printf("List of vectors:\n");
  printf("================\n");
  printListVector(basis->rest,numOfVars); */

  return(basis->rest);
}
/* ----------------------------------------------------------------- */
