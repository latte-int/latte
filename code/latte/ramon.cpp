#include "myheader.h"
#include "barvinok/Cone.h"
#include <fstream.h>
#include <NTL/ZZ.h>
#include <NTL/vec_ZZ.h>
/* ----------------------------------------------------------------- */
listVector* appendVectorToListVector(vector v, listVector *REST) {
  // listVector *LIST;
  listVector * LIST = new listVector;  
  // LIST=(listVector *)malloc(sizeof(listVector));
  LIST->first = v;
  LIST->rest = REST;
  return (LIST);
}
/* ----------------------------------------------------------------- */
vector createVector(int numOfVars) {
  vector w;

  w.SetLength(numOfVars);
  return (w);
}
/* ----------------------------------------------------------------- */
vector* createArrayVector(int numOfVectors) {
  vector* w;

  w = new vector[numOfVectors+1];

//    w = (vector*)calloc(sizeof(vector*),(numOfVectors+1));
  if (w==0) exit(0);
  return (w);
}
/* ----------------------------------------------------------------- */
listVector* createListVector(vector v) {
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
vector copyVector(vector v, int numOfVars) {
  int i;
  vector w;

  w = createVector(numOfVars);
  for (i=0; i<numOfVars; i++) w[i] = v[i];
  return (w);
}
/* ----------------------------------------------------------------- */
vector addVector(vector v, vector w, int numOfVars) {
  int i;

  for (i=0; i<numOfVars; i++) v[i]=v[i]+w[i];
  return (v);
}
/* ----------------------------------------------------------------- */
vector subVector(vector v, vector w, int numOfVars) {
  int i;

  for (i=0; i<numOfVars; i++) v[i]=v[i]-w[i];
  return (v);
}
/* ----------------------------------------------------------------- */
vector negativeVector(vector v, int numOfVars) {
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
vector* transformListVectorToArrayVector(listVector *L, vector* A) {
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
listVector* transformArrayVectorToListVector(vector *A, int numOfVectors) {
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
int isVectorEqualToVector(vector v, vector w, int numOfVars) {
  int i;

//    if ((v==0) || (w==0)) return (0);
  for (i=0; i<numOfVars; i++) if (!(v[i]==w[i])) return (0);
  return (1);
}
/* ----------------------------------------------------------------- */
int isVectorInListVector(vector v, listVector* LIST, int numOfVars) {
  vector w;

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
  vector b;

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
  vector b;

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
