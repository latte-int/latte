#include <fstream>
#include "ramon.h"

/* ----------------------------------------------------------------- */
listVector* appendVectorToListVector(const vec_ZZ &v, listVector *REST) {
  // listVector *LIST;
  listVector * LIST = new listVector(v);  
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
void freeListVector( listVector* p )
{
  while (p != NULL) {
    listVector *rest = p->rest;
    delete p;
    p = rest;
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
int isVectorEqualToVector(vec_ZZ v, vec_ZZ w, int numOfVars) {
  int i;

//    if ((v==0) || (w==0)) return (0);
  for (i=0; i<numOfVars; i++) if (!(v[i]==w[i])) return (0);
  return (1);
}
/* NOTE: This function changes the LIST pointer. */ 
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

/*
 * This function does NOT change the LIST pointer 
 */
int
isVectorInListVector(const vec_ZZ & v, listVector *list) {
   int len = lengthListVector(list);
   int numOfVars = v.length();
   listVector *tmp_list = list;

   while (tmp_list) {
      if (!isVectorEqualToVector(v, tmp_list->first, numOfVars)) {
         return (0);
      }
      tmp_list = tmp_list->rest;
   }
   return (1);
}

/* ----------------------------------------------------------------- */
listVector* readListVector(char *fileName, int* numOfVars) {
  int i,j,numOfVectors;
  listVector *basis, *endBasis;
  vec_ZZ b;

  /* Reads numOfVars, numOfVectors, list of vectors. */

  //setbuf(stdout,0);
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
/*
 * Testing whether or not two lists contain exactly the same
 * points. This is a debugging diagnostic function.
 */
int
isEqual(listVector *first, listVector *second) {
   int first_len = lengthListVector(first);
   int second_len = lengthListVector(second);
   int numOfVars = first->first.length();
   listVector *tmp_list = first;
    
   if (first_len != second_len) {
      return (0);
   }
  
   /*
    * since the two lists have the same length, if every point in one list
    * is contained in the other list, the two lists must be equal.
    */ 
   while (tmp_list) {
      if (!isVectorInListVector(tmp_list->first, second)) {
         return (0);
      }
      tmp_list = tmp_list->rest;
   }
   return (1);
}
