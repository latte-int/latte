#include "convert.h"
#include "ramon.h"

listVector*
transformArrayBigVectorToListVector(const mat_ZZ &A, int numOfVectors,
				    int numOfVars)
{
  int i;
  vec_ZZ v;
  listVector *L, *endL;

  v=createVector(numOfVars);
  L=createListVector(v);
  endL=L;

  for (i=0; i<numOfVectors; i++) {
    v=A[i];
    endL->rest = createListVector(v);
    endL = endL->rest;
  }

  return (L->rest);
}

mat_ZZ
createConeDecMatrix(listCone *cones, int numOfRays, int numOfVars)
{
  int i;
  mat_ZZ mat;
  listVector *tmp;

  mat.SetDims(numOfRays, numOfVars);

  tmp=cones->rays;
  for (i=0; i<numOfRays; i++) {
    mat[i]=copyVector(tmp->first,numOfVars);
    tmp=tmp->rest;
  }
  //removeListVector(cones->rays);
  return (mat);
}

