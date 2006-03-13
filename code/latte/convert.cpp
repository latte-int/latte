#include <cassert>
#include "convert.h"
#include "ramon.h"

listVector*
transformArrayBigVectorToListVector(const mat_ZZ &A, int numOfVectors,
				    int numOfVars)
{
  int i;
  vec_ZZ v;
  listVector *L = NULL, **endL;

  endL = &L;
  for (i=0; i<numOfVectors; i++) {
    v=A[i];
    *endL = createListVector(v);
    endL = &(*endL)->rest;
  }

  return L;
}

mat_ZZ
createConeDecMatrix(const listCone *cone, int numOfRays, int numOfVars)
{
  int i;
  mat_ZZ mat;
  listVector *tmp;

  mat.SetDims(numOfRays, numOfVars);

  tmp=cone->rays;
  for (i=0; i<numOfRays; i++) {
    mat[i]=copyVector(tmp->first,numOfVars);
    tmp=tmp->rest;
  }
  //removeListVector(cone->rays);
  return (mat);
}

mat_ZZ
createFacetMatrix(const listCone *cone, int numOfFacets, int numOfVars)
{
  int i;
  mat_ZZ mat;
  listVector *tmp;

  mat.SetDims(numOfFacets, numOfVars);

  tmp=cone->facets;
  for (i=0; i<numOfFacets; i++) {
    ZZ multiplier, remainder;
    DivRem(multiplier, remainder,
	   cone->determinant, cone->facet_divisors[i]);
    assert(IsZero(remainder));
    mat[i] = copyVector(tmp->first,numOfVars) * multiplier;
    tmp=tmp->rest;
  }
  return (mat);
}

