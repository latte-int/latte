#include <cassert>
#include "convert.h"
#include "ramon.h"

listVector*
transformArrayBigVectorToListVector(const mat_ZZ &A, int numOfVectors,
				    int numOfVars)
{
  int i;
  listVector *L = NULL, **endL;

  endL = &L;
  for (i=0; i<numOfVectors; i++) {
    *endL = createListVector(A[i]);
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

mat_ZZ
createFacetMatrix2(const listCone *cone, int numOfFacets, int numOfVars)
{
  int i;
  mat_ZZ mat;
  listVector *tmp;

  mat.SetDims(numOfFacets, numOfVars);

  tmp=cone->facets;
  for (i=0; i<numOfFacets; i++) {
    ZZ multiplier, remainder;
    DivRem(multiplier, remainder,
	   abs(cone->determinant), cone->facet_divisors[i]);
    assert(IsZero(remainder));
    mat[i] = copyVector(tmp->first,numOfVars) * multiplier;
    tmp=tmp->rest;
  }
  return (mat);
}

/* latte to NTL conversions */

/* converts a latte listVector to a mat_ZZ */ 
mat_ZZ
convert_listVector_to_mat_ZZ(listVector *list) {
   listVector *tmp_list = list;
   int rows = tmp_list->first.length();
   int cols = lengthListVector(tmp_list);
   int cur_col = 0;
   mat_ZZ m;
                                                                                
   m.SetDims(cols, rows);
                                                                                
   /* add columns as rows and then take the transpose */
   while (tmp_list) {
      m[cur_col] = copyVector(tmp_list->first, rows);
      cur_col++;
      tmp_list = tmp_list->rest;
   }
   return (transpose(m));
}

