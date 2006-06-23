#include <iostream>
using namespace std;

#include <cassert>

#include "Irrational.h"
#include "dual.h"
#include "convert.h"
#include "print.h"

static ZZ
lcm(const ZZ& a, const ZZ& b)
{
  return a * ( b / GCD(a, b));
}

void
irrationalizeCone(listCone *cone, int numOfVars)
{
  if (cone->facets == NULL)
    computeDetAndFacetsOfSimplicialCone(cone, numOfVars);
  cerr << "irrationalizeCone: Warning: This code is based on the irrationalization"
       << endl
       << "described in math.CO/0603308v1, which is WRONG."
       << endl
       << "Code needs to be fixed according to math.CO/0603308v2."
       << endl;
#ifdef DEBUG_IRRATIONAL
  printCone(cone, numOfVars);
#endif
  ZZ vertex_denominator;
  vec_ZZ vertex_numerator
    = scaleRationalVectorToInteger(cone->vertex, numOfVars,
				   vertex_denominator);
#ifdef DEBUG_IRRATIONAL
  cout << "vertex = " << vertex_numerator << " / "
       << vertex_denominator << endl;
#endif
  // Compute vertex multipliers, multiplied by the determinant
  mat_ZZ dual = createFacetMatrix(cone, numOfVars, numOfVars);
#ifdef DEBUG_IRRATIONAL
  cout << "dual (-B^{-1}) = " << dual;
#endif
  vec_ZZ scaled_D_lambda = dual * vertex_numerator;
  // Round down to go to the bottom of the known stability region
  int i;
  for (i = 0; i<numOfVars; i++) {
    ZZ l = scaled_D_lambda[i];
    if (l > 0)
      l = l / vertex_denominator;
    else
      l = - ((-l) / vertex_denominator); // FIXME: Right?
    scaled_D_lambda[i] = l;
  }
#ifdef DEBUG_IRRATIONAL
  cout << "bottom multipliers: " << scaled_D_lambda << endl;
#endif
  // Now move to the center; read it as the numerator over 2*D
  ZZ D = abs(cone->determinant);
  for (i = 0; i<numOfVars; i++)
    scaled_D_lambda[i] = scaled_D_lambda[i] * 2 + 1;
  vec_ZZ center_numerator
    =  (-transpose(createConeDecMatrix(cone, numOfVars, numOfVars))
	* scaled_D_lambda);
  ZZ center_denominator = 2 * D;
#ifdef DEBUG_IRRATIONAL
  cout << "--center--> " << center_numerator << " / "
       << center_denominator << endl; 
#endif
  // The actual irrationalization.
  ZZ M = 2 * power(D, numOfVars + 1);
  vec_ZZ irrationalizer_numerator;
  irrationalizer_numerator.SetLength(numOfVars);
  ZZ irrationalizer_denominator;
  ZZ TwoM_power;
  TwoM_power = 1;
  for (i = numOfVars-1; i>=0; i--) {
    irrationalizer_numerator[i] = TwoM_power;
    TwoM_power *= (2 * M);
  }
  irrationalizer_denominator = 2 * TwoM_power;
  ZZ common_denominator = lcm(center_denominator, irrationalizer_denominator);
  // Store the new vertex
  rationalVector *new_vertex = createRationalVector(numOfVars);
  for (i = 0; i<numOfVars; i++) {
    new_vertex->enumerator[i]
      = center_numerator[i] * (common_denominator / center_denominator)
      + irrationalizer_numerator[i] * (common_denominator / irrationalizer_denominator);
    new_vertex->denominator[i] = common_denominator;
  }
#ifdef DEBUG_IRRATIONAL
  cout << "--irrationalize--> ";
  printRationalVector(new_vertex, numOfVars);
#endif
  canonicalizeRationalVector(new_vertex, numOfVars);
#ifdef DEBUG_IRRATIONAL
  cout << "--canonicalize---> ";
  printRationalVector(new_vertex, numOfVars);
#endif
  assertConesIntegerEquivalent(cone, cone->vertex, numOfVars);
  assertConesIntegerEquivalent(cone, new_vertex, numOfVars);
  delete cone->vertex;
  cone->vertex = new_vertex;
  assert(isConeIrrational(cone, numOfVars));
}

void
irrationalizeCones(listCone *cones, int numOfVars)
{
  listCone *cone;
  for (cone = cones; cone != NULL; cone=cone->rest)
    irrationalizeCone(cone, numOfVars);
}

bool
isConeIrrational(listCone *cone, int numOfVars)
{
  /* The affine hulls of the proper faces do not contain any integer
     points if and only if the scalar product of the integrally
     scaled, primitive dual basis vectors with the rational vertex is
     non-integral. */
  ZZ vertex_denominator;
  vec_ZZ vertex_numerator
    = scaleRationalVectorToInteger(cone->vertex, numOfVars,
				   vertex_denominator);
  ZZ scaled_sp;
  listVector *facet;
  for (facet = cone->facets; facet; facet = facet->rest) {
    InnerProduct(scaled_sp, vertex_numerator, facet->first);
    if (divide(scaled_sp, vertex_denominator))
      return false;
  }
  return true;
}

void
checkConeIrrational(listCone *cone, int numOfVars)
{
  if (not isConeIrrational(cone, numOfVars)) {
    static NotIrrationalException notirrational;
    throw notirrational;
  }
}

void
assertConesIntegerEquivalent(listCone *cone1, rationalVector *new_vertex,
			     int numOfVars)
{
  ZZ vertex1_denominator;
  vec_ZZ vertex1_numerator
    = scaleRationalVectorToInteger(cone1->vertex, numOfVars,
				   vertex1_denominator);
  ZZ vertex2_denominator;
  vec_ZZ vertex2_numerator
    = scaleRationalVectorToInteger(new_vertex, numOfVars,
				   vertex2_denominator);
  ZZ scaled_sp_1, scaled_sp_2;
  ZZ interval_1, interval_2;
  listVector *facet1, *facet2;
  for (facet1 = cone1->facets; facet1; facet1 = facet1->rest) {
    InnerProduct(scaled_sp_1, vertex1_numerator, facet1->first);
    InnerProduct(scaled_sp_2, vertex2_numerator, facet1->first);
    div(interval_1, scaled_sp_1, vertex1_denominator);
    div(interval_2, scaled_sp_2, vertex2_denominator);
    assert(interval_1 == interval_2);
  }
}

