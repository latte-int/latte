#include <iostream>
using namespace std;

#undef DEBUG_IRRATIONAL

#include <cassert>

#include "Irrational.h"
#include "dual.h"
#include "convert.h"
#include "print.h"

rationalVector *
computeConeStabilityCube(listCone *cone, int numOfVars, bool simplicial,
			 ZZ &length_numerator, ZZ &length_denominator)
{
  if (!simplicial) {
    cerr << "Computing stability region for non-simplicial cones is not implemented yet.";
    abort();
  }
  if (cone->facets == NULL)
    computeDetAndFacetsOfSimplicialCone(cone, numOfVars);

  ZZ vertex_denominator;
  vec_ZZ vertex_numerator
    = scaleRationalVectorToInteger(cone->vertex, numOfVars,
				   vertex_denominator);
#ifdef DEBUG_IRRATIONAL
  cout << "vertex = " << vertex_numerator << " / "
       << vertex_denominator << endl;
#endif
  // Compute vertex multipliers, multiplied by the determinant
  mat_ZZ dual = createFacetMatrix2(cone, numOfVars, numOfVars);
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
      l = - ((-l) / vertex_denominator);
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
  // Compute stability length
  length_numerator = 1;
  ZZ dual_1_norm, max_dual_1_norm;
  max_dual_1_norm = 0;
  for (i = 0; i<numOfVars; i++) {
    int j;
    dual_1_norm = 0;
    for (j = 0; j<numOfVars; j++)
      dual_1_norm += abs(dual[i][j]);
    if (dual_1_norm > max_dual_1_norm)
      max_dual_1_norm = dual_1_norm;
  }
  length_denominator = 2 * max_dual_1_norm;
  // Check whether at least the new vertex is OK.
  // (We do not check the length.)
  rationalVector *v = createRationalVector(numOfVars);
  for (i = 0; i<numOfVars; i++) {
    v->enumerator[i] = center_numerator[i];
    v->denominator[i] = center_denominator;
  }
  assertConesIntegerEquivalent(cone, v, numOfVars,
			       "Not integer equivalent with center");
  
  return v;
}

static ZZ
lcm(const ZZ& a, const ZZ& b)
{
  return a * ( b / GCD(a, b));
}

void
irrationalizeCone(listCone *cone, int numOfVars)
{
  ZZ center_denominator;
  vec_ZZ center_numerator;
  ZZ stability_length_numerator, stability_length_denominator;
  {
    rationalVector *stability_center
      = computeConeStabilityCube(cone, numOfVars, true,
				 stability_length_numerator,
				 stability_length_denominator);
    center_numerator =
      scaleRationalVectorToInteger(stability_center, numOfVars,
				   center_denominator);
    delete stability_center;
  }
  // The actual irrationalization.
  ZZ M;
  // Value from math.CO/0603308 v1 (wrong):
  //M = 2 * power(D, numOfVars + 1);
  // Value from v2:
  {
    int barvinok_depth
      = 1 + (int) ceil(log((double) NumBits(cone->determinant))
		       / log(1 + 1.0 / (numOfVars - 1)));
#ifdef DEBUG_IRRATIONAL
    cout << "Estimated tree depth " << barvinok_depth << endl;
#endif
    ZZ C;
    C = 0;
    int i;
    listVector *facet;
    for (facet = cone->facets; facet; facet = facet->rest) {
      for (i = 0; i<numOfVars; i++) {
	if (abs(facet->first[i]) > C)
	  C = abs(facet->first[i]);
      }
    }
    ZZ factorial;
    factorial = 1;
    for (i = 2; i<=numOfVars - 1; i++)
      factorial *= i;
    ZZ d;
    d = numOfVars;
    M = 2 * power(power(d, barvinok_depth) * C,
		  numOfVars - 1);
  }
  vec_ZZ irrationalizer_numerator;
  irrationalizer_numerator.SetLength(numOfVars);
  ZZ irrationalizer_denominator;
  ZZ TwoM_power;
  TwoM_power = 1;
  int i;
  for (i = numOfVars-1; i>=0; i--) {
    irrationalizer_numerator[i] = TwoM_power * stability_length_numerator;
    TwoM_power *= (2 * M);
  }
  irrationalizer_denominator = 2 * (TwoM_power / M) * stability_length_denominator;

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
  assertConesIntegerEquivalent(cone, cone->vertex, numOfVars,
			       "cone and cone not integer equivalent");
  assertConesIntegerEquivalent(cone, new_vertex, numOfVars,
			       "Not integer equivalent with new_vertex");
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
			     int numOfVars, const char *message)
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
  listVector *facet1;
  for (facet1 = cone1->facets; facet1; facet1 = facet1->rest) {
    InnerProduct(scaled_sp_1, vertex1_numerator, facet1->first);
    InnerProduct(scaled_sp_2, vertex2_numerator, facet1->first);
    div(interval_1, scaled_sp_1, vertex1_denominator);
    div(interval_2, scaled_sp_2, vertex2_denominator);
    if (interval_1 != interval_2) {
      cerr << message << endl;
      assert(interval_1 == interval_2);
    }
  }
}

