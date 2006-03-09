#include <iostream>
using namespace std;

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
  computeDetAndFacetsOfSimplicialCone(cone, numOfVars);
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
  for (i = 0; i<numOfVars; i++) {
    cone->vertex->enumerator[i]
      = center_numerator[i] * (common_denominator / center_denominator)
      + irrationalizer_numerator[i] * (common_denominator / irrationalizer_denominator);
    cone->vertex->denominator[i] = common_denominator;
  }
#ifdef DEBUG_IRRATIONAL
  cout << "--irrationalize--> ";
  printRationalVector(cone->vertex, numOfVars);
#endif
  canonicalizeRationalVector(cone->vertex, numOfVars);
#ifdef DEBUG_IRRATIONAL
  cout << "--canonicalize---> ";
  printRationalVector(cone->vertex, numOfVars);
#endif
}

void
irrationalizeCones(listCone *cones, int numOfVars)
{
  listCone *cone;
  for (cone = cones; cone != NULL; cone=cone->rest)
    irrationalizeCone(cone, numOfVars);
}

