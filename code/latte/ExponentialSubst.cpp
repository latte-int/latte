#include "ExponentialSubst.h"
#include "todd/todd.h"
#include "todd/gmp_pow.h"

Integer
scalar_power(const vec_ZZ &generic_vector,
	     const vec_ZZ &point,
	     int exponent)
{
  Integer inner;
  InnerProduct(inner, generic_vector, point);
  return power(inner, exponent);
}

Integer
sum_of_scalar_powers(const vec_ZZ &generic_vector,
		     listVector *points,
		     int exponent)
{
  listVector *point;
  Integer result;
  result = 0;
  for (point = points; point != NULL; point = point->rest) {
    result += scalar_power(generic_vector, point->first, exponent);
  }
  return result;
}

Integer
computeExponentialResidue_Single(listCone *cone, int numOfVars)
{
  cerr << "computeExponentialResidue_Single: Not implemented." << endl;
  abort();
#if 0
  mpz_class *ray_scalar_products = new mpz_class[numOfVars];
  mpz_class prod_ray_scalar_products;
  vec_ZZ generic_vector; // FIXME
  prod_ray_scalar_products = 1;
  {
    listVector *ray;
    int k;
    for (k = 0, ray = cone->rays; ray != NULL; k++, ray = ray->rest) {
      ZZ inner;
      InnerProduct(inner, generic_vector, ray->first);
      ray_scalar_products[k] = inner;
      prod_ray_scalar_products *= ray_scalar_products[k];
    }
    assert(k == numOfVars);
  }
  int k;
  ZZ k_factorial;
  mpq_t result;
  result = 0;
  for (k = 0, k_factorial = 1; k<=numOfVars; k++, k_factorial *= k) {
    Integer sum = sum_of_scalar_powers(generic_vector,
				       cone->latticePoints, k);
    mpq_class td = todd(numOfVars, k, ray_scalar_products);
    td /= prod_ray_scalar_products;
    result += sum * td;
  }
  return cone->coefficient * result;
#endif
}

Integer
computeExponentialResidue(listCone *cones, int numOfVars)
{
  listCone *cone;
  Integer result;
  result = 0;
  for (cone = cones; cone != NULL; cone = cone->rest) 
    result += computeExponentialResidue_Single(cone, numOfVars);
  return result;
}
