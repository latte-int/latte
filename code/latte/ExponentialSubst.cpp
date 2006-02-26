#include <cassert>
#include "ExponentialSubst.h"
#include "todd/todd.h"
#include "todd/gmp_pow.h"
#include "latte_gmp.h"
#include "todd/todd-expansion.h"

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

#define MODULUS 1000000000
static vec_ZZ
guess_generic_vector(int numOfVars)
{
  vec_ZZ result;
  result.SetLength(numOfVars);
  int i;
  for (i = 0;  i < numOfVars; i++)
    result[i] = (rand () % MODULUS) * ((rand() % 2) * 2 - 1);
  return result;
}

struct NotGenericException {};
NotGenericException not_generic;

mpq_class
computeExponentialResidue_Single(const vec_ZZ &generic_vector,
				 listCone *cone, int numOfVars)
{
  // Compute dimension; can be smaller than numOfVars
  int dimension = 0;
  listVector *ray;
  for (ray = cone->rays; ray != NULL; ray = ray->rest)
    dimension++;
  vector<mpz_class> ray_scalar_products(dimension);
  mpz_class prod_ray_scalar_products;
  prod_ray_scalar_products = 1;
  cout << "Scalar products: ";
  {
    int k;
    for (k = 0, ray = cone->rays; ray != NULL; k++, ray = ray->rest) {
      ZZ inner;
      InnerProduct(inner, generic_vector, ray->first);
      ray_scalar_products[k] = convert_ZZ_to_mpz(inner);
      cout << ray_scalar_products[k] << " ";
      if (ray_scalar_products[k] == 0)
	throw not_generic;
      prod_ray_scalar_products *= ray_scalar_products[k];
    }
  }
  cout << endl;
  int k;
  mpz_class k_factorial;
  mpq_class result;
  mpq_vector todds = evaluate_todd(ray_scalar_products);
  result = 0;
  for (k = 0, k_factorial = 1; k<=dimension; k++, k_factorial *= k) {
    Integer sum = sum_of_scalar_powers(generic_vector,
				       cone->latticePoints, k);
    mpq_class td = todds[dimension - k];
    cerr << "scalar product still: ";
    for (int l = 0; l<dimension; l++)
      cerr << ray_scalar_products[l];
    cerr << endl;
    mpq_class td_maple = todd(dimension, dimension - k, ray_scalar_products);
    if (td != td_maple) {
      cerr << "Todd discrepancy: " << td << " (us) vs. "
	   << td_maple << " (Maple)" << endl;
    }
    td /= prod_ray_scalar_products;
    result += convert_ZZ_to_mpz(sum) * td / k_factorial;
  }
  cout << "Cone contributes: "
       << cone->coefficient << " * " << result << endl;
  return cone->coefficient * result;
}

Integer
computeExponentialResidue(listCone *cones, int numOfVars)
{
  listCone *cone;
  do {
    vec_ZZ generic_vector = guess_generic_vector(numOfVars);
    mpq_class result;
    result = 0;
    try {
      for (cone = cones; cone != NULL; cone = cone->rest) 
	result += computeExponentialResidue_Single(generic_vector, cone, numOfVars);
      //cout << "Result: " << result << endl;
      assert(result.get_den()==1);
      return convert_mpz_to_ZZ(result.get_num());
    }
    catch (NotGenericException) {};
  } while (1);
}