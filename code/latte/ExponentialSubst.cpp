#include <cassert>
#include "ExponentialSubst.h"
#include "todd/gmp_pow.h"
#include "latte_gmp.h"
#include "todd/todd-expansion.h"
#include "genFunction/piped.h"
#include "dual.h"
#include "cone.h"
#include "barvinok/ConeDecom.h"
#include "print.h"

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

mpq_vector
computeExponentialResidueWeights(const vec_ZZ &generic_vector,
				 mpz_class &prod_ray_scalar_products,
				 const listCone *cone, int numOfVars)
  throw(NotGenericException)
{
  // Compute dimension; can be smaller than numOfVars
  int dimension = 0;
  listVector *ray;
  for (ray = cone->rays; ray != NULL; ray = ray->rest)
    dimension++;
  vector<mpz_class> ray_scalar_products(dimension);
  //mpz_class prod_ray_scalar_products;
  prod_ray_scalar_products = 1;
  {
    int k;
    for (k = 0, ray = cone->rays; ray != NULL; k++, ray = ray->rest) {
      ZZ inner;
      InnerProduct(inner, generic_vector, ray->first);
      ray_scalar_products[k] = convert_ZZ_to_mpz(inner);
      if (ray_scalar_products[k] == 0) {
	static NotGenericException not_generic;
	throw not_generic;
      }
      prod_ray_scalar_products *= ray_scalar_products[k];
    }
  }
  int k;
  mpz_class k_factorial;
  mpq_vector weights(dimension + 1);
  mpq_vector todds = evaluate_todd(ray_scalar_products);
  for (k = 0, k_factorial = 1; k<=dimension; k++, k_factorial *= k) {
    mpq_class td = todds[dimension - k];
    td /= prod_ray_scalar_products;
    weights[k] = td / k_factorial;
  }
  return weights;
}

mpq_vector
computeExponentialResidueWeights(const vec_ZZ &generic_vector,
				 const listCone *cone, int numOfVars)
  throw(NotGenericException)
{
  mpz_class prod_ray_scalar_products;
  return computeExponentialResidueWeights(generic_vector,
					  prod_ray_scalar_products,
					  cone,
					  numOfVars);
}

mpq_class
computeExponentialResidue_Single(const vec_ZZ &generic_vector,
				 const listCone *cone, int numOfVars)
{
  mpq_vector weights
    = computeExponentialResidueWeights(generic_vector, cone, numOfVars);
  int dimension = weights.size() - 1;
  int k;
  mpq_class result = 0;
  for (k = 0; k<=dimension; k++) {
    Integer sum = sum_of_scalar_powers(generic_vector,
				       cone->latticePoints, k);
    result += convert_ZZ_to_mpz(sum) * weights[k];
  }
//   cout << "Cone contributes: "
//        << cone->coefficient << " * " << result << endl;
  return cone->coefficient * result;
}

Integer
computeExponentialResidue(const listCone *cones, int numOfVars)
{
  const listCone *cone;
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

int Exponential_Single_Cone_Parameters::ConsumeCone(listCone *cone)
{
  assert(cone->rest == NULL);
  computePointsInParallelepiped(cone, Number_of_Variables);
  int status = 1;
  try {
    result += computeExponentialResidue_Single(generic_vector,
					       cone, Number_of_Variables);
  } catch (NotGenericException) {
    status = -1;
  }
  if (Total_Uni_Cones % 1000 == 0) {
    gmp_printf("Fun fact: Number of lattice points currently %g\n",
	       result.get_d());
  }
  freeListCone(cone);
  return status;
}

void Exponential_Single_Cone_Parameters::InitializeComputation()
{
  Generic_Vector_Single_Cone_Parameters::InitializeComputation();
  result = 0;
}

Integer
decomposeAndComputeExponentialResidue(listCone *cones, 
				      Exponential_Single_Cone_Parameters &param)
{
  //printListCone(cones, param.Number_of_Variables);
  barvinokDecomposition_List(cones, param);
  assert(param.result.get_den()==1);
  return convert_mpz_to_ZZ(param.result.get_num());
}
