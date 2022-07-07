/* ExponentialSubst.cpp -- the exponential substitution

   Copyright 2006 Matthias Koeppe

   This file is part of LattE.
   
   LattE is free software; you can redistribute it and/or modify it
   under the terms of the version 2 of the GNU General Public License
   as published by the Free Software Foundation.

   LattE is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with LattE; if not, write to the Free Software Foundation,
   Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
*/

#include <cassert>
#include "ExponentialSubst.h"
#include "todd/gmp_pow.h"
#include "latte_gmp.h"
#include "todd/todd-expansion.h"
#include "genFunction/piped.h"
#include "dual.h"
#include "cone.h"
#include "print.h"
#include "latte_ntl_integer.h"

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
{
  mpz_class prod_ray_scalar_products;
  return computeExponentialResidueWeights(generic_vector,
					  prod_ray_scalar_products,
					  cone,
					  numOfVars);
}

vec_ZZ
compute_sums_of_scalar_powers(listCone *cone,
			      int numOfVars,
			      const vec_ZZ &generic_vector,
			      BarvinokParameters *params)
{
  computeLatticePointsScalarProducts(cone, numOfVars, generic_vector, params);
  vec_ZZ sum;
  int dimension = numOfVars;
  sum.SetLength(dimension + 1);
  int i;
  int num_points = cone->lattice_points_scalar_products.length();
  for (i = 0; i<num_points; i++) {
    Integer inner = cone->lattice_points_scalar_products[i];
    Integer scalar_power;
    scalar_power = 1;
    int k;
    for (k = 0; k<=dimension; k++) {
      sum[k] += scalar_power;
      scalar_power *= inner;
    }
  }
  return sum;
}

mpz_vector
compute_sums_of_scalar_powers_mpz(listCone *cone,
				  int numOfVars,
				  const vec_ZZ &generic_vector,
				  BarvinokParameters *params)
{
  vec_ZZ sums_of_scalar_powers_zz
    = compute_sums_of_scalar_powers(cone, numOfVars, generic_vector, params);
  mpz_vector sums_of_scalar_powers(numOfVars + 1);
  {
    int i;
    for (i = 0; i<=numOfVars; i++)
      sums_of_scalar_powers[i]
	= convert_ZZ_to_mpz(sums_of_scalar_powers_zz[i]);
  }
  return sums_of_scalar_powers;
}

mpq_class
computeExponentialResidue_Single(const vec_ZZ &generic_vector,
				 listCone *cone, int numOfVars,
				 BarvinokParameters *params)
{
  mpq_vector weights
    = computeExponentialResidueWeights(generic_vector, cone, numOfVars);
  int dimension = weights.size() - 1;
  int k;
  mpq_class result = 0;
#if 1
  /* Equivalent, but faster code: */
  computeLatticePointsScalarProducts(cone, numOfVars, generic_vector, params);
  mpz_vector sum = compute_sums_of_scalar_powers_mpz(cone, numOfVars, generic_vector, params);
  for (k = 0; k<=dimension; k++)
    result += sum[k] * weights[k];
#else
  computePointsInParallelepiped(cone, numOfVars);
  for (k = 0; k<=dimension; k++) {
    Integer sum = sum_of_scalar_powers(generic_vector,
				       cone->latticePoints, k);
    result += convert_ZZ_to_mpz(sum) * weights[k];
  }
#endif
//   cerr << "Cone contributes: "
//        << cone->coefficient << " * " << result << endl;
  return cone->coefficient * result;
}

Integer
computeExponentialResidue(listCone *cones, int numOfVars, BarvinokParameters *params)
{
  listCone *cone;
  do {
    vec_ZZ generic_vector = guess_generic_vector(numOfVars);
    mpq_class result;
    result = 0;
    try {
      for (cone = cones; cone != NULL; cone = cone->rest) 
	result += computeExponentialResidue_Single(generic_vector, cone, numOfVars, params);
      //cerr << "Result: " << result << endl;
      assert(result.get_den()==1);
      return convert_mpz_to_ZZ(result.get_num());
    }
    catch (NotGenericException) {};
    cerr << "New generic vector..." << endl;
  } while (1);
}

int Exponential_Single_Cone_Parameters::ConsumeCone(listCone *cone)
{
  assert(cone->rest == NULL);
  int status = 1;
  try {
    result += computeExponentialResidue_Single(generic_vector,
					       cone, Number_of_Variables, this);
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
#if 0
  gmp_printf("Fun fact: Number of lattice points %Qd = %g\n",
	     param.result.get_mpq_t(), param.result.get_d());
#endif
  assert(param.result.get_den()==1);
  return convert_mpz_to_ZZ(param.result.get_num());
}
