/// Experimental code.

#include <cassert>
#include "dual.h"
#include "print.h"
#include "ExponentialSubst.h"
#include "ExponentialApprox.h"
#include "latte_random.h"
#include "genFunction/piped.h"

void Write_Exponential_Sample_Formula_Single_Cone_Parameters::InitializeComputation()
{
  Exponential_Single_Cone_Parameters::InitializeComputation();
  stream << "*** Computation with new generic vector " << endl;
}

int
Write_Exponential_Sample_Formula_Single_Cone_Parameters::ConsumeCone(listCone *cone)
{
  assert(cone->rest == NULL);
  try {
    mpq_vector weights
      = computeExponentialResidueWeights(generic_vector,
					 cone,
					 Number_of_Variables);
    int dimension = weights.size() - 1;
  
    /*** Compute bounds for the phi function. ***/

    /* Scalar products of the rays. */
    vec_ZZ ray_scalar_products(INIT_SIZE, Number_of_Variables);
    {
      int j;
      listVector *ray;
      for (j = 0, ray = cone->rays; ray != NULL; j++, ray = ray->rest) {
	ZZ inner;
	InnerProduct(inner, generic_vector, ray->first);
	ray_scalar_products[j] = inner;
      }
      assert(j == Number_of_Variables);
    }

    /* Scalar product of the vertex. */
    ZZ vertex_divisor;
    vec_ZZ scaled_vertex
      = scaleRationalVectorToInteger(cone->vertex, Number_of_Variables,
				     vertex_divisor);
    ZZ scaled_vertex_scalar_product;
    ZZ vertex_scalar_product_lower, vertex_scalar_product_upper;
    InnerProduct(scaled_vertex_scalar_product,
		 generic_vector, scaled_vertex);
    ZZ remainder;
    DivRem(vertex_scalar_product_lower, remainder,
	   scaled_vertex_scalar_product, vertex_divisor);
    if (IsZero(remainder))
      vertex_scalar_product_upper = vertex_scalar_product_lower;
    else
      vertex_scalar_product_upper = vertex_scalar_product_lower + 1;

    /* Bounds for <c,x> on v+Pi */
    ZZ lower_bound, upper_bound;
    lower_bound = vertex_scalar_product_lower;
    upper_bound = vertex_scalar_product_upper;
    {
      int j;
      for (j = 0; j<Number_of_Variables; j++) {
	switch (sign(ray_scalar_products[j])) {
	case -1:
	  lower_bound += ray_scalar_products[j];
	  break;
	case +1:
	  upper_bound += ray_scalar_products[j];
	  break;
	}
      }
    }

    /* Bounds for <c,x>^k on v+Pi */
    vec_ZZ lower_bounds(INIT_SIZE, Number_of_Variables + 1);
    vec_ZZ upper_bounds(INIT_SIZE, Number_of_Variables + 1);
    {
      int k;
      ZZ lower_k, upper_k;
      lower_k = 1;
      upper_k = 1;
      for (k = 0; k<=Number_of_Variables; /*EMPTY*/) {
	if (k % 2 == 0) {
	  /* even case */
	  if (sign(upper_bound) == -1) {
	    lower_bounds[k] = upper_k;
	    upper_bounds[k] = lower_k;
	  }
	  else if (sign(lower_bound) == +1) {
	    lower_bounds[k] = lower_k;
	    upper_bounds[k] = upper_k;
	  }
	  else {
	    lower_bounds[k] = 0;
	    upper_bounds[k] = (lower_k > upper_k) ? lower_k : upper_k;
	  }
	}
	else {
	  /* odd case */
	  lower_bounds[k] = lower_k;
	  upper_bounds[k] = upper_k;
	}
	k++;
	lower_k *= lower_bound;
	upper_k *= upper_bound;
      }
    }

    /* Total bounds */
    vector<double> total_lower_bounds(Number_of_Variables + 1);
    vector<double> total_upper_bounds(Number_of_Variables + 1);

    int k;
    for (k = 0; k<=Number_of_Variables; k++) {
      ZZ l = lower_bounds[k] * abs(cone->determinant);
      double lower_contrib = convert_ZZ_to_mpz(l).get_d() * weights[k].get_d();
      ZZ u = upper_bounds[k] * abs(cone->determinant);
      double upper_contrib = convert_ZZ_to_mpz(u).get_d() * weights[k].get_d();

      if (l < u) {
	total_lower_bounds[k] += lower_contrib;
	total_upper_bounds[k] += upper_contrib;
      }
      else {
	total_lower_bounds[k] += upper_contrib;
	total_upper_bounds[k] += lower_contrib;
      }
      total_lower_bound += total_lower_bounds[k];
      total_upper_bound += total_upper_bounds[k];
    }

    /* Output */
    
    //printConeToFile(stream, cone, Number_of_Variables);
    
    stream << "*** Cone with index " << cone->determinant << endl;
    stream << "Approximate Weights: " << endl << "  ";
    mpq_vector::const_iterator i;
    for (i = weights.begin(); i!=weights.end(); ++i) {
      double weight = (*i).get_d();
      stream << weight << " ";
    }
    stream << endl;
    stream << "Lower bounds of k! phi_k: " << endl << "  ";
    for (k = 0; k<=Number_of_Variables; k++)
      stream << convert_ZZ_to_mpz(lower_bounds[k]).get_d() << " ";
    stream << endl;
    stream << "Upper bounds of k! phi_k: " << endl << "  ";
    for (k = 0; k<=Number_of_Variables; k++)
      stream << convert_ZZ_to_mpz(upper_bounds[k]).get_d() << " ";
    stream << endl;
    stream << "Lower bounds of contributions: " << endl << "  ";
    for (k = 0; k<=Number_of_Variables; k++)
      stream << total_lower_bounds[k] << " ";
    stream << endl;
    stream << "Upper bounds of contributions: " << endl << "  ";
    for (k = 0; k<=Number_of_Variables; k++)
      stream << total_upper_bounds[k] << " ";
    stream << endl << endl;

    /* Sample */

    {
      int i;
#define SAMPLE_LIMIT 1000000
      if (abs(cone->determinant) > SAMPLE_LIMIT) {
	cerr << "Refusing to handle a cone with index more than "
	     << SAMPLE_LIMIT << ". Sorry." << endl;
	exit(1);
      }
      PointsInParallelepipedGenerator generator(cone, Number_of_Variables);
      vector<int> max_multipliers = generator.GetMaxMultipliers();
      
      // For the beginning, see what kind of approximate we get when
      // we sample as many points as the parallelepiped has.
      // If this is already bad... 
      int num_samples = to_int(cone->determinant);
      int length = max_multipliers.size();
      int *multipliers = new int[length];
      mpq_class sum = 0;
      for (i = 0; i<num_samples; i++) {
	unsigned int j;
	for (j = 0; j<max_multipliers.size(); j++)
	  multipliers[j] = uniform_random_number(0, max_multipliers[j] - 1);
	vec_ZZ lattice_point = generator.GeneratePoint(multipliers);
	for (k = 0; k<=dimension; k++) {
	  sum += convert_ZZ_to_mpz(scalar_power(generic_vector,
						lattice_point, k))
	    * weights[k];
	}
      }
      result += cone->coefficient * mpq_class(convert_ZZ_to_mpz(abs(cone->determinant)),
					      num_samples)
	* sum;
      delete[] multipliers;
    }
    return 1;
  }
  catch (NotGenericException) {
    return -1;
  }
}

void
decomposeAndWriteExponentialSampleFormula(listCone *cones,
					  Write_Exponential_Sample_Formula_Single_Cone_Parameters &param)
{
  barvinokDecomposition_List(cones, param);
  cout << endl << "*** Lower bound: " << param.total_lower_bound << endl;
  cout << endl << "*** Upper bound: " << param.total_upper_bound << endl;
  cout << endl << "*** Estimate obtained by sampling: "
       << param.result.get_d() << endl;
}
