/* ExponentialApprox.cpp -- Sample lattice points in fundamental
   parallelepipeds, using the exponential substitution

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
#include "dual.h"
#include "print.h"
#include "ExponentialSubst.h"
#include "ExponentialApprox.h"
#include "latte_random.h"
#include "genFunction/piped.h"
#include "genFunction/IntCombEnum.h"

void Write_Exponential_Sample_Formula_Single_Cone_Parameters::InitializeComputation()
{
  Exponential_Single_Cone_Parameters::InitializeComputation();
  approximate_result = 0;
  total_lower_bound = 0.0;
  total_upper_bound = 0.0;
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
      = scaleRationalVectorToInteger(cone->vertex->vertex, Number_of_Variables,
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
    double this_total_lower_bound = 0.0;
    double this_total_upper_bound = 0.0;
    for (k = 0; k<=Number_of_Variables; k++) {
      ZZ l = cone->coefficient * lower_bounds[k] * abs(cone->determinant);
      double lower_contrib = convert_ZZ_to_mpz(l).get_d() * weights[k].get_d();
      ZZ u = cone->coefficient * upper_bounds[k] * abs(cone->determinant);
      double upper_contrib = convert_ZZ_to_mpz(u).get_d() * weights[k].get_d();

      if (lower_contrib < upper_contrib) {
	total_lower_bounds[k] += lower_contrib;
	total_upper_bounds[k] += upper_contrib;
      }
      else {
	total_lower_bounds[k] += upper_contrib;
	total_upper_bounds[k] += lower_contrib;
      }
      this_total_lower_bound += total_lower_bounds[k];
      this_total_upper_bound += total_upper_bounds[k];
    }
    total_lower_bound += this_total_lower_bound;
    total_upper_bound += this_total_upper_bound;

    /* Output */
    
    //printConeToFile(stream, cone, Number_of_Variables);
    
    stream << "*** Cone with index " << abs(cone->determinant) << endl;
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
    stream << endl;
    stream << "Total lower bound: " << this_total_lower_bound << endl;
    stream << "Total upper bound: " << this_total_upper_bound << endl;

    /* Sample */

    {
      PointsInParallelepipedGenerator generator(cone, Number_of_Variables);
      const vec_ZZ &max_multipliers = generator.GetMaxMultipliers();
      
      // For the beginning, see what kind of approximate we get when
      // we sample as many points as the parallelepiped has.
      // If this is already bad... 
      int actual_num_samples;
      if (abs(cone->determinant) == 1) {
	/* don't bother with unimodular cones */
	actual_num_samples = 1;
      }
      else {
	if (num_samples >= 1) {
	  actual_num_samples = num_samples;
	}
	else {
	  actual_num_samples = (int) (sampling_factor * to_int(abs(cone->determinant)));
	  if (actual_num_samples == 0) actual_num_samples = 1;
#define SAMPLE_LIMIT 1000000
	  if (actual_num_samples > SAMPLE_LIMIT) {
	    cerr << "Refusing to sample more than "
		 << SAMPLE_LIMIT << "points. Sorry." << endl;
	    exit(1);
	  }
	}
      }
      int length = max_multipliers.length();
      vec_ZZ multipliers;
      multipliers.SetLength(length);
      mpq_class sum = 0;
      // Sampling
      int i;
      for (i = 0; i<actual_num_samples; i++) {
	unsigned int j;
	for (j = 0; j<length; j++)
	  multipliers[j] = uniform_random_number(ZZ::zero(), max_multipliers[j] - 1);
	vec_ZZ lattice_point = generator.GeneratePoint(multipliers);
	for (k = 0; k<=dimension; k++) {
	  sum += convert_ZZ_to_mpz(scalar_power(generic_vector,
						lattice_point, k))
	    * weights[k];
	}
      }
      mpq_class scale_factor(convert_ZZ_to_mpz(abs(cone->determinant)),
			     actual_num_samples);
      scale_factor.canonicalize();
      mpq_class contribution = cone->coefficient * scale_factor * sum;
      approximate_result += contribution;
      stream << "Sampled " << actual_num_samples << " points, contribution: "
	     << contribution.get_d() << endl << "  ";
#ifdef DO_EXACT_COMPUTATION_TOO
      // Exact computation, for comparison
      {
	sum = 0;
	int *n = new int[length];
	int i;
	for (i = 0; i<length; i++)
	  n[i] = max_multipliers[i];
	IntCombEnum iter_comb(n, length);
	iter_comb.decrementUpperBound();
	int *next;
	while((next = iter_comb.getNext())) {
	  vec_ZZ lattice_point = generator.GeneratePoint(next);
	  for (k = 0; k<=dimension; k++) {
	    sum += convert_ZZ_to_mpz(scalar_power(generic_vector,
						  lattice_point, k))
	      * weights[k];
	  }
	}
	delete[] n;
	mpq_class contribution = cone->coefficient * sum;
	result += contribution;
	stream << "Exact contribution: " << contribution.get_d() << endl << endl;
      }
#endif
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
  cout << "*** Lower bound: " << param.total_lower_bound << endl;
  cout << "*** Upper bound: " << param.total_upper_bound << endl;
  cout << "*** Estimate obtained by sampling: "
       << param.approximate_result.get_d() << endl;
#ifdef DO_EXACT_COMPUTATION_TOO
  cout << "*** Exact answer: "
       << param.result.get_d() << endl;
#endif
}
