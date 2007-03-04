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
#include "triangulation/triangulate.h"
#include "Irrational.h"

ConeApproximationData::ConeApproximationData(listCone *a_cone,
					     const Write_Exponential_Sample_Formula_Single_Cone_Parameters &params)
  : cone(a_cone),
    ray_scalar_products(INIT_SIZE, params.Number_of_Variables)
{
  int Number_of_Variables = params.Number_of_Variables;

  weights
    = computeExponentialResidueWeights(params.generic_vector,
				       cone,
				       Number_of_Variables);
  int dimension = weights.size() - 1;
  
  /*** Compute bounds for the phi function. ***/

  /* Scalar products of the rays. */
  {
    int j;
    listVector *ray;
    for (j = 0, ray = cone->rays; ray != NULL; j++, ray = ray->rest) {
      ZZ inner;
      InnerProduct(inner, params.generic_vector, ray->first);
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
	       params.generic_vector, scaled_vertex);
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
  mpq_vector total_lower_bounds(Number_of_Variables + 1);
  mpq_vector total_upper_bounds(Number_of_Variables + 1);

  this_total_lower_bound = 0;
  this_total_upper_bound = 0;
  int k;
  for (k = 0; k<=Number_of_Variables; k++) {
    ZZ l = cone->coefficient * lower_bounds[k] * abs(cone->determinant);
    mpq_class lower_contrib = convert_ZZ_to_mpz(l) * weights[k];
    ZZ u = cone->coefficient * upper_bounds[k] * abs(cone->determinant);
    mpq_class upper_contrib = convert_ZZ_to_mpz(u) * weights[k];

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
}

double ConeApproximationData::GetWeight()
{
  mpq_class difference = this_total_upper_bound - this_total_lower_bound;
  return difference.get_d();
}

void Write_Exponential_Sample_Formula_Single_Cone_Parameters::InitializeComputation()
{
  Exponential_Single_Cone_Parameters::InitializeComputation();
  approximate_result = 0;
  total_lower_bound = 0;
  total_upper_bound = 0;
  stream << "*** Computation with new generic vector " << endl;
}

int
Write_Exponential_Sample_Formula_Single_Cone_Parameters::ConsumeCone(listCone *cone)
{
  // We receive triangulated cones.
  assert(cone->rest == NULL);
  
  switch (decomposition) {
  case BarvinokParameters::IrrationalPrimalDecomposition:
    // Do Barvinok decomposition on the primal cones.
    dualizeBackCones(cone, Number_of_Variables);
    irrationalizeCone(cone, Number_of_Variables);
    break;
  case BarvinokParameters::IrrationalAllPrimalDecomposition:
    // FIXME: triangulateCone may feed us simplicial cones for which
    // the determinant has not been computed yet, but the facets are!
    // But computeDetAndFacetsOfSimplicialCone will complain about
    // existing facets...
    freeListVector(cone->facets);
    cone->facets = NULL;
    computeDetAndFacetsOfSimplicialCone(cone, Number_of_Variables);
    break;
  }

  EvaluateCone(cone);
  return 1;
}

void
Write_Exponential_Sample_Formula_Single_Cone_Parameters::EvaluateCone(listCone *cone)
{
  if (abs(cone->determinant) <= max_determinant) {
    /* Small cone, so do exact evaluation. */
    mpq_class result
      = computeExponentialResidue_Single(generic_vector, cone, Number_of_Variables);
    total_lower_bound += result;
    total_upper_bound += result;
    cout << "* [" << abs(cone->determinant) << "], ";
    freeCone(cone);
  }
  else {
    // Compute approximation data and add to heap.
    ConeApproximationData *data = new ConeApproximationData(cone, *this);
    total_lower_bound += data->this_total_lower_bound;
    total_upper_bound += data->this_total_upper_bound;
    heap_insert_dblweight(cone_heap, data, data->GetWeight());
    //cout << "Received cone with bound difference: " << data->GetWeight() << endl;
    cout << data->GetWeight() << " [" << abs(cone->determinant) << "], ";
  }
}

#if 0
{
  try {
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
#endif

void
Write_Exponential_Sample_Formula_Single_Cone_Parameters::ShowStats()
{
  mpq_class diff = total_upper_bound - total_lower_bound;
  cout << "L: " << total_lower_bound.get_d()
       << " U: " << total_upper_bound.get_d()
       << " D: " << diff.get_d()
       << " #: " << heap_num_fill(cone_heap) << ". ";
}

void
decomposeAndWriteExponentialSampleFormula(listCone *cones,
					  Write_Exponential_Sample_Formula_Single_Cone_Parameters &param)
{
  param.InitializeComputation();
  param.cone_heap = heap_alloc(HEAP_DBL_MAX, /*FIXME: fixed size!*/ 100000,
			       NULL);
  cout << "Initial cones: ";
  {
    listCone *cone;
    for (cone = cones; cone!=NULL; cone = cone->rest) {
      triangulateCone(cone, param.Number_of_Variables, &param, param);
      // calls ConsumeCone, filling the heap
    }
  }
  cout << endl;
  while (heap_nonempty_p(param.cone_heap)) {
    param.ShowStats();
    listCone *cone;
    {
      ConeApproximationData *data = (ConeApproximationData*) heap_top(param.cone_heap);
      heap_pop(param.cone_heap);
      param.total_lower_bound -= data->this_total_lower_bound;
      param.total_upper_bound -= data->this_total_upper_bound;
      
      cout << "Cone bound diff " << data->GetWeight()
	   << " [Det " << abs(data->cone->determinant) << "]"
	   << " --> " << flush;
      cone = data->cone;
      delete data;
    }
    vector <listCone *> Cones(param.Number_of_Variables);
    vec_ZZ Dets;
    Dets.SetLength(param.Number_of_Variables);
    bool barvinok_success
      = barvinokStep(cone, Cones, Dets, param.Number_of_Variables, &param);
    freeCone(cone);
    int i;
    for (i = 0; i<param.Number_of_Variables; i++) {
      if (Cones[i] != NULL)
	param.EvaluateCone(Cones[i]);
    }
    cout << endl;
  }
  cout << "*** Lower bound: " << param.total_lower_bound.get_d() << endl;
  cout << "*** Upper bound: " << param.total_upper_bound.get_d() << endl;
}

#if 0
  cout << "*** Estimate obtained by sampling: "
       << param.approximate_result.get_d() << endl;
#ifdef DO_EXACT_COMPUTATION_TOO
  cout << "*** Exact answer: "
       << param.result.get_d() << endl;
#endif
#endif
