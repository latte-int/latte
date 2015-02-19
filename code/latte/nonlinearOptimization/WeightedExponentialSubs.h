/*
 * WeightedExponentialSubs.h
 *
 *  Created on: Apr 24, 2014
 *      Author: bedutra
 */

#ifndef WEIGHTEDEXPONENTIALSUBS_H_
#define WEIGHTEDEXPONENTIALSUBS_H_

#include "latte_gmp.h"
#include "barvinok/dec.h"
#include "barvinok/barvinok.h"
#include "integration/burstTrie.h"
#include "integration/PolyTrie.h"
#include "nonlinearOptimization/WeightedCountingBuffer.h"

/**
 * Cone consumer.
 * For finding \sum_{x \in cone} <x, linForm>^linFormPow
 *
 */
class Weighted_Exponential_Single_Cone_Parameters
  : public Generic_Vector_Single_Cone_Parameters {
public:
  mpq_class result;				// result of \sum_{x \in box} <x, linForm>^linFormPow
  vec_ZZ linForm;				// has length n (ambient dimension), but real dimension is Number_of_Variables
  int linFormPow;
  WeightedCountingBuffer *wcb;
  Weighted_Exponential_Single_Cone_Parameters() :
    result(0) {};
  Weighted_Exponential_Single_Cone_Parameters(const BarvinokParameters &params, const vec_ZZ & linform) :
    Generic_Vector_Single_Cone_Parameters(params),
    result(0),
    linForm(linform) {};
  virtual void InitializeComputation();
  virtual int ConsumeCone(listCone *cone);
};


class Weighted_Exponential_Single_Cone_Parameters_BranchBound
  : public Generic_Vector_Single_Cone_Parameters {
public:
  int index;
  WeightedExponentialTable *table;
  vec_ZZ linForm;				// has length n (ambient dimension), but real dimension is Number_of_Variables
  int linFormPow;
  WeightedCountingBuffer *wcb;
  Weighted_Exponential_Single_Cone_Parameters_BranchBound() {};
  Weighted_Exponential_Single_Cone_Parameters_BranchBound(const BarvinokParameters &params, const vec_ZZ & linform) :
    Generic_Vector_Single_Cone_Parameters(params),
    linForm(linform) {};
  virtual void InitializeComputation();
  virtual int ConsumeCone(listCone *cone);
};

/**
 * Computes the weighted lattice point count where the weight function is sum of powers of linear forms.
 *
 * @param cone: list of simple tangent cones, they might not be unimodular.
 * @param params: Number_of_Variables needs to be set
 * @param originalLinearForm: sum of linear forms
 * @throw LattException::bug_NotImplementedHere
 * @return the sum of (l_1x_1 + \cdots + l_dx_d)^M for the integer points of a Polytope.
 *  That is, we find
 *
 *  \sum_{l} \sum_{x \in P \cap \Z^d} (l_1x_1 + \cdots + l_dx_d)^M as the coefficient of t^M in the series expansion of
 *
 *  M! *  \sum_{tangent cones} \epsilon_i *    \sum_{v \in \Pi} exp(t * <v, l>)
 *                                             --------------------------------
 *                                             \prod_{rays} (1 - exp(t* <r_i, l>))
 *
 */
mpq_class computeWeightedExponentialResidue(WeightedCountingBuffer &wcb, listCone *cone, BarvinokParameters * params, linFormSum &originalLinearForm);

/**
 * Same as computeWeightedExponentialResidue but only processes one linear form
 * @param generic_vector: should be used if the linear form is orthogonal with one of the rays...but this is not implemented yet. So this parameter is meaningless.
 * @throw NotGenericException
 */
mpq_class computeWeightedExponentialResidue_singleForm(WeightedCountingBuffer &wcb, listCone *cone, BarvinokParameters * params, const vec_ZZ &linFormCoeffs, const vec_ZZ &generic_vector, int M);


#endif /* WEIGHTEDEXPONENTIALSUBS_H_ */
