/*
 * WeightedCountingBuffer.h
 *
 *  Created on: Dec 7, 2014
 *      Author: bedutra
 */

#ifndef WEIGHTEDCOUNTINGBUFFER_H_
#define WEIGHTEDCOUNTINGBUFFER_H_

#include <vector>
#include "latte_gmp.h"
#include "MpqClassLazy.h"
#include "nonlinearOptimization/PolynomialMap.h"

using namespace std;


typedef struct WeightedCountingBuffer
{
	vector<mpq_class_lazy> weights;

	vector<mpq_class_lazy> todds;
	vector<mpq_class_lazy> factorMaster;
	vector<mpq_class_lazy> factor;
	vector<mpq_class_lazy>  result;


} WeightedCountingBuffer;

class WeightedExponentialTable
{
public:
	int linFormPow;
	vector<vector<mpq_class_lazy> > weights;
	vec_ZZ linForm;
	PolynomialMap sPoly;
	WeightedExponentialTable * next;
};


#endif /* WEIGHTEDCOUNTINGBUFFER_H_ */
