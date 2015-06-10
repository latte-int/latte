/*
 * BoxOptimizationContinuous.h
 *
 *  Created on: Mar 31, 2015
 *      Author: bedutra
 */

#ifndef BOXOPTIMIZATIONCONTINUOUS_H_
#define BOXOPTIMIZATIONCONTINUOUS_H_


#include "latte_gmp.h"
#include "barvinok/dec.h"
#include "barvinok/barvinok.h"
#include "integration/burstTrie.h"
#include "integration/PolyTrie.h"
#include "integration/multiply.h"
#include "nonlinearOptimization/WeightedCountingBuffer.h"
#include "nonlinearOptimization/PolynomialMap.h"


#include <map>

typedef RR RRorQQ;
typedef vec_RR vec_RRorQQ;


//typedef RationalNTL RRorQQ;
//typedef vec_RationalNTL vec_RRorQQ;

class BoxOptimizationContinuous
{
private:

public:
	RR U;
	RR L;
	RRorQQ V; //volume of the box in its affine dimension.
	int zero = 0;

	RRorQQ lipschitz;
	RRorQQ M;

	typedef enum AlgoType
	{
	    naturalSummation
	} AlgoType;


	monomialSum originalPolynomial; //d + 1 variables. (f + s)
	monomialSum currentPolynomial; //(f + s)^currentPower
	int currentPower;

	vec_RRorQQ lowerBound, upperBound;

	//BurstTrie<PolynomialMap, ZZ> *theTrie; //todo: constructor + if(theTrie) then delete old when calling new
	//PolynomialMap currentMap;
	vec_RRorQQ currentMap;
	int currentMapPower;
	double fmaxLowerBound;


	BoxOptimizationContinuous();
	~BoxOptimizationContinuous();
	void setPolynomial(const monomialSum & poly);
	void setBounds(const vec_RRorQQ &lowBound, const vec_RRorQQ &upBound);
	void setPower(int k);
	//void decomposePoly(AlgoType algoType);
	void findSPolynomial(AlgoType algoType, const vec_RRorQQ &lowerBound, const vec_RRorQQ & upperBound);
	int findRange(int nitr);
	void printSpolynomial() const;
	void printStats() const;
	void setFmaxLowerbound(const double & fmaxLB);


	RR maximumLowerBound();
	RR maximumUpperbound();

	//RR sampleLowerBound(monomialSum &poly, const vec_RRorQQ & point);

	BoxOptimizationContinuous & operator=(const BoxOptimizationContinuous & rhs);
	RR evalSpoly(const RR & s) const;
};

#endif /* BOXOPTIMIZATIONCONTINUOUS_H_ */
