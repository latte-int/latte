/*
 * BoxOptimization.h
 *
 *  Created on: Apr 23, 2014
 *      Author: bedutra
 */

#ifndef BOXOPTIMIZATION_H_
#define BOXOPTIMIZATION_H_

#include "latte_gmp.h"
#include "barvinok/dec.h"
#include "barvinok/barvinok.h"
#include "integration/burstTrie.h"
#include "integration/PolyTrie.h"
#include "integration/multiply.h"
#include <map>
#include <ostream>


class PolynomialMap
{
public:
	std::map<int, RationalNTL> terms;
	PolynomialMap & operator+=(const PolynomialMap &rhs);
	bool  operator==(const int rhs);
	friend std::ostream & operator<<(std::ostream& os, const PolynomialMap & rhs);
	RR eval(const RR & s) const;
	void mult(const RationalNTL &rhs);
	void print(ostream &out);
};



class BoxOptimization
{
public:
	RR U;
	RR L;
	ZZ N;
	monomialSum originalPolynomial; //d + 1 variables. (f + s)
	monomialSum currentPolynomial; //(f + s)^currentPower
	int currentPower;

	vec_ZZ lowerBound, upperBound;

	//BurstTrie<PolynomialMap, ZZ> *theTrie; //todo: constructor + if(theTrie) then delete old when calling new
	PolynomialMap currentMap;
	int currentMapPower;

	BoxOptimization();
	void setPolynomial(const vec_ZZ &lowBound, const vec_ZZ &upBound, const monomialSum & poly);
	void setPower(int k);
	void findRange(int itr);
	void findNewUpperbound();
	void findNewLowerbound();

	RR maximumUpperbound();
	RR maximumLowerBound();
};


class WeightedBoxProducer: public ConeProducer {
private:
	vec_ZZ lowerBound, upperBound; //!< assumes the dimension of the box is lowerBound[i] < upperBound[i]
public:
	WeightedBoxProducer(const vec_ZZ & lowerb, const vec_ZZ & upperb);
	virtual ~WeightedBoxProducer();

	/**
	 * Produces the tangent cones of a full dimensional integer-vertex box.
	 */
	virtual void Produce(ConeConsumer &consumer);
};


/**
 * @param lowerBound, upperBound of a box. Note that this bounding box must be full dimensional, meaning lowerBound[i] < upperBound[i]
 * @return the weighted lattice point count
 */
mpq_class computeWeightedCountingBox(const vec_ZZ &lowerBound, const vec_ZZ &upperBound, const linFormSum &originalLinearForm);


/**
 * likewise, but for just one linear form.
 */
mpq_class computeWeightedCountingBox_singleForm(const vec_ZZ &lowerBound, const vec_ZZ &upperBound, const ZZ* linFormExps, const int degree, const RationalNTL & coef);

#endif /* BOXOPTIMIZATION_H_ */
