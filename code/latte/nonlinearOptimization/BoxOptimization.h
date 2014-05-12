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


#endif /* BOXOPTIMIZATION_H_ */
