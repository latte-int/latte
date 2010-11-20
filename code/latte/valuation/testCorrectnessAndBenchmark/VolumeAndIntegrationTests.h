/*
 * VolumeAndIntegrationTests.h
 *
 *  Created on: Nov 19, 2010
 *      Author: bedutra
 */

#ifndef VOLUMEANDINTEGRATIONTESTS_H_
#define VOLUMEANDINTEGRATIONTESTS_H_

#include "buildPolytopes/BuildRandomPolytope.h"
#include "buildPolytopes/BuildHypersimplexEdgePolytope.h"
#include "rational.h"
#include "../valuation.h"

namespace VolumeTests
{

void printVolumeTest(const RationalNTL &correctVolumeAnswer,
		const Valuation::ValuationContainer & valuationResults, const string &file,
		const string &comments);

void runOneTest(int ambientDim, int numPoints);
void runTests();

void runHyperSimplexTests();
void runBirkhoffTests();

}//namespace VolmeTests



#endif /* VOLUMEANDINTEGRATIONTESTS_H_ */
