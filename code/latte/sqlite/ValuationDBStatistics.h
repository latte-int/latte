/*
 * ValuationDBStatistics.h
 *
 *  Created on: Feb 11, 2011
 *      Author: bedutra
 */

#ifndef VALUATIONDBSTATISTICS_H_
#define VALUATIONDBSTATISTICS_H_

struct ValuationDBStatistics
{
	//description parameters
	int dim;
	int vertexCount;
	int degree;
	bool useDual;

	//output parameters.
	double avgTriangulationTime;
	double avgLawrenceTime;

	double minTriangulationTime;
	double minLawrenceTime;

	double maxTriangulationTime;
	double maxLawrenceTime;

	int totalFinishedTriangulationTestCases;
	int totalFinishedLawrenceTestCases;

	int totalTestCases; //includes finished and not finished. same value for Lawrence and triangulation

	bool manuallyLimitedLawrence; //true = this test case has a "-2"
	bool manuallyLimitedTriangulation;
};



#endif /* VALUATIONDBSTATISTICS_H_ */
