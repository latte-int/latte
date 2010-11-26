/*
 * VolumeAndIntegrationTests.h
 *
 *  Created on: Nov 19, 2010
 *      Author: Brandon Dutra and Gregory Pinto
 */

#ifndef VOLUMEANDINTEGRATIONTESTS_H_
#define VOLUMEANDINTEGRATIONTESTS_H_

//these headers are needed by the findAllLatteFilesInDirectory function.
#include <dirent.h>
#include <errno.h>
#include <string.h>


#include "buildPolytopes/BuildRandomPolytope.h"
#include "buildPolytopes/BuildHypersimplexEdgePolytope.h"
#include "rational.h"
#include "../valuation.h"


//These functions are used to populate data tables for an integration paper we are working on.
namespace IntegrationPaper
{
	//makes a random polynomial in dim variables of a set degree.
	string makeMonomial(const int dim, int totalDegree);

	//makes many random monomials
	string makePolynomial(const int dim, const int total, const int numMonomials);

	//Giving a path to a directory, will return a list of all the .latte files in that directory. This is NOT recursive and only works on unix boxes I think.
	void findAllLatteFilesInDirectory(const string &dir, vector<string> &latteFiles);


	//integrates each file with a polynomial/monomial and save it to a log file.
	void integrateFiles(const string &logFileName, const vector<string> &files, const int dim, const int polynomialDegree);

}//IntegrationPaper

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
