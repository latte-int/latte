/*
 * VolumeAndIntegrationTests.cpp
 *
 *  Created on: Nov 19, 2010
 *      Author: Brandon Dutra and Gregory Pinto
 */

#include "VolumeAndIntegrationTests.h"






/**
 * Checks to see if the triangulation and lawrence volume equal the expected volume.
 */
void VolumeTests::printVolumeTest(const RationalNTL &correctVolumeAnswer,
		const Valuation::ValuationContainer & valuationResults, const string &file,
		const string &comments)
{
	RationalNTL lawrence, triangulate;
	int found = 0;

	for(int i = 0; i < valuationResults.answers.size(); ++i)
		if ( valuationResults.answers[i].valuationType == Valuation::ValuationData::volumeLawrence)
		{
			++found;
			lawrence = valuationResults.answers[i].answer;
		}
		else if ( valuationResults.answers[i].valuationType == Valuation::ValuationData::volumeTriangulation)
		{
			++found;
			triangulate = valuationResults.answers[i].answer;
		}
	if ( found != 2)
	{
		cerr << "VolumeTests::printVolumeTest: could not find lawrence and triangulation volumes." << endl;
		exit(1);
	}

	if (correctVolumeAnswer != lawrence || correctVolumeAnswer
			!= triangulate)
	{
		cerr << "******* ERROR ******" << endl;
		cerr << "correct answer: " << correctVolumeAnswer << endl;
		cerr << "lawrence: " << lawrence << endl;
		cerr << "triangulate: " << triangulate << endl;
		cerr << "see file " << file.c_str() << endl;
		exit(1); //dont' delete the latte file.
	}//if error
	else
		cerr << comments.c_str() << " CORRECT!" << endl;
}//printVolumeTest

/**
 * Calls polymake to make a random interger (or rational) vertex polytope, and then makes the latte file.
 * The latte file is then passed into mainValuationDriver() to find the volume
 *
 * We cannot check our volume with polymake for low-dimensional polytopes.
 */
void VolumeTests::runOneTest(int ambientDim, int numPoints)
{
	const char * argv[] =
	{ "runTests()", "--valuation=volume", "--all", 0 };
	stringstream comments;
	comments << "Making random integer polytope with " << numPoints
			<< " points in R^" << ambientDim << " for volume testing";

	BuildRandomPolytope buildPolytope;
	buildPolytope.setIntegerPoints(false); //make random rational points.
	buildPolytope.makePoints(ambientDim, numPoints);
	buildPolytope.buildLatteHRepFile();

	string file = buildPolytope.getLatteHRepFile();

	char * sFile = new char[file.size() + 1];
	strcpy(sFile, file.c_str());
	argv[3] = sFile;
	Valuation::mainValuationDriver(argv, 4);
	delete[] sFile;
	buildPolytope.deleteLatteHRepFile();
	buildPolytope.deletePolymakeFile();
}//RunOneTest

/**
 * Runs many random tests by calling runOneTest
 */
void VolumeTests::runTests()
{
	int startAmbientDim = 6, endAmbientDim = 50;
	int pointStepSize = 5;

	for (int ambientDim = startAmbientDim; ambientDim < endAmbientDim; ambientDim
			= ambientDim + 3)
	{
		for (int numberPoints = startAmbientDim / 2; numberPoints
				< startAmbientDim / 4 + startAmbientDim; numberPoints
				= numberPoints + pointStepSize)
			runOneTest(ambientDim, numberPoints);

	}//for ambientDim

}//runTests

/**
 * Finds the volume of hypersimplex polytopes and checks for correctness.
 */
void VolumeTests::runHyperSimplexTests()
{
	const char * argv[] =
	{ "runHyperSimplexTests()", "--valuation=volume", "--all", 0 };
	//   n  k  num/denom
	int hyperSimplexData[][4] =
	{ {4, 1, 1, 6},
	 {4, 2, 2, 3},
	 {5, 1, 1, 24},
	 {5, 2, 11, 24},
	 {6, 1, 1, 120},
	 {6, 2, 13, 60},
	 {6, 3, 11, 20},
	 {7, 1, 1, 720},
	 {7, 2, 19, 240},
	 {7, 3, 151, 360},
	 /*{8, 1, 1, 5040},
	 {8, 2, 1, 42},
	 {8, 3, 397, 1680},
	 {8, 4, 151, 315},
	 {9, 1, 1, 40320},
	 {9, 2, 247, 40320},
	 {9, 3, 477, 4480},
	 {9, 4, 15619, 40320},
	 {10, 1, 1, 362880},
	 {10, 2, 251, 181440},
	 {10, 3, 913, 22680},
	 {10, 4, 44117, 181440},
	{ 10, 5, 15619, 36288 },

			{ 11, 1, 1, 3628800 },
			{ 11, 2, 1013, 3628800 },
			{ 11, 3, 299, 22680 },
			{ 11, 4, 56899, 453600 },
			{ 11, 5, 655177, 1814400 }, //these have never been tested...they took to long.
			{ 12, 1, 1, 39916800 },
			{ 12, 2, 509, 9979200 },
			{ 12, 3, 50879, 13305600 },
			{ 12, 4, 1093, 19800 },
			{ 12, 5, 1623019, 6652800 },
			{ 12, 6, 655177, 1663200 }
	*/};//hyperSimplexData

	int numberTestCases = 34;

	Valuation::ValuationContainer volumeAnswer;
	for (int i = 0; i < numberTestCases; ++i)
	{
		stringstream comments;
		BuildHypersimplexEdgePolytope hyperSimplex;
		hyperSimplex.generatePoints(hyperSimplexData[i][0],
				hyperSimplexData[i][1]);

		comments << "finding volume of Hypersimplex(" << hyperSimplexData[i][0]
				<< ", " << hyperSimplexData[i][1] << ")";
		hyperSimplex.buildPolymakeFile();
		hyperSimplex.buildLatteHRepFile();



		string file = hyperSimplex.getLatteHRepFile();

		char * sFile = new char[file.size() + 1];
		strcpy(sFile, file.c_str());
		argv[3] = sFile;
		volumeAnswer = Valuation::mainValuationDriver(argv, 4);
		delete[] sFile;
		hyperSimplex.deletePolymakeFile();
		hyperSimplex.deleteLatteHRepFile();
		RationalNTL correctVolumeAnswer(hyperSimplexData[i][2],
				hyperSimplexData[i][3]);
		VolumeTests::printVolumeTest(correctVolumeAnswer, volumeAnswer, file, comments.str());
	}//for i.
}//runHyperSimplexTests


/**
 * Finds the volume of Birkhoff polytopes and checks for correctness.
 */
void VolumeTests::runBirkhoffTests()
{

	string birkhoff[] =
	{ "../../../../EXAMPLES/birkhoff/birkhoff-5.latte",
			"../../../../EXAMPLES/birkhoff/birkhoff-6.latte",
			"birkhoff7.latte.vrep" };
	string birkhoffVolume[][2] =
	{
	{ "188723", "836911595520" }, //5

			{ "9700106723", "10258736801144832000000" }, //6

			{ "225762910421308831", "4709491654300668677115504230400000000" } //7
	};
	int numberTestCases = 3;

	Valuation::ValuationContainer volumeAnswer;
	const char * argv[] =
	{ "runBirkhoffTests()", "--valuation=volume", "--all", 0 };

	for (int i = 0; i < numberTestCases; ++i)
	{
		char * sFile = new char[birkhoff[i].length() + 1];
		strcpy(sFile, birkhoff[i].c_str());
		argv[3] = sFile;
		if ( i != 2)
			volumeAnswer = Valuation::mainValuationDriver(argv, 4);
		else
		{
			const char * argv2[5];
			argv2[0] = argv[0];
			argv2[1] = argv[1];
			argv2[2] = argv[2];
			argv2[3] = "--vrep";
			argv2[4] = sFile;
			volumeAnswer = Valuation::mainValuationDriver(argv2, 5);
		}//need to handel the v-rep file differently.
		delete[] sFile;

		RationalNTL correctVolumeAnswer(birkhoffVolume[i][0],
				birkhoffVolume[i][1]);
		VolumeTests::printVolumeTest(correctVolumeAnswer, volumeAnswer, string(birkhoff[i]),
				string("testing ") + string(birkhoff[i]));
	}//for ever file in the directory
}//runBirkhoffTests

