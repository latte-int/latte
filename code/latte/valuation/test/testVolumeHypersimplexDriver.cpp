/*
 * testVolumeHypersimplex.cpp
 *
 *  Created on: Aug 10, 2010
 *      Author: Brandon Dutra and Gregory Pinto
 *  Computes the volume of the standard simplex in many dimensions.
 */

#include "VolumeAndIntegrationTests.h"


int main(int argc, char *argv[])
{
	VolumeTests::runHyperSimplexTests();
	//VolumeTests::runBirkhoffTests();
	return 0;
}
//main()
