/*
 * BuildRandomPolytope.cpp
 *
 *  Created on: June 15, 2010
 *      Author: Brandon Dutra and Gregory Pinto
 */
#include "BuildRandomPolytope.h"
#include <iostream>
#include <cstdlib> //for RAND_MAX

using namespace std;

BuildRandomPolytope::BuildRandomPolytope(): BuildPolytope()
{
	_maxInteger = 50;
	integerPoints = true;
	_probNegative = .5;
}//buildRandomPolytope

void BuildRandomPolytope::makePoints(int ambient_dim, int numberPoints, int maxInt, double probNeg)
{
	_maxInteger = maxInt;
	_probNegative = probNeg;
	makePoints(ambient_dim, numberPoints);
}

/**
 * makes numberPoints many points in Q^(ambient_dim)
 */
void BuildRandomPolytope::makePoints(int ambient_dim, int numberPoints)
{
	assert(numberPoints >= ambient_dim +1);
	ambientDim = ambient_dim;

	clearPoints();//base class
	for (int i = 0; i < numberPoints; ++i)
	{
		vector<mpq_class> onePoint;

		onePoint.resize(ambientDim);
		for(int j = 0; j < ambientDim; ++j)
		{
			if ( integerPoints == true)
				onePoint[j] = mpq_class(rand()%_maxInteger, 1);
			else
				onePoint[j] = mpq_class(rand()%_maxInteger, (rand()%_maxInteger) +1);

			if ( rand() < RAND_MAX * _probNegative)
				onePoint[j] *= -1;//make negative.
		}//build one point.
		addPoint(onePoint);//base class.

		//for(int j = 0; j <(int)onePoint.size(); ++j)
		//	cout << onePoint[j] << ", ";
		//cout << endl;
	}//for each point.
}//makePoints
