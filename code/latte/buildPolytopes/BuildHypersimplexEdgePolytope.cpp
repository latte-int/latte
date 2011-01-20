/*
 * BuildHypersimplexEdgePolytope.cpp
 *
 *  Created on: June 15, 2010
 *      Author: bedutra
 *
 *  Makes Hypersimplex edge polytopes
 */
#include "BuildHypersimplexEdgePolytope.h"
#include <iostream>

using namespace std;


BuildHypersimplexEdgePolytope::BuildHypersimplexEdgePolytope(): BuildPolytope()
{
}//BuildHypersimplexEdgePolytope


/**
 * To start things off, call with currentPoint = 11111....11000....00 (all ones, then all zeros) and addCurrent = true.
 * @parm list: a list of points found so far.
 * @parm currentPoint: the current point being processed.
 * @parm base:	keep elements in currentPont form 0 to base-1 alone
 * @parm addCurrent: checks to see if the currentPoint should be added to the list. The default is false.
 */
void BuildHypersimplexEdgePolytope::addToPoints(
	vector< vector<mpq_class> > &list,
	vector<mpq_class> currentPoint,
	int base, 
	bool addCurrent)
{
	if ( addCurrent == true)
		list.push_back(currentPoint);
	
	if ( base + 1 >= ambientDim)
		return; //numOnes must be = 1 and we are done.
	
	int lastOne = -1;
	//find the index of the last one or go the the end of the current point.
	for(int k = base; k < ambientDim && lastOne == -1; ++k)
		if ( currentPoint[k] == 0)
			lastOne = k - 1;
	if ( lastOne == -1)
	{
		if ( addCurrent == false)
			list.push_back(currentPoint);
		return; 
	}//if at end  ex: 000111
	//example: say current point = 110011100, assume base = 4, then lastOne = 6
	
	if ( lastOne - base == 0)
	{
		//example: if current point = 01011010000 and base= 6, then move the last one over a bunch of times
		for(int k = lastOne; k < ambientDim -1; ++k)
		{
			currentPoint[k] = 0;
			currentPoint[k + 1] = 1;
			list.push_back(currentPoint);
			
			//cout << "just added:";
			//for(int s = 0; s < (int) currentPoint.size(); ++s ) cout << currentPoint[s] << ' ';
			//cout << endl;
		}//for each shift
	
	}//move last one over and stop recursion
	else
	{
		addToPoints(list, currentPoint, base + 1); //don't add the current point again.
		//ex: if currentpoint = 001100 and base = 2, then find all permutations of "100" and prefix it to "001"

		//cout << "just finished processing point";
		//for(int s = 0; s < (int) currentPoint.size(); ++s ) cout << currentPoint[s] << ' ';
		//cout << " , going to start shifting" << endl;
			
		
		int numberOfShifts = (ambientDim - lastOne - 1);
		for(int k = 0; k < numberOfShifts; ++k)
		{
			currentPoint[base] = 0;
			for(int j = base + 1; j <= lastOne + 1; ++j)
				currentPoint[j] = 1;
				
			//cout << "shifted point:";
			//for(int s = 0; s < (int) currentPoint.size(); ++s ) cout << currentPoint[s] << ' ';
			//cout << endl;
			
			base = base + 1;
			lastOne = lastOne + 1;
			addToPoints(list, currentPoint, base + 1, true);
			//ex: if current point = 00111000 and base = 2,
			//    call addToPoints with base = 3, 4, 5 and
			//    currentPoint = 00011100, 00001110, 00000111
		}//for each shift, recurse.
	}//keep the recurison.
}//addToPoints

/**
 * Makes the vertices/points of the hypersimplex.
 * @parm ambient_dim: the ambient dim! What did you expect?
 * @numones: The number of 1's each vertex should have.
 */
void BuildHypersimplexEdgePolytope::generatePoints(int ambient_dim, int numones)
{
	ambientDim = ambient_dim;
	_numOnes = numones;

	clearPoints(); //base class
	vector<mpq_class> starter;
	vector<vector<mpq_class> > nonHomogenizedPoints;

	for(int k = 0; k < _numOnes; ++k)
		starter.push_back(1);
	for(int k = _numOnes; k < ambientDim; ++k)
		starter.push_back(0); // ex: starter = 1110000

	addToPoints(nonHomogenizedPoints, starter, 0, true);
	//nonHomogenizedPoints now contains all the points we need.

	for(int k = 0; k < (int) nonHomogenizedPoints.size(); ++k)
		addPoint(nonHomogenizedPoints[k]);//base class.
}




