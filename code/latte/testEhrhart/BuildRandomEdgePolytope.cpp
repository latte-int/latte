/*
 * PolynomialInterpolation.cpp
 *
 *  Created on: June 15, 2010
 *      Author: bedutra
 */
#include "BuildRandomEdgePolytope.h"
#include <iostream>

using namespace std;


BuildRandomEdgePolytope::BuildRandomEdgePolytope(int ambient_dim, int numones): BuildRandomPolytope(ambient_dim), numOnes(numones)
{

}//buildRandomEdgePolytope


void BuildRandomEdgePolytope::addToPoints(
	vector< vector<int> > &list,
	vector<int> currentPoint, 
	int base, 
	bool addCurrent) const
{
	if ( addCurrent == true)
		list.push_back(currentPoint);
	
	if ( base + 1 >= ambientDim)
		return; //numOnes must be = 1 and we are done.
	
	int lastOne = -1;
	for(int k = base; k < ambientDim && lastOne == -1; ++k)
		if ( currentPoint[k] == 0)
			lastOne = k - 1;
	if ( lastOne == -1)
	{
		if ( addCurrent == false)
			list.push_back(currentPoint);
		return; 
	}//if at end  ex: 000111
	
	
	
	if ( lastOne - base == 0)
	{
		for(int k = lastOne; k < ambientDim -1; ++k)
		{
			currentPoint[k] = 0;
			currentPoint[k + 1] = 1;
			list.push_back(currentPoint);
			
			cout << "just added:";
			for(int s = 0; s < (int) currentPoint.size(); ++s ) cout << currentPoint[s] << ' ';
			cout << endl;
		}//for each shift
	
	}//move last one over and stop recursion
	else
	{
		addToPoints(list, currentPoint, base + 1);
		cout << "just finished processing point";
		for(int s = 0; s < (int) currentPoint.size(); ++s ) cout << currentPoint[s] << ' ';
		cout << " , going to start shifting" << endl;
			
		
		int numberOfShifts = (ambientDim - lastOne - 1);
		for(int k = 0; k < numberOfShifts; ++k)
		{
			currentPoint[base] = 0;
			for(int j = base + 1; j <= lastOne + 1; ++j)
				currentPoint[j] = 1;
				
			cout << "shifted point:";
			for(int s = 0; s < (int) currentPoint.size(); ++s ) cout << currentPoint[s] << ' ';
			cout << endl;
			
			base = base + 1;
			lastOne = lastOne + 1;
			addToPoints(list, currentPoint, base + 1, true);

		}//for each shift, recurse.
	}//keep the recurison.
	

}//addToPoints


void BuildRandomEdgePolytope::buildPolymakeFile()
{
	ofstream file;
	
	file.open(fileName.c_str());
	file << "POINTS" << endl;
	
	vector< vector<int> > points;
	
	vector<int> starter;
	for(int k = 0; k < numOnes; ++k)
		starter.push_back(1);
	for(int k = numOnes; k < ambientDim; ++k)
		starter.push_back(0);
	
	addToPoints(points, starter, 0, true);
	for(int k = 0; k < (int) points.size(); k++)
	{
		for(int i = 0; i < ambientDim; ++i)
			cout << points[k][i] << " ";
		cout << endl;
	}//for k
	
	
	
	file.close();
	
}//BuildRandomPolytope



