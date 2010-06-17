/*
 * Driver
 *
 *  Created on: June 15, 2010
 *      Author: Brandon
 */
#include <iostream>
#include "BuildRandomPolytope.h"
#include "BuildRandomEdgePolytope.h"
#include "GraphMaker.h"
#include "BuildGraphPolytope.h"



void doRandom()
{
	int size, points;
	cout << "Enter: size numPoints >> ";
	cin >> size >> points;
	BuildRandomPolytope rPoly(size);
	
	rPoly.buildPolymakeFile(points);
	//rPoly.callPolymake();
	rPoly.findEhrhardPolynomial();
}


void doEdge()
{
	int n, k;
	cout << "Enter: n (vector lenght) k (number of ones) >> ";
	cin >> n >> k;
	BuildRandomEdgePolytope rPoly(n, k);
	
	rPoly.buildPolymakeFile();
	//rPoly.callPolymake();
	rPoly.findEhrhardPolynomial();

}


void doGraphs()
{
	GraphMaker g;
	BuildGraphPolytope gPoly;
	int numVertx, offset;
	string type;

	cout << "graph type (l, cir, circ, or check, randC) >> ";
	cin >> type;

	if ( type == "l")
	{
		cout << " size >> ";
		cin >> numVertx;
		g.makeLinearGraph(numVertx);
	}//make edge poly
	else if (type == "cir")
	{
		cout << " size >> ";
		cin >> numVertx;
		g.makeCircleGraph(numVertx);
	}//if circle
	else if ( type == "circ")
	{
		cout << " size offset >> ";
		cin >> numVertx >> offset;
		g.makeCircleWithCenter(numVertx, offset);
	}//if circle with center
	else if ( type == "check")
	{
		cout << " row col >> ";
		cin >> numVertx >> offset; //I didn't want to bother making a better-named input var.
		g.makeCheckerboard(numVertx, offset);
	}//if checkerboard.
	else if ( type == "randC")
	{
		cout << "size edgeCount >> ";
		cin >> numVertx >> offset;
		g.makeRandomConnectedGraph(numVertx, offset);
	}//if random connected graph.



	g.printEdges();
	cout << "poly type (e or se) >> ";
	cin >> type;

	if (type == "e")
	{
		gPoly.buildPolymakeFile(g.getEdges(), BuildGraphPolytope::EDGE);
	}//edge poly
	else if (type == "se")
	{
		gPoly.buildPolymakeFile(g.getEdges(), BuildGraphPolytope::SYMMETRIC_EDGE);
	}//symmetric edge poly.

	gPoly.findEhrhardPolynomial();
}//doGraphs

int main()
{

	//doRandom();
	//doEdge();
	doGraphs();

	
	return 0;
}//main

