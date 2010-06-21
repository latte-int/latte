/*
 * Driver
 *
 *  Created on: June 15, 2010
 *      Author: Brandon
 */
#include <iostream>
#include <sstream>
#include "BuildRandomPolytope.h"
#include "BuildRandomEdgePolytope.h"
#include "GraphMaker.h"
#include "BuildGraphPolytope.h"



void doRandom()
{
	int size, points;
	stringstream comments;

	cout << "Enter: size numPoints >> ";
	cin >> size >> points;
	BuildRandomPolytope rPoly(size);
	
	rPoly.buildPolymakeFile(points);
	//rPoly.callPolymake();
	comments << "Random edge polytope with " << points << " random points of size " << size << ". ";
	rPoly.setComments(comments.str());
	rPoly.findEhrhardPolynomial();
}


void doEdge()
{
	int n, k;
	stringstream comments;

	cout << "Enter: n (vector lenght) k (number of ones) >> ";
	cin >> n >> k;
	BuildRandomEdgePolytope rPoly(n, k);
	
	rPoly.buildPolymakeFile();
	//rPoly.callPolymake();

	comments << "Hypersimplex with n=" << n << ", and k=" << k << ". ";
	rPoly.setComments(comments.str());
	rPoly.findEhrhardPolynomial();

}


void doGraphs()
{
	GraphMaker g;
	BuildGraphPolytope gPoly;
	int numVertx, offset;
	string type;
	stringstream comments;

	cout << "graph type (l, cir, circ, or check, randC) >> ";
	cin >> type;

	if ( type == "l")
	{
		cout << " size >> ";
		cin >> numVertx;
		g.makeLinearGraph(numVertx);
		comments << "List graph with" << numVertx << " nodes. ";
	}//make edge poly
	else if (type == "cir")
	{
		cout << " size >> ";
		cin >> numVertx;
		g.makeCircleGraph(numVertx);
		comments << "Circle graph width " << numVertx << " nodes. ";
	}//if circle
	else if ( type == "circ")
	{
		cout << " size offset >> ";
		cin >> numVertx >> offset;
		g.makeCircleWithCenter(numVertx, offset);
		comments << "Center-Circle graph with " << numVertx << " nodes and offset " << offset << ". ";
	}//if circle with center
	else if ( type == "check")
	{
		cout << " row col >> ";
		cin >> numVertx >> offset; //I didn't want to bother making a better-named input var.
		g.makeCheckerboard(numVertx, offset);
		comments << "Checkerboard graph of size " << numVertx << " by " << offset << ". ";
	}//if checkerboard.
	else if ( type == "randC")
	{
		cout << "size edgeCount >> ";
		cin >> numVertx >> offset;
		g.makeRandomConnectedGraph(numVertx, offset);
		comments << "Random graph with " << numVertx << " nodes and " << offset << " many edges. ";
	}//if random connected graph.
	else
	{
		cout << "ops, that is not a graph type, calling exit" << endl;
		exit(1);
	}//else stop



	g.printEdges();
	cout << "poly type (e or se) >> ";
	cin >> type;

	if (type == "e")
	{
		gPoly.buildPolymakeFile(g.getEdges(), BuildGraphPolytope::EDGE);
		comments << " Edge polytope encoding used. ";
	}//edge poly
	else if (type == "se")
	{
		gPoly.buildPolymakeFile(g.getEdges(), BuildGraphPolytope::SYMMETRIC_EDGE);
		comments << "Symmetric edge polytope encoding used. ";
	}//symmetric edge poly.
	else
	{
		cout << "ops, that is not a polytope type, calling exit" << endl;
		exit(1);
	}//else stop

	gPoly.setComments(comments.str());
	gPoly.findEhrhardPolynomial();
}//doGraphs

void doAuto()
{
	stringstream comments;

/*	
	for(int i =  20; i <= 40; ++i)
	{
		stringstream comments;
		BuildGraphPolytope gPoly;
		GraphMaker g;
		
		comments << "doing line graph of size " << i << " with edge encoding" ;
		cout << "***********************************\n" << comments.str() << endl;
		g.makeLinearGraph(i);
		gPoly.buildPolymakeFile(g.getEdges(), BuildGraphPolytope::EDGE);
		gPoly.setComments(comments.str());
		gPoly.findEhrhardPolynomial();
	}//line graphs.
*/	


	for(int i = 15; i <= 30; ++i)
		for(int j = i; j <= i+10; ++j)
		{
			stringstream comments;
			BuildGraphPolytope gPoly;
			GraphMaker g;
			
			comments << "doing random graph of size " << i << " and edge " << j << " with edge encoding" ;
			cout << "***********************************\n" << comments.str() << endl;
			g.makeRandomConnectedGraph(i, j);
			gPoly.buildPolymakeFile(g.getEdges(), BuildGraphPolytope::EDGE);
			gPoly.setComments(comments.str());
			gPoly.findEhrhardPolynomial();
		}//rand graph
		
		
		

	for(int i = 15; i <= 30; ++i)
		for(int j = i; j <= i+10; ++j)
		{
			stringstream comments;
			BuildGraphPolytope gPoly;
			GraphMaker g;
			
			comments << "doing random graph of size " << i << " and edge " << j << " with symmetric edge encoding" ;
			cout << "***********************************\n" << comments.str() << endl;
			g.makeRandomConnectedGraph(i, j);
			gPoly.buildPolymakeFile(g.getEdges(), BuildGraphPolytope::SYMMETRIC_EDGE);
			gPoly.setComments(comments.str());
			gPoly.findEhrhardPolynomial();
		}//rand graph		
}//doAuto

int main()
{
	string type;
	
	doAuto();
	return 0;

	cout << "run type: (rand, edge, graph) >> ";
	cin >> type;
	if ( type == "rand")
		doRandom(); //call polymake on random points
	else if ( type == "edge")
		doEdge();   //call polymake on hypersimplices
	else if ( type == "graph")
		doGraphs(); //call polymake on graph polytopes
	else if ( type == "auto")
		doAuto();
	else	
		exit(1); //error.

	
	return 0;
}//main

