/*
 * Driver
 *
 *  Created on: June 15, 2010
 *      Author: Brandon
 */
#include <iostream>
#include <sstream>
#include "BuildRandomPolytope.h"
#include "BuildHypersimplexEdgePolytope.h"
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
	BuildHypersimplexEdgePolytope rPoly(n, k);
	
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

	cout << "graph type \n(l, cir, circ, or check, rand[CD], petersen, petersenFun, my, kneser) >> ";
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
		comments << "Random connected graph with " << numVertx << " nodes and " << offset << " many edges. ";
	}//if random connected graph.
	else if (type == "randD")
	{
		cout << "size edgeCount >> ";
		cin >> numVertx >> offset;
		g.makeRandomDisconnectedGraph(numVertx, offset);
		comments << "Random disconnected graph with " << numVertx << " nodes and " << offset << " many edges. ";
	}//if random dis-connnected graph.
	else if ( type == "petersen")
	{
		g.makePetersenGraph();
		comments << "Petersen graph.";
	}//if petersen
	else if (type == "petersenFun")
	{
		cout << "num >> ";
		cin >> numVertx;
		g.makePetersenFunGraph(numVertx);
		comments << "A Fun Peterson Graph";
	}//if petersen fun (need to read source, this funtion will change now and then)
	else if (type == "my")
	{
		g.makeYourOwnGraph();
		comments << "User defined edge-graph. ";
	}//if user defined graph.
	else if (type == "kneser")
	{
		cout << " setSize subSetSize >> ";
		cin >> numVertx >> offset;
		g.makeKneserGraph(numVertx, offset);
		comments << "Kneser graph with " << offset << " in " << numVertx << ".";
	}// if kneser
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

//try two disconnected petersen graphs.
	stringstream comments2;
	BuildGraphPolytope gPoly2;
	GraphMaker g2;
	g2.makePetersenFunGraph(2);
	comments2 << "A Fun Peterson Graph";
	gPoly2.buildPolymakeFile(g2.getEdges(), BuildGraphPolytope::EDGE);
	gPoly2.setComments(comments2.str());
	gPoly2.findEhrhardPolynomial();

	for(int i = 5; i <= 20; i = i + 1)
		for(int j = 1; j <= i/2; j = j + 1)
		{
			stringstream comments;
			BuildGraphPolytope gPoly;
			GraphMaker g;

			comments << "Kneser graph with " << j << " < " << i << ". ";
			cout << "***********************************\n" << comments.str() << endl;
			g.makeKneserGraph(i, j);
			g.printEdges();
			gPoly.buildPolymakeFile(g.getEdges(), BuildGraphPolytope::EDGE);
			gPoly.setComments(comments.str());
			gPoly.findEhrhardPolynomial();
		}//rand graph


	//already did i=[15,30], j=[i, i+10] untill i = 17, j= 24 for connected edge graphs.
	
	
	for(int i = 35; i <= 40; i = i + 5)
		for(int j = i+3; j <= i+9; j = j + 3)
		{
			stringstream comments;
			BuildGraphPolytope gPoly;
			GraphMaker g;
			
			comments << "Random disconnected graph with " << i << " nodes and " << j << " edges. ";
			cout << "***********************************\n" << comments.str() << endl;
			g.makeRandomDisconnectedGraph(i, j);
			g.printEdges();
			gPoly.buildPolymakeFile(g.getEdges(), BuildGraphPolytope::EDGE);
			gPoly.setComments(comments.str());
			gPoly.findEhrhardPolynomial();
		}//rand graph


/*
	for(int i = 3; i <= 30; i = i + 5)
		for(int j = i; j <= i+0; ++j)
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
*/		


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

