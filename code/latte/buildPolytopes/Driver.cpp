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


void findEhrhartPolynomial(const string &file)
{
	system((string("./count --ehrhart-polynomial --maxdet=1000000 ") + file).c_str());
}

void doRandom()
{
	int size, points;
	stringstream comments;

	cout << "Enter: size numPoints >> ";
	cin >> size >> points;
	BuildRandomPolytope rPoly;
	rPoly.makePoints(size, points);

	rPoly.buildLatteHRepFile();

	findEhrhartPolynomial(rPoly.getLatteHRepFile());
	rPoly.deleteLatteHRepFile();
	rPoly.deletePolymakeFile();
}


void doHypersimplex()
{
	int n, k;
	stringstream comments;

	cout << "Enter: n (vector lenght) k (number of ones) >> ";
	cin >> n >> k;
	BuildHypersimplexEdgePolytope rPoly;
	rPoly.generatePoints(n, k);
	
	rPoly.buildLatteHRepFile();

	findEhrhartPolynomial(rPoly.getLatteHRepFile());

	rPoly.deleteLatteHRepFile();
	rPoly.deletePolymakeFile();
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



	//g.printEdges();
	cout << "edge type: regular edge or symmetric edge (e or se) >> ";
	cin >> type;

	if (type == "e")
	{
		gPoly.findEdgePolytope(g.getEdges());
		comments << " Edge polytope encoding used. ";
	}//edge poly
	else if (type == "se")
	{
		gPoly.findSymmetricEdgePolytope(g.getEdges());
		comments << "Symmetric edge polytope encoding used. ";
	}//symmetric edge poly.
	else
	{
		cout << "ops, that is not a polytope type, calling exit" << endl;
		exit(1);
	}//else stop

	gPoly.buildLatteHRepFile();
	findEhrhartPolynomial(gPoly.getLatteHRepFile());

	gPoly.deleteLatteHRepFile();
	gPoly.deletePolymakeFile();
}//doGraphs



int main()
{
	string type;


	cout << "run type: (rand, hypersimplex, graph) >> ";
	cin >> type;
	if ( type == "rand")
		doRandom(); //call polymake on random points
	else if ( type == "hypersimplex")
		doHypersimplex();   //call polymake on hypersimplices
	else if ( type == "graph")
		doGraphs(); //call polymake on graph polytopes
	else	
		exit(1); //error.

	
	return 0;
}//main

