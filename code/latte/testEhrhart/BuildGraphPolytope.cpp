/*
 * BuildGraphPolytope.cpp
 *
 *  Created on: June 17, 2010
 *      Author: bedutra
 */
#include "BuildGraphPolytope.h"
#include <iostream>

using namespace std;


BuildGraphPolytope::BuildGraphPolytope(): BuildRandomPolytope(0)
{
	PolytopeComments = "A graph polytope. "; //caller should call setComments for a better label.
}//buildRandomPolytope

void BuildGraphPolytope::buildPolymakeFile(const vector< vector<int> > &edges, PolytopeType polytopeType)
{

	switch(polytopeType)
	{
		case EDGE: findEdgePolytope(edges); break;
		case SYMMETRIC_EDGE: findSymmetricEdgePolytope(edges); break;
		default:
			cout << "BuildGraphPolytope::buildPolymakeFile(): ops, unknown PolyTopeType" << endl;
			break;
	}//switch


	ofstream file;

	file.open(fileName.c_str());
	file << "POINTS" << endl;

	for(int k = 0; k < points.size(); ++k)
	{
		file << 1 << ' ';
		for(int vectorElm = 0; vectorElm < ambientDim; ++vectorElm)
		{
			file << points[k][vectorElm] << ' ';
			//cout << points[k][vectorElm] << ' ';
		}//vectorElm
		file << endl;
		//cout << endl;
	}//for k
	file.close();
}//BuildRandomPolytope


void BuildGraphPolytope::findEdgePolytope(const vector< vector<int> > &edges)
{
	ambientDim = edges.size();
	points.clear();

	for(int k = 0; k < edges.size(); ++k)
	{
		for(int j = 0; j < edges[k].size(); ++j)
		{
			vector<int> oneEdge(ambientDim, 0);
			oneEdge[k] = 1;
			oneEdge[edges[k][j]] = 1;
			points.push_back(oneEdge);
		}//for each edge.
	}//for each row
}//findEdgePolytope


void BuildGraphPolytope::findSymmetricEdgePolytope(const vector< vector<int> > &edges)
{
	ambientDim = edges.size();
	points.clear();

	for(int k = 0; k < edges.size(); ++k)
	{
		for(int j = 0; j < edges[k].size(); ++j)
		{
			vector<int> oneEdge(ambientDim, 0);//make an "ambientDim' long vector filled with zeros.
			oneEdge[k] = 1;
			oneEdge[edges[k][j]] = -1;
			points.push_back(oneEdge);
			oneEdge[k] = -1;
			oneEdge[edges[k][j]] = 1;
			points.push_back(oneEdge);
		}//for each edge.
	}//for each row
}//findSymmetricEdgePolytope


