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


/**
 * @parm edges: an adj. list of edges.
 * @parm polytopeType: defines which encoding should be used.
 * Assumes: edges is non empty, otherwise will make an empty polymake file.
 */
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


/**
 * add points using edge encoding: if vertex i, j are connected, add point e_i + e_j
 *   where e_i, e_j are std. normal directional vectors in R^(number of vertex in graph).
 */
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



/**
 * add points using symmetric edge encoding: if vertex i, j are connected, add point e_i - e_j and e_j - e_i
 *   where e_i, e_j are std. normal directional vectors in R^(number of vertex in graph).
 */
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


