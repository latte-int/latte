/*
 * BuildGraphPolytope.cpp
 *
 *  Created on: June 17, 2010
 *      Author: bedutra
 */
#include "BuildGraphPolytope.h"
#include <iostream>

using namespace std;


BuildGraphPolytope::BuildGraphPolytope(): BuildPolytope()
{
}//buildRandomPolytope


/**
 * add points using edge encoding: if vertex i, j are connected, add point e_i + e_j
 *   where e_i, e_j are std. normal directional vectors in R^(number of vertex in graph).
 *  @parm edges: an adj. list of edges. assuems an edge is only listed once.
 */
void BuildGraphPolytope::findEdgePolytope(const vector< vector<int> > &edges)
{
	ambientDim = edges.size();
	points.clear();

	for(int k = 0; k < edges.size(); ++k)
	{
		for(int j = 0; j < edges[k].size(); ++j)
		{
			vector<mpq_class> oneEdge(ambientDim, 0);
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
			vector<mpq_class> oneEdge(ambientDim, 0);//make an "ambientDim' long vector filled with zeros.
			oneEdge[k] = 1;
			oneEdge[edges[k][j]] = -1;
			points.push_back(oneEdge);
			oneEdge[k] = -1;
			oneEdge[edges[k][j]] = 1;
			points.push_back(oneEdge);
		}//for each edge.
	}//for each row
}//findSymmetricEdgePolytope


