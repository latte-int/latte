/*
 * GraphMaker.cpp
 *
 *  Created on: June 16, 2010
 *      Author: bedutra
 */

#include <iostream>
#include <algorithm>

#include "GraphMaker.h"

using namespace std;


GraphMaker::GraphMaker()
{
	numVertex = 0; 
	srand(time(0));
}//graphMaker()



/**
 * Adds the edge (v1, v2) or (v2, v1) to the graph if the edge does not already exist.
 * @return bool: true of new edge added, false if edge already exists.
 */
bool GraphMaker::addEdgeInOrder(const int v1, const int v2)
{
	if (v1 <= v2 )
	{
		if (find(edges[v1].begin(), edges[v1].end(), v2) == edges[v1].end())
		{
			edges[v1].push_back(v2);
			return true;
		}//not found. so add it.
		else
			return false;
	}
	else
		return addEdgeInOrder(v2, v1);
}//addEdgeInOrder

/**
 * Returns a const reference to the edges structure.
 */
const vector< vector<int>  > & GraphMaker::getEdges() const
{
	return edges;
}//getEdges


/**
 * Makes a rectangle graph:
 *   0  1  2  3
 *   4  5  6  7
 *   8  9 10 11
 *   Each node is connected to its neighbors to the top/bottom and left/right with no wrapping (4 and 7 are not connected but 4 and 5 are).
 */
void GraphMaker::makeCheckerboard(const int row, const int col)
{

	if ( row <= 1 || col <= 1)
	{
		cout << "makeLinearGraph(): please give a row/col larger than 1" << endl;
		return;
	}//if crapy size.

	numVertex = row * col;
	edges.clear();
	edges.resize(numVertex);
	for(int k = 0; k < numVertex; ++k) edges[k].clear();

	//How to think about edges: consider a row x col matrix (with index starting at zero!)
	//	we save this 2D matrix in 1D array using row-first order. So the (a, b) element of
	//  the matrix (again, index starts at 0) is saved in location row*a + b

	//first process a (row-1)x(col -1) grid and then fix the right and lower edges.
	for(int currentRow = 0; currentRow < row - 1; ++currentRow )
	{
		for(int currentCol = 0; currentCol < col - 1; ++currentCol)
		{
			int rowIndex = currentRow * col;

			edges[rowIndex + currentCol].push_back(rowIndex + currentCol + 1);
			edges[rowIndex + currentCol].push_back((currentRow + 1)*col + currentCol);
		}//for currentCol
	}//for currentRow

	//now fix the right edge
	for(int currentRow = 0; currentRow < row -1; ++currentRow)
	{
		edges[currentRow*col + col - 1].push_back((currentRow + 1)*col + col -1);
	}//for currentRow

	//now fix the bottom edge.
	for(int currentCol = 0; currentCol < col -1; ++currentCol)
	{
		edges[(row -1)*col + currentCol].push_back((row - 1)*col + currentCol + 1);
	}//for currentCol
}//makeCheckerboard()


/**
 *	Make a graph in the form 0-1-2-3-4-5-0 (so it loops around) and place a edge
 *  	from every offset many to an other center node number 6.
 *
 *  @parm size: number of vertices in the graph
 *  @parm offset: over offset many vertices will be connected to the center node (if offset = 1,
 *  	then every node on the outer circle will be connected to the center node!)
 */
void GraphMaker::makeCircleWithCenter(const int size, const int offset)
{

	if ( size <= 3)
	{
		cout << "makeLinearGraph(): please give a size larger than 3" << endl;
		return;
	}//if crapy size.

	numVertex = size;
	edges.clear();
	edges.resize(numVertex);
	for(int k = 0; k < numVertex; ++k) edges[k].clear();

//cout << "got here 1" << endl;
//cout << "edges.size()=" << edges.size() << endl;
//cout << "edges.cap() =" << edges.capacity() << endl;
	for(int k = 0; k < numVertex - 2; k++)
	{
		//cout << "k=" << k << "numVertex=" << numVertex << endl;
		cout << "edges[k].zize() = " << edges[k].size() << endl;
		edges[k].push_back(k + 1);
	}
//cout << "got here 1.2" << endl;
	edges[0].push_back(numVertex -2);
//cout << "got here 2" << endl;
	for(int k = 0; k < numVertex - 1; k++)
	{
		if ( k % offset == 0)
		{
			edges[k].push_back(numVertex - 1);
		}//add edge.
	}//for each node on the circle.

}//makeCircleWithCenter()

/**
 * Make a graph in the form 0-1-2-3-4-5-0 (so it loops around)
 */
void GraphMaker::makeCircleGraph(const int size)
{

	if ( size <= 2)
	{
		cout << "makeLinearGraph(): please give a size larger than 2" << endl;
		return;
	}//if crapy size.

	numVertex = size;
	edges.clear();
	edges.resize(numVertex);
	for(int k = 0; k < numVertex; ++k) edges[k].clear();

	for(int k = 0; k < numVertex - 1; k++)
	{
		edges[k].push_back(k + 1);
	}

	edges[0].push_back(numVertex -1);
}//makeCircleGraph()


/**
 * Make a graph in the form: 0-1-2-3-4-5
 */
void GraphMaker::makeLinearGraph(const int size)
{
	if ( size <= 1)
	{
		cout << "makeLinearGraph(): please give a size larger than 1" << endl;
		return;
	}//if crapy size.
	
	numVertex = size;
	edges.clear();
	edges.resize(numVertex);
	for(int k = 0; k < numVertex; ++k) edges[k].clear();
	
	for(int k = 0; k < numVertex - 1; k++)
	{		
		edges[k].push_back(k + 1);
	}

}//makeLinearGraph


void GraphMaker::makeRandomConnectedGraph(const int size, const int edgeCount)
{
	int v1, v2, currentEdgeCount = 0;

	if ( size <= 2 || size > edgeCount+1 || edgeCount > size*(size -1)/2)
	{
		cout << "makeLinearGraph(): please give a size larger than 2 or an edgeCount >= size or you have too many edges" << endl;
		return;
	}//if bad size.

	numVertex = size;
	edges.clear();
	edges.resize(numVertex);
	for(int k = 0; k < numVertex; ++k) edges[k].clear();

	//first, make a random min. spanning tree.
	makeRandomSpanningTree();
	currentEdgeCount = numVertex -1;

	cout << "spanning tree:" << endl;
	printEdges();
	//cout << "now filling rest of quota" << endl;

	while ( currentEdgeCount < edgeCount)
	{
		//pick two vertices and try them.
		v1 = rand() % numVertex;
		v2 = rand() % numVertex;
		if ( v1 == v2)
			continue; //try again.

		//cout << "going to try v1, v2=" << v1 << ", " << v2 << endl;
		//printEdges();
		//cout << endl;
		if ( addEdgeInOrder(v1, v2) == true)
			++currentEdgeCount;

	}//while, keep adding more edges.

	//cout << "full graph" << endl;
	//printEdges();
	
}//makeRandomConnectedGraph()


void GraphMaker::makeRandomSpanningTree()
{
	vector<int> unusedVertices(numVertex-1);
	vector<int> usedVertices;
	int currentUsedVertex, currentUnusedVertexIndex, lastUnused = numVertex - 2;

	for(int i = 0; i < numVertex-1; ++i)
		unusedVertices[i] = i;
	usedVertices.push_back(numVertex -1);

	while(lastUnused >= 0)
	{
		currentUsedVertex = usedVertices[rand() % usedVertices.size()];
		currentUnusedVertexIndex = rand() % (lastUnused +1);

		//move the vertex we wish to add to the end of the unusedVertices vector (by end, I mean to lastUnused).
		// 	lastUnused will be decremented later, so we will never see it again!
		swap(unusedVertices[currentUnusedVertexIndex], unusedVertices[lastUnused]);
		usedVertices.push_back(unusedVertices[lastUnused]);

		//now add the edge (currentusedVertex, unusedVertices[lastUnused]) to the graph, saving the smaller in the outer vertex.
		addEdgeInOrder(currentUsedVertex, unusedVertices[lastUnused]);
		/*
		if ( currentUsedVertex <= unusedVertices[lastUnused])
		{
			edges[currentUnusedVertex].push_back(unusedVertices[lastUnused]);
		}//add (currentusedVertex, unusedVertices[lastUnused]) to edges
		else
		{
			edges[unusedVertices[lastUnused]].push_back(currentUsedVertex);
		}//add (unusedVertices[lastUnused], currentusedVertex) to edges
		*/
		--lastUnused;
	}//while
}//makeRandomSpanningTree()


void GraphMaker::printEdges() const
{
	cout << "numVertex=" << numVertex << endl;
	for(int i = 0; i < numVertex; ++i)
	{
		for(int k = 0; k < (int) edges[i].size(); ++k)
			cout << "( " << i << ", " << edges[i][k] << ")" << endl;
	}//for i

}//printEdges

