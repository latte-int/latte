/*
 * GraphMaker.cpp
 *
 *  Created on: June 16, 2010
 *      Author: bedutra
 */



#include "GraphMaker.h"
#include "BuildHypersimplexEdgePolytope.h" //needed to make Kneser Graph.

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
 *	Kneser graph KG(n,k) is the graph whose vertices correspond to the k-element subsets of a set of n elements, and where two vertices are connected if and only if the two corresponding sets are disjoint
 */
void GraphMaker::makeKneserGraph(const int setSize, const int subSetSize)
{
	BuildHypersimplexEdgePolytope hPoly;
	vector< vector<mpq_class> > points;
	vector<mpq_class> startingPoint;
	bool match;

	hPoly.generatePoints(setSize, subSetSize);

	if (subSetSize > (setSize +1)/2)
	{
		cout << "subsetSize too large, no edges can exist" << endl;
		return;
	}//too large.

	numVertex = nchoosek(setSize, subSetSize);
	edges.clear();
	edges.resize(numVertex);

	if ( setSize <= subSetSize)
	{
		cout << "subset size is too big" << endl;
		return;
	}//if too big.

	//set the starting point to 11111....11000....0.
	startingPoint.resize(setSize);
	for(int i = 0; i < subSetSize; ++i)
		startingPoint[i] = 1;
	for(int i = subSetSize; i < setSize; ++i)
		startingPoint[i] = 0;

	hPoly.addToPoints(points, startingPoint, 0, true);
	//ok, points now contains all subsets of size subSetSize from a set of setSize.

//	for(int i = 0; i < points.size(); ++i)
//	{
//		cout << i << ") ";
//		for(int k = 0; k < setSize; ++k)
//			cout << " " << points[i][k];
//		cout << endl;
//	}//print points out.

	for(size_t i = 0; i < points.size(); ++i)
		for(size_t k = i + 1; k < points.size(); ++k)
		{
			match = false;
			for(int curElement = 0; curElement < setSize && ! match; ++curElement)
				if ( points[i][curElement] == 1 && points[i][curElement] == points[k][curElement])
					match = true;

			if ( match == false)
			{
				//cout << "match = " << match << endl;
				edges[i].push_back(k);
			}
		}//for k. check to see if point[i] and point[k] intercept.
}//makeKneserGraph



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


/**
 * Makes the classical well-known Petersen Graph (10 vertices) starting at vertex 0.
 */
void GraphMaker::makePetersenGraph()
{
	numVertex = 10;
	edges.clear();
	edges.resize(numVertex);

	makePetersenSubGraph(0);
}//makePertersonGraph

/**
 * This function could change over time.
 */
void GraphMaker::makePetersenFunGraph(const int num)
{
	edges.clear();
	numVertex = 10 * num;
	edges.resize(10 * num);

	//make a couble disjoing petersen graphs.
	for(int i = 0; i < num; ++i)
		makePetersenSubGraph(i*10);

}//makePetersenFunGraph

/**
 *  Makes a petersen graph starting at vertex startVertex
 *  This is helpful if you want to add a Petersen graph to an existing graph.
 */
void GraphMaker::makePetersenSubGraph(const int startVertex)
{
	for(int i = startVertex; i < startVertex + 4; ++i)
		edges[i].push_back(i+1);
	edges[startVertex].push_back(startVertex + 4);

	for(int i = startVertex; i < startVertex + 5; ++i)
		edges[i].push_back(i + 5);
	for(int i = startVertex + 5; i < startVertex + 8; ++i)
		edges[i].push_back(i + 2);
	for(int i = startVertex + 5; i < startVertex + 7; ++i)
		edges[i].push_back(i + 3);
}//makePertersonSubGraph

/**
 * Will make two (simple) graphs of size/2 with edgeCount/2 edges that are not connected to each other.
 * These two graphs may or may not further be connected.
 */
void GraphMaker::makeRandomDisconnectedGraph(const int size, const int edgeCount)
{

	int vertexG1, vertexG2, edgeG1, edgeG2; //size and edgeCount of graph 1 and 2
	int v1, v2, currentEdgeCount;

	if ( size < 4)
	{
		cout << "please give a size larger than 4";
		return;
	}

	numVertex = size;
	edges.clear();
	edges.resize(numVertex);

	//find the size and edgeCount of the first graph.
	vertexG1 = (int)((size + 1 )/2); //round up.
	vertexG2 = size / 2;

	//find the size and edgeCount of the 2nd graph.
	edgeG1 = (int)((edgeCount + 1 )/2); //round up.
	edgeG2 = edgeCount / 2;

	cout << vertexG1 << "::" << edgeG1 << ", " << vertexG2 << "::" << edgeG2 << endl;
	currentEdgeCount = 0;
	while ( currentEdgeCount < edgeG1)
	{
		v1 = rand() % vertexG1;
		v2 = rand() % vertexG1;
		if ( v1 == v2) continue; //try again. Keep the graph simple.

		if ( addEdgeInOrder(v1, v2) == true)
		{
			//cout << "1) v1=" << v1 << ", v2=" << v2 << endl;
			++currentEdgeCount;
		}
	}//while not done making graph.

//cout << "got here" << endl;
	currentEdgeCount = 0;
	while ( currentEdgeCount < edgeG2)
	{
		v1 = (rand() % vertexG2) + vertexG1;
		v2 = (rand() % vertexG2) + vertexG1 ;

		if ( v1 == v2) continue; //try again. Keep the graph simple.

		if ( addEdgeInOrder(v1, v2) == true)
		{
			//cout << "2) v1=" << v1 << ", v2=" << v2 << endl;
			++currentEdgeCount;
		}
	}//while not done making graph.

}//makeRandomDisconnectedGraph

/**
 * Makes a random connected simple graph.
 */
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


/**
 * Make sure the graph is connected.
 */
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
		--lastUnused;
	}//while
}//makeRandomSpanningTree()


/*
 * Ask the user for input. Does not require the end graph to be connected or simple.
 */
void GraphMaker::makeYourOwnGraph()
{
	int v1, v2;

	cout << "One past the largest graph vertex number >> ";
	cin >> numVertex;

	edges.clear();
	edges.resize(numVertex);
	for(int i = 0; i < numVertex; ++i) edges[i].clear();

	while (true)
	{
		cout << "Enter -1 or vertex_1 vertex_2 >> ";
		cin >> v1;
		if ( v1 == -1)
			break;
		cin >> v2;

		if ( v1 >= numVertex || v2 >= numVertex || v1 < 0 || v2 < 0)
		{
			cout << "vertex number too large or too small" << endl;
			continue;
		}//if too large/small

		if ( false == addEdgeInOrder(v1, v2))
			cout << "That edge already exists" << endl;
	}//add another edge.

}//makeYourOwnGraph()


/**
 * returns "n choose k"
 */
int GraphMaker::nchoosek(const int n1, const int k1)
{
	mpz_class n(n1), k(k1), topProduct(1), bottomProduct(1);
	mpq_class ans;
	if ( n1 < k1 )
	{
		cout << "nchoosek() bad input" << endl;
		return -1;
	}//error

	for(mpz_class i(0); i < k; ++i)
	{
		topProduct *= ( n - i);
	}
	for(mpz_class i(1); i <= k; ++i)
	{
		bottomProduct *= i;
	}

	ans = mpq_class(topProduct, bottomProduct);
	ans.canonicalize();

	return ((int) ans.get_num().get_si());
}//nchoosek()

/**
 * Prints the current edges.
 */
void GraphMaker::printEdges() const
{
	cout << "numVertex=" << numVertex << endl;
	for(int i = 0; i < numVertex; ++i)
	{
		for(int k = 0; k < (int) edges[i].size(); ++k)
			cout << "( " << i << ", " << edges[i][k] << ")" << endl;
	}//for i

}//printEdges

