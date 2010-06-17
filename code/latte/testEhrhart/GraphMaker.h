/*
 *
 *
 *  Created on: June 16, 2010
 *      Author: bedutra
 */

#ifndef GRAPH_MAKER_H_
#define GRAPH_MAKER_H_

#include <iomanip>
#include <string>
#include <fstream>
#include <ctime>
#include <vector>
#include <sstream>
#include <cstdlib>
#include "gmp.h"
#include <gmpxx.h>



using namespace std;
/**
 * Assumes the graph is bi-directional.
 */
class GraphMaker
{
protected:
	vector< vector<int> > edges; //the lower indexed-node is saved in the outer vector. example: edge (5, 1) is saved in edges[1]
	int numVertex; 
public:
    GraphMaker();
	
	//A B C D E F G H I J K L M N O P Q R S T U V W X Y Z
	bool addEdgeInOrder(const int v1, const int v2);

    const vector< vector<int>  > & getEdges() const;
    void makeCheckerboard(const int row, const int col);
	void makeCircleWithCenter(const int size, const int offset);
	void makeCircleGraph(const int size);
	void makeLinearGraph(const int size);
	void makeRandomConnectedGraph(const int size, const int edgeCount);
	void makeRandomSpanningTree();
	
	void printEdges() const;
	
	


};//GraphMaker



#endif /* GRAPH_MAKER_H_ */

