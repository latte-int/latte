/*
 *
 *
 *  Created on: June 15, 2010
 *      Author: bedutra
 *
 */

#ifndef BUILD_HYPERSIMPLEX_EDGE_POLYTOPE_H_
#define BUILD_HYPERSIMPLEX_EDGE_POLYTOPE_H_

#include "BuildPolytope.h"



using namespace std;

class BuildHypersimplexEdgePolytope: public BuildPolytope
{
private:
	int _numOnes;


public:
	BuildHypersimplexEdgePolytope();
	
	//A B C D E F G H I J K L M N O P Q R S T U V W X Y Z

	void addToPoints(vector< vector<mpq_class> > &list, vector<mpq_class> currentPoint, int base, bool addCurrent = false);
	void generatePoints(int ambient_dim, int numones);

	
};//BuildRandomEdgePolytope



#endif /* BUILD_HYPERSIMPLEX_EDGE_POLYTOPE_H_ */

