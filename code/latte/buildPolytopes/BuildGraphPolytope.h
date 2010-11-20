/*
 *
 *
 *  Created on: June 17, 2010
 *      Author: bedutra
 */

#ifndef BUILD_GRAPH_POLYTOPE_H_
#define BUILD_GRAPH_POLYTOPE_H_


#include "BuildRandomPolytope.h"



using namespace std;

class BuildGraphPolytope: public BuildRandomPolytope
{
private:
	vector< vector<int> > points; //list of points to add.

public:
	BuildGraphPolytope();

	enum PolytopeType { EDGE, SYMMETRIC_EDGE}; //encoding types.


	//A B C D E F G H I J K L M N O P Q R S T U V W X Y Z

	void buildPolymakeFile(const vector< vector<int> > &edges, PolytopeType polytopeType);

	void findEdgePolytope(const vector< vector<int> > &edges);
	void findSymmetricEdgePolytope(const vector< vector<int> > &edges);
};//BuildRandomPolytope



#endif /* BUILD_GRAPH_POLYTOPE_H_ */

