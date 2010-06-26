/*
 *
 *
 *  Created on: June 15, 2010
 *      Author: bedutra
 */

#ifndef BUILD_RANDOM_POLYTOPE_H_
#define BUILD_RANDOM_POLYTOPE_H_

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

class BuildRandomPolytope
{
protected:
    int ambientDim; 		//dim of space the poly lives in.
    string fileName; 		//file name for polymake
    string latteFile;		//file name for latte
    string PolytopeComments;//describes the current polytope (example: "edge poly of a line graph of 30 nodes")
    int maxInteger;			//when making random points, take elements in  |x| in [0, maxInteger)
    double probNegative;	//when making random points, probNegative percent of them should be negative.
    int dim;				//dim of the polytope found by polymake
    vector< vector<mpq_class> > facets; //facets of the polytope found by polymake
    int numAffineHull;		//number of affine hull facets found by polymake (saved at the end of the facets vector)
public:
    BuildRandomPolytope(int ambient_dim);
	
	//A B C D E F G H I J K L M N O P Q R S T U V W X Y Z

	void buildPolymakeFile(const int numPoints);
	void callPolymake();
	void convertFacetEquations();
	void findEhrhardPolynomial();
	void findFacetEquations();


	void printFacetEquationsForLattE();

	void setComments(const string& newComments);
};//BuildRandomPolytope



#endif /* BUILD_RANDOM_POLYTOPE_H_ */
