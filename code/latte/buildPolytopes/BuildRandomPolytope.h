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
#include <cassert>
#include "BuildPolytope.h"



using namespace std;

class BuildRandomPolytope: public BuildPolytope
{
protected:
    int _maxInteger;			//when making random points, take elements in  |x| in [0, maxInteger)
    double _probNegative;	//when making random points, probNegative percent of them should be negative.

public:
    BuildRandomPolytope();
	
	//A B C D E F G H I J K L M N O P Q R S T U V W X Y Z

	void makePoints(int ambient_dim, int numberPoints, int maxInt, double probNeg);
	void makePoints(int ambient_dim, int numberPoints);

};//BuildRandomPolytope



#endif /* BUILD_RANDOM_POLYTOPE_H_ */

