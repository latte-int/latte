/*
 * Driver
 *
 *  Created on: June 15, 2010
 *      Author: Brandon
 */
#include "BuildRandomPolytope.h"
#include "BuildRandomEdgePolytope.h"
#include <iostream>


void doRandom()
{
	BuildRandomPolytope rPoly(6);
	
	rPoly.buildPolymakeFile(3);
	//rPoly.callPolymake();
	rPoly.findEhrhardPolynomial();
}


void doEdge()
{
	BuildRandomEdgePolytope rPoly(6,3);
	
	rPoly.buildPolymakeFile();
	//rPoly.callPolymake();
	//rPoly.findEhrhardPolynomial();

}
int main()
{

	//doRandom();
	doEdge();

	
	return 0;
}//main

