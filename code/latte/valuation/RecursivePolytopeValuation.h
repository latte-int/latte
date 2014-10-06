/*
 * RecursivePolytopeValuation.h
 *
 *  Created on: May 22, 2011
 *      Author: bedutra
 */

#ifndef RECURSIVEPOLYTOPEVALUATION_H_
#define RECURSIVEPOLYTOPEVALUATION_H_

#include "PolytopeValuation.h"

/**
 * DO NOT USE THIS CLASS: it is experimental/partly finished
 */
class RecursivePolytopeValuation {
public:
	RecursivePolytopeValuation();
	virtual ~RecursivePolytopeValuation();


	RationalNTL findVolume(ReadPolyhedronDataRecursive & readPolyhedron, BarvinokParameters * parm);
	RationalNTL findVolume_recursive(ReadPolyhedronDataRecursive & readPolyhedron, BarvinokParameters * parm, int power, vec_ZZ &l);

	void setMaxRecursiveLevel(int);
	void setMinDimension(int); //at what dimension is the integration "easy" and stokes need not  be used?
private:
	int maxRecursiveLevel;
	int recursiveLevel;
	int minAffineDimension;

};

#endif /* RECURSIVEPOLYTOPEVALUATION_H_ */
