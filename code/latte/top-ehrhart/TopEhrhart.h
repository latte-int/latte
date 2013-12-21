/*
 * TopEhrhart.h
 *
 *  Created on: Jan 2, 2012
 *      Author: Koeppe, Brandon
 */

#ifndef TOPEHRHART_H_
#define TOPEHRHART_H_


#include <cstdlib>
#include <iostream>
#include <string>
#include <sstream>

#include "barvinok/barvinok.h"
#include "ReadPolyhedron.h"
#include "triangulation/triangulate.h"
#include "convert.h"
#include "print.h"
#include "gnulib/progname.h"
#include "barvinok/dec.h"
#include <NTL/vec_ZZ.h>
#include <NTL/vec_RR.h>
#include "rational.h"
#include "cone.h"

/* Integration Headers */
#include "integration/PolyTrie.h"
#include "integration/newIntegration.h"

/* for maple */
#include "latte_relocatable.h"
#include "latte_system.h"

class TopEhrhart
{
private:
	BarvinokParameters &parameters; //Barvinok Parameters.
	Polyhedron * poly;				//The polyhedron, vertexRayCones or PolytopeAsOneCone points to the polyhedron's cones.
	int numTopCoefficients;
	bool realDilations;
	string saveTopEhrhartPolynomial;
public:
	TopEhrhart(Polyhedron * polyhedron, BarvinokParameters & para, int numTopCoeff, bool real, string savePolynomial = "-1");


	void computeTopEhrhartPolynomial(const monomialSum & polynomial); //the weight function is a polynomial
	void computeTopEhrhartPolynomial(const linFormSum  & linForm);    //the weight function is a linear form sum
	void computeTopEhrhartPolynomial();								  //the weight function is 1 (no weight)
	virtual ~TopEhrhart();
};


#endif
