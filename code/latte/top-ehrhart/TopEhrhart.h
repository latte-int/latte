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
public:
	TopEhrhart(Polyhedron * polyhedron, BarvinokParameters & para, int numTopCoeff, bool real);


	void computeTopEhrhartPolynomial(const monomialSum & polynomial);
	void computeTopEhrhartPolynomial(const linFormSum  & linForm);
	virtual ~TopEhrhart();
};

#endif /* TOPEHRHART_H_ */
