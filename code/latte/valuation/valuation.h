/*
 * valuation.h
 *
 *  Created on: Jul 21, 2010
 *      Author: bedutra
 */

#ifndef VALUATION_H_
#define VALUATION_H_

#include <cstdlib>
#include <iostream>
#include <string>
#include <cassert>
#include <vector>
#include <string.h>
#include <stdio.h>
#include <cstdlib>
#include <ctime>
#include <NTL/vec_ZZ.h>
#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <NTL/mat_ZZ.h>

#include "valuation/PolytopeValuation.h"
#include "valuation/RecursivePolytopeValuation.h"

#include "CheckEmpty.h"
#include "Polyhedron.h"
/* START COUNT INCLUDES */
#include "config.h"
#include "latte_ntl_integer.h"
#include "barvinok/dec.h"
#include "barvinok/barvinok.h"
#include "barvinok/Triangulation.h"
#include "vertices/cdd.h"
#include "genFunction/maple.h"
#include "genFunction/piped.h"
#include "cone.h"
#include "dual.h"
#include "RudyResNTL.h"
#include "Residue.h"
#include "Grobner.h"
#include "preprocess.h"
#include "print.h"
#include "ramon.h"
#include "rational.h"
#include "timing.h"
#include "flags.h"
#include "ExponentialSubst.h"
#include "latte_random.h"
#include "Irrational.h"
#include "ExponentialEhrhart.h"
#include "triangulation/triangulate.h"
#include "genFunction/matrix_ops.h"
#ifdef HAVE_EXPERIMENTS
#include "ExponentialApprox.h"
#include "TrivialSubst.h"
#endif

#include "banner.h"
#include "convert.h"
#include "latte_system.h"
#include "Polyhedron.h"
#include "ReadPolyhedron.h"
#include "ProjectUp.h"

#include "gnulib/progname.h"

/* END COUNT INCLUDES */

using namespace std;






namespace Valuation
{
class ValuationData
{
public:
	//Notes:
	//1) volumeLawrence, volumeTriangulation, and integrateTriangulation timers
	//start from when the tangent cones are computed to the volume/integral computation.
	//2) entireValuation timer starts from when mainValuationDriver is called to when it finishes.
	//Also "answer" is meaningless for "entireValuation".

	enum ValuationType { unknown, volumeLawrence, volumeTriangulation, integrateTriangulation, integrateLawrence, entireValuation};
	ValuationType valuationType;
	RationalNTL answer;
	Timer timer;

	ValuationData();
};


class ValuationContainer
{

public:
	vector<ValuationData> answers;

	//Adds a ValuationData to the vector
	void add(const ValuationData & d );

	//Prints the stats: valuation type, valuation time, total program time, etc.
	void printResults(ostream & out) const;
};


typedef enum {inputPolynomial, inputLinearForm, nothing} IntegrandType;

ValuationContainer computeVolume(Polyhedron * poly,
		BarvinokParameters &myParameters, const char *valuationAlg,
		const char * print);

ValuationContainer computeIntegral(Polyhedron *poly,
		BarvinokParameters &myParameters, const char *valuationAlg,
		const char * polynomialString, const IntegrandType);

ValuationContainer mainValuationDriver(const char *argv[], int argc);

void polyhedronToCones(const char valuationAlg[], Polyhedron *Poly, BarvinokParameters * params);

static void usage(const char *progname);

}//namespace Valuation





#endif /* VALUATION_H_ */
