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
#include <cstring>
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
#include "top-ehrhart/TopEhrhart.h"

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
	//1) for all volumeLawrence, volumeTriangulation, and integrateTriangulation timers
	//start from when the tangent cones are computed to the volume/integral computation.
	//2) entireValuation timer starts from when mainValuationDriver is called to when it finishes.
	//Also "answer" is meaningless for "entireValuation".
								//dummy, saves the total computational time.
	PolytopeValuation::ValuationAlgorithm valuationType;
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

class IntegrationInput
{
public:
	//These bool's

	typedef enum {inputVolume, inputPolynomial, inputLinearForm, inputProductLinearForm, nothing} IntegrandType;
	// FIXME: maybe rename inputVolume to inputUnweighted.  And what's
	// the difference to `nothing'? --mkoeppe
	// We could rename inputVolume. If it is set, it means the user typed "./exe --valuation=volume file.latte"
	// If nothing is set, then the user typed "./exe file.latte" (didn't give a valuation command).
	// Also, these enum's are prefixed with "input" because if a user gives a polynomial, we
	// do not know (yet) if they want to decompose to linear forms or think of the polynomial
	// as a product of linear forms--Brandon
	IntegrandType integrandType;
	string fileName;
	string integrand;

	//volume algorithms
	bool volumeCone;									//volume using the cone method.
	bool volumeTriangulation;							//volume using triangulation

	//integration algorithms.
	bool integratePolynomial;
	bool integratePolynomialAsLinearFormTriangulation; 	//decompose polynomial to LF, use triangulation.
	bool integratePolynomialAsLinearFormCone;			//decompose polynomial to LF, use cone method.

	bool integratePolynomialAsPLFTriangulation; 		//decompost polynomial to PLF, use triangulation.

	bool integrateLinearFormTriangulation;				//integrate linear forms using triangulation
	bool integrateLinearFormCone;						//integrate linear forms using cone method

	bool integrateProductLinearFormsTriangulation;		//integrate product of linear forms using triangulation.

	//Ehrhart algorithms.
	bool topEhrhart;					//compute top Ehrhart coefficients using cone method only.
	int numEhrhartCoefficients;
	bool realDilations;
	string saveTopEhrhartPolynomial;
	bool unweightedCounting;

	
	IntegrationInput();
	void processUserInput();	//set the above integration/volume algorithmic booleans.

	//Command line options. These are used by processUserInput()
	//  to set up which algorithms will be used.
	bool valuationVolume;		//--valuation=volume
	bool valuationIntegrate;	//etc
	bool valuationEhrhart;

	bool useTangentCones;		//--cone-decompose
	bool useTriangulation;		//--triangulate

	bool polynomialAsPLF;		//--polynomial-as-plf


};

ValuationContainer computeVolume(Polyhedron * poly,
		BarvinokParameters &myParameters, const IntegrationInput &intInput,
		const char * print);

//Does integration on specific integrands.
ValuationContainer computeIntegral(Polyhedron *poly,
		BarvinokParameters &myParameters, const IntegrationInput &intInput);
ValuationContainer computeIntegralPolynomial(Polyhedron *poly,
		BarvinokParameters &myParameters, const IntegrationInput & intInput);
ValuationContainer computeIntegralLinearForm(Polyhedron *poly,
		BarvinokParameters &myParameters, const IntegrationInput & intInput);
ValuationContainer computeIntegralProductLinearForm(Polyhedron *poly,
		BarvinokParameters &myParameters, const IntegrationInput & intInput);
void 			   computeTopEhrhart(Polyhedron *poly,
		BarvinokParameters &myParameters, const IntegrationInput & intInput); //Computes top weighted Ehrhart coefficients


ValuationContainer mainValuationDriver(const char *argv[], int argc);

void polyhedronToCones(const IntegrationInput &intInput, Polyhedron *Poly, BarvinokParameters * params);

static void usage(const char *progname);

}//namespace Valuation





#endif /* VALUATION_H_ */
