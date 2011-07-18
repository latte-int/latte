/*
 * count.h
 *
 *  Created on: Jul 16, 2011
 *      Author: bedutra
 */

#ifndef COUNT_H_
#define COUNT_H_

#include <string.h>
#include <stdio.h>
#include <cassert>

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
#include "gnulib/pathmax.h"

class CountAnswerContainer
{
public:
	vec_ZZ seriesExpansion;
	ZZ numLaticePoints;
	string multivariateGenFunctionFileName;
	mpq_vector ehrhart_coefficients;
	void checkPolynomial();
};

CountAnswerContainer mainCountDriver(int argc, char *argv[]);

#endif /* COUNT_H_ */
