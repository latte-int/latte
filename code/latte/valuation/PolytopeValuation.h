/*
 * PolytopeValuation.h
 *
 *  Created on: Jun 25, 2010
 *      Author: bedutra
 */

#ifndef POLYTOPEVALUATION_H_
#define POLYTOPEVALUATION_H_

#include <cstdlib>
#include <iostream>
#include <string>

#include "barvinok/barvinok.h"
#include "ReadPolyhedron.h"
#include "triangulation/triangulate.h"
#include "convert.h"
#include "print.h"
#include "gnulib/progname.h"
#include "barvinok/dec.h"
#include "valuation/PolytopeValuation.h"
#include <NTL/vec_ZZ.h>
#include "rational.h"
#include "cone.h"

using namespace std;




class PolytopeValuation
{

	Polyhedron *poly;				//lists of vertex-ray pairs.
	BarvinokParameters *parameters; //Barvinok Parameters.
	listCone * polytopeAsOneCone;	//From poly, create one code with vertex=[0,0...0], rays={[1, v] | v is a vertex of the polytope}
	listCone * triangulatedPoly;	//The triangulation of polytopeAsOneCone.



public:
	PolytopeValuation(Polyhedron *p, BarvinokParameters *bp);
	virtual ~PolytopeValuation();


	// A B C D E F G H I J K L M N O P Q R S T U V W X Y Z

	void convertToOneCone(); //convert from poly to polytopeAsOneCone
	RR findDetermiantForVolume(const listCone * oneSimplex) const;
	RR findVolume();		 //finds the volume of the Polyhedron.
	ZZ static factorial(const int n);
	void triangulatePolytopeCone();  //convert polytopeAsOneCone to triangulatedPoly

};

#endif /* POLYTOPEVALUATION_H_ */
