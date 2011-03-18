/*
 * PolytopeValuation.h
 *
 *  Created on: Jun 25, 2010
 *      Author: Brandon Dutra and Gregory Pinto
 */

#ifndef POLYTOPEVALUATION_H_
#define POLYTOPEVALUATION_H_

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
#include "valuation/PolytopeValuation.h"
#include <NTL/vec_ZZ.h>
#include <NTL/vec_RR.h>
#include "rational.h"
#include "cone.h"

/* Integration Headers */
#include "integration/PolyTrie.h"
#include "integration/newIntegration.h"


using namespace std;

/*
 * Constructing:
 * 1) Pass in a Polyhedron *.
 *
 * Finding Volume:
 * 1) Call findVolume(DeterminantVolume) if you want to use the Determinant method. Will convert first
 * 		convert  vertexRayCones into a triangulation if needed.
 * 2) Call findVolume(LawrenceVolume) if you want to use the Lawrence method. Will convert first
 * 		convert  vertexRayCones into a triangulation if needed using decomposeCones
 *
 */


class PolytopeValuation
{
private:
	BarvinokParameters &parameters; //Barvinok Parameters.
	Polyhedron * poly;				//The polyhedron, vertexRayCones or PolytopeAsOneCone points to the polyhedron's cones.
	listCone * vertexRayCones;		//list of  vertex-ray pairs.
	listCone * polytopeAsOneCone;	//From poly, create one code with vertex=[0,0...0], rays={[1, v] | v is a vertex of the polytope}
	listCone * triangulatedPoly;	//The triangulation of polytopeAsOneCone.
	int numOfVars, numOfVarsOneCone;
	bool freeVertexRayCones, freePolytopeAsOneCone, freeTriangulatedPoly; //denotes if we made these objects (and should free them) or if they were passed in.


	// A B C D E F G H I J K L M N O P Q R S T U V W X Y Z
	void convertToOneCone(); //convert from poly to polytopeAsOneCone. Then late you can triangulate the polytope into simpleces
	void dilatePolytopeOneCone(const ZZ & factor); //dilates polytope by changing the vertices.
	void dilatePolytopeVertexRays(const RationalNTL & factor); //dilates polytope by changing the vertices.
	void dilatePolynomialToLinearForms(linFormSum &linearForms, const monomialSum& originalPolynomial, const ZZ &dilationFactor); //given a polynomial and a dilation factor, replaces x_i^k --> factor^k*x_i^k
	ZZ findDilationFactorOneCone() const;
	ZZ findDilationFactorVertexRays() const;
	RationalNTL findIntegralUsingTriangulation(linFormSum &forms) const; //computes the integral over every simplex
	RationalNTL findIntegralUsingLawrence(linFormSum &forms) const;      //computes the integral over every simple cone.
	RationalNTL findVolumeUsingDeterminant(const listCone * oneSimplex) const; //computes the volume over every simplex
	RationalNTL findVolumeUsingLawrence();								 //computes the volume over every simple cone.
	void triangulatePolytopeCone();  			//convert polytopeAsOneCone to triangulatedPoly (triangulates the polytope)
	void triangulatePolytopeVertexRayCone();	//convert vertexRayCones to triangulatedPoly using decomposeCones (triangulates the vertex cones)


public:
	typedef enum {DeterminantVolume, LawrenceVolume} VolumeType;
	typedef enum {TriangulationIntegration, LawrenceIntegration} IntegrationType;
	typedef enum {VertexRayCones, TriangulatedCones} ConeType;

	PolytopeValuation(Polyhedron *p, BarvinokParameters &bp);
	virtual ~PolytopeValuation();


	// A B C D E F G H I J K L M N O P Q R S T U V W X Y Z
	RationalNTL findIntegral(const monomialSum& polynomial, const IntegrationType integrationType);
	RationalNTL findVolume(const VolumeType v);	//finds the volume of the Polyhedron.
	ZZ static factorial(const int n);			//computes n!
	ZZ static lcm(const ZZ &a, const ZZ & b);
	void printLawrenceVolumeFunction();			//Finds the Lawrence rational function for the volume. triangulates vertexRayCones if needed.

};

#endif /* POLYTOPEVALUATION_H_ */
