/*
 * PolytopeValuation.cpp
 *
 *  Created on: Jun 25, 2010
 *      Author: bedutra
 */

#include "PolytopeValuation.h"

using namespace std;

/**
 * Does not keep a local copy of the polyhedron, only a pointer.
 * Polyhedron *poly contains the list of vertex-ray pairs.
 */
PolytopeValuation::PolytopeValuation(Polyhedron *p, BarvinokParameters *bp) :
	poly(p), parameters(bp), polytopeAsOneCone(NULL), triangulatedPoly(NULL)
{
	// TODO Auto-generated constructor stub

}

PolytopeValuation::~PolytopeValuation()
{
	//don't free parameters or poly, because we did not make them!

	if (polytopeAsOneCone) freeListCone(polytopeAsOneCone);
	if (triangulatedPoly) freeListCone(triangulatedPoly);
}

void PolytopeValuation::convertToOneCone()
{
	if ( ! poly )
	{
		cout << "Polyhedron* is not defined" << endl;
		exit(1);
	}//error.

	listCone * oneCone = new listCone();
	oneCone->coefficient = 1;
	oneCone->determinant = 0;
	oneCone->subspace_generators = NULL;
	oneCone->dual_determinant = 0;
	oneCone->facets = NULL;
	oneCone->equalities = NULL;
	oneCone->latticePoints = NULL;
	oneCone->rest = NULL;

	//set to zero vector of numofvars + 1 size.
	oneCone->vertex = new Vertex();
	oneCone->vertex->vertex = new rationalVector(poly->numOfVars + 1);

	oneCone->rays = new listVector;
	oneCone->rays->rest = 0;

	//now add the vertex to the rays list with a leading 1: (1, old poly cone vertex).
	//The first entry in masterList should be ignored because masterList->first = masterList->rest->first.
	listVector * masterList = new listVector;


	//FIXME: rays are integer, but vertex are rational. Will need to chagne this.
	//cout << "BUILDING THE RAYS" << endl;
	for (listCone * currentCone = poly->cones; currentCone; currentCone
			= currentCone->rest)
	{
		vec_ZZ buildRay; //buildRay = [1, old-vertex]
		ZZ nume, denom;
		buildRay.SetLength(poly->numOfVars + 1);
		buildRay[0] = 1;


		for (int i = 0; i < poly->numOfVars; ++i)
		{

			currentCone->vertex->vertex->getEntry(i, nume, denom); //get vertex[i] = nume/denome.
			if (denom != 1 || denom == 0)
			{
				cerr
						<< "Converting from rational vector to vector failed because denom. = "
						<< denom << endl;
				exit(1);
			}//if error.

			buildRay[i + 1] = nume;
		}//for i

		//cout << buildRay << endl;

		masterList->first = buildRay;
		masterList = appendVectorToListVector(buildRay, masterList);
	}//for currentCone

	//	cout << "END  BUILDING THE RAYS" << endl;

	oneCone->rest = 0;
	oneCone->rays = masterList->rest; //ignore masterList->first, so just return the rest and NOT masterList.

	if ( polytopeAsOneCone) delete polytopeAsOneCone;
	polytopeAsOneCone = oneCone; //finally, save the result.
}//convertToOneCone




RR PolytopeValuation::findDetermiantForVolume(const listCone * oneSimplex) const
{
	int i, numOfRays;
	mat_ZZ mat;
	listVector *tmp;

	numOfRays = lengthListVector(oneSimplex->rays);

	mat.SetDims(numOfRays - 1, poly->numOfVars);

	tmp = oneSimplex->rays;
	vec_ZZ startingRay = tmp->first;
	for (i = 1; i < numOfRays; i++)
	{
		tmp = tmp->rest;
		vec_ZZ endingRay = tmp->first - startingRay; //FIXME: must convert to vect_RR by dividing by leading element (to get the form [1 old vertex])

		for(int j = 1; j < poly->numOfVars + 1 ; ++j)
			mat[i - 1][j - 1] = endingRay[j];
	}

	return to_RR(determinant(mat))/to_RR(factorial(numOfRays -1));
}//findDetermiantForVolume

RR PolytopeValuation::findVolume()
{
	RR sum;
	sum = 0;

	convertToOneCone();
	triangulatePolytopeCone();

	for (listCone * oneSimplex = triangulatedPoly; oneSimplex; oneSimplex = oneSimplex->rest)
	{
		sum += abs(findDetermiantForVolume(oneSimplex));
	}

	return sum;

}//findVolume



ZZ  PolytopeValuation::factorial(const int n)
{
	ZZ product;
	product = 1;
	for(int i = n; i > 1; --i)
		product *= i;
	return product;
}//factorial


void PolytopeValuation::triangulatePolytopeCone()
{
	parameters->Number_of_Variables = poly->numOfVars + 1;
	triangulatedPoly = triangulateCone(polytopeAsOneCone, poly->numOfVars + 1, parameters);
}//triangulateCone()
