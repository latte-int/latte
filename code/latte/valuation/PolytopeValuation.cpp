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
	vertexRayCones(p->cones), parameters(bp), polytopeAsOneCone(NULL), triangulatedPoly(NULL)
{
	numOfVars = p->numOfVars;
}


PolytopeValuation::PolytopeValuation(listCone *vertexraycones, int numofvars, BarvinokParameters *bp) :
	vertexRayCones(vertexraycones), parameters(bp), polytopeAsOneCone(NULL), triangulatedPoly(NULL), numOfVars(numofvars)
{
}


PolytopeValuation::~PolytopeValuation()
{
	//don't free parameters or vertexRayCones, because we did not make them!

	if (polytopeAsOneCone) freeListCone(polytopeAsOneCone);
	if (triangulatedPoly) freeListCone(triangulatedPoly);
}

void PolytopeValuation::convertToOneCone()
{
	if ( ! vertexRayCones )
	{
		cout << "vertexRayCones* is not defined" << endl;
		exit(1);
	}//error.
	if ( polytopeAsOneCone)
		return ; //already did computation.

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
	oneCone->vertex->vertex = new rationalVector(numOfVars + 1);

	oneCone->rays = new listVector;
	oneCone->rays->rest = 0;

	//now add the vertex to the rays list with a leading 1: (1, old poly cone vertex).
	//The first entry in masterList should be ignored because masterList->first = masterList->rest->first.
	listVector * masterList = new listVector;



	for (listCone * currentCone = vertexRayCones; currentCone; currentCone
			= currentCone->rest)
	{
		vec_ZZ buildRay; //buildRay = [1, old-vertex]
		ZZ nume, denom;
		buildRay.SetLength(numOfVars + 1);


		ZZ scaleFactor; //not used, but need to pass into scaleRationalVectorToInteger().
		vec_ZZ integerVertex = scaleRationalVectorToInteger(currentCone->vertex->vertex, numOfVars, scaleFactor);

		buildRay[0] = scaleFactor; //1 * scaleFactor.
		for (int i = 0; i < numOfVars; ++i)
		{

			//currentCone->vertex->vertex->getEntry(i, nume, denom); //get vertex[i] = nume/denome.
			//if (denom != 1 || denom == 0)
			//{
			//	cerr
			//			<< "Converting from rational vector to vector failed because denom. = "
			//			<< denom << endl;
			//	exit(1);
			//}//if error.

			//buildRay[i + 1] = nume;
			buildRay[i + 1] = integerVertex[i];
		}//for i

		//cout << buildRay << endl;

		masterList->first = buildRay;
		masterList = appendVectorToListVector(buildRay, masterList);
	}//for currentCone

		//cout << "END  BUILDING THE RAYS" << endl;

	oneCone->rest = 0;
	oneCone->rays = masterList->rest; //ignore masterList->first, so just return the rest and NOT masterList.


	polytopeAsOneCone = oneCone; //finally, save the result.
}//convertToOneCone


RationalNTL PolytopeValuation::findDetermiantForVolume(const listCone * oneSimplex) const
{
	int i, numOfRays;
	mat_ZZ mat;


	vec_ZZ head;
	vec_ZZ tail;
	ZZ numerator, denominator;
	numerator = 1;
	denominator = 1;


	head.SetLength(numOfVars + 1);
	tail.SetLength(numOfVars + 1);


	numOfRays = lengthListVector(oneSimplex->rays);

	mat.SetDims(numOfRays, numOfVars+1);


	i = 0;
	for(listVector * currentRay = oneSimplex->rays; currentRay; currentRay = currentRay->rest)
	{
		for (int k = 0; k < numOfVars +1; ++k)
			mat[i][k] = ((currentRay->first)[k]);
		denominator *= (currentRay->first)[0];
		++i;
	}//for currentRay

	//cout << "numerator = " << numerator << endl;
	numerator =  abs(determinant(mat));
	denominator *= factorial(numOfRays -1);
	//cout << mat << " = " << determinant(mat) << "\n./." << factorial(numOfRays -1) << endl;
	return RationalNTL(numerator, denominator);
}//findDetermiantForVolume


RR PolytopeValuation::findDetermiantForVolume_old(const listCone * oneSimplex) const
{
	int i, numOfRays;
	mat_RR mat;
	listVector *tmp;

	vec_RR head;
	vec_RR tail;

	head.SetLength(numOfVars + 1);
	tail.SetLength(numOfVars + 1);


	numOfRays = lengthListVector(oneSimplex->rays);

	mat.SetDims(numOfRays - 1, numOfVars);


	tmp = oneSimplex->rays;
	vec_ZZ startingRay = tmp->first;
	for (i = 1; i < numOfRays; i++)
	{
		tmp = tmp->rest;

		for(int k = 0; k < (tmp->first).length(); ++k)
		{
			head[k] = to_RR((tmp->first)[k]);
			tail[k] = to_RR(startingRay[k]);
		}//convert from ZZ to RR vectors.

		//RR inverse1, inverse2;
		//inverse1 = 1 / to_RR((tmp->first)[0]);
		//inverse2 = 1 / to_RR(startingRay[0]);

		//If head and tail are two vectors, head - tail is a vector from tail to head
		// so endingRay is a vector from the starting vertex (corresponding to startingRay) to all the other vertices.
		head = (head)*(1 / to_RR((tmp->first)[0]));
		tail = (tail)*(1 / to_RR(startingRay[0])); //must convert to vect_RR by dividing by leading element (to get the form [1 old vertex])

		vec_RR endingRay =  head - tail;


		for(int j = 1; j < numOfVars + 1 ; ++j)
			mat[i - 1][j - 1] = endingRay[j];
	}//for i.

	//cout << mat << " = " << determinant(mat) << "\n./." << factorial(numOfRays -1) << endl;
	return abs((determinant(mat))/to_RR((factorial(numOfRays -1))));
}//findDetermiantForVolume

RationalNTL PolytopeValuation::findVolume()
{
	RationalNTL sum;

	convertToOneCone();
	triangulatePolytopeCone();

	int i = 0;
	for (listCone * oneSimplex = triangulatedPoly; oneSimplex; oneSimplex = oneSimplex->rest)
	{
		sum.add(findDetermiantForVolume(oneSimplex));
		//cout << i++ << " sub-volume " << findDetermiantForVolume(oneSimplex) << endl << endl;
	}

	return sum;

}//findVolume

RR PolytopeValuation::findVolume_old()
{
	RR sum;

	convertToOneCone();
	triangulatePolytopeCone();

	int i = 0;
	for (listCone * oneSimplex = triangulatedPoly; oneSimplex; oneSimplex = oneSimplex->rest)
	{
		sum += (findDetermiantForVolume_old(oneSimplex));
		//cout << i++ << " sub-volume_old " << findDetermiantForVolume_old(oneSimplex) << endl << endl;
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
	if ( triangulatedPoly)
		return ; //allready did computation.
	if ( polytopeAsOneCone == NULL)
	{
		cout << "PolytopeValuation::triangulatePolytopeCone(): there is no code to triangulate" << endl;
		exit(1);
	}


	parameters->Number_of_Variables = numOfVars + 1;
	triangulatedPoly = triangulateCone(polytopeAsOneCone, numOfVars + 1, parameters);
}//triangulateCone()
