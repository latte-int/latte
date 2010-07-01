/*
 * Driver.cpp
 *
 *  Created on: Jun 24, 2010
 *      Author: bedutra
 */

#include <cstdlib>
#include <iostream>
#include <string>
#include <cassert>

#include "barvinok/barvinok.h"
#include "ReadPolyhedron.h"
#include "triangulation/triangulate.h"
#include "convert.h"
#include "print.h"
#include "gnulib/progname.h"
#include "barvinok/dec.h"
#include "valuation/PolytopeValuation.h"
#include "barvinok/dec.h"
#include "rational.h"
#include <cstdlib>
#include <ctime>
#include "cone.h"
#include <NTL/vec_ZZ.h>
#include <NTL/ZZ.h>
#include <NTL/mat_ZZ.h>
#include "testEhrhart/BuildRandomPolytope.h"


using namespace std;


BarvinokParameters parameters;
ReadPolyhedronData read_polyhedron_data;
string output_filename;








listCone * polytopeToCone(Polyhedron * poly)
{

	int numberOfVertices = 0;
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
	//listVector * newRay = new listVector;
	//listVector *masterList = newRay;
	listVector * masterList = new listVector;

	//masterList->rest = NULL;

	//cout << "BUILDING THE RAYS" << endl;
	for (listCone * currentCone = poly->cones; currentCone; currentCone
			= currentCone->rest)
	{

		//newRay->rest = 0;
		vec_ZZ buildRay;
		ZZ nume, denom;
		buildRay.SetLength(poly->numOfVars + 1);
		buildRay[0] = 1;

		++numberOfVertices;
		for (int i = 0; i < poly->numOfVars; ++i)
		{

			currentCone->vertex->vertex->getEntry(i, nume, denom);
			if (denom != 1 && denom != 0)
			{
				cerr
						<< "Converting from rational vector to vector failed because 1 != "
						<< denom << endl;
				exit(1);
			}//if error.

			buildRay[i + 1] = nume;
		}//for i

		//cout << buildRay << endl;
		//listVector *appendVectorToListVector(const vec_ZZ &, listVector*);
		masterList->first = buildRay;

		masterList = appendVectorToListVector(buildRay, masterList);
		//newRay->first = buildRay;

		//newRay->rest = masterList->rest;


		//masterList = newRay;

		//buildRay.
		//(newray->first).
	}//for currentCone
	//cout << "END  BUILDING THE RAYS" << endl;

	oneCone->rest = 0;
	oneCone->rays = masterList->rest;

	return oneCone;

}

ZZ factorial(const int n)
{
	ZZ product;
	product = 1;
	for (int i = n; i > 1; --i)
		product *= i;
	return product;
}

RR findDetermiantForVolume(listCone *cone, const int numOfVars)
{
	int i, numOfRays;
	mat_ZZ mat;
	listVector *tmp;

	numOfRays = lengthListVector(cone->rays);

	mat.SetDims(numOfRays - 1, numOfVars - 1);

	tmp = cone->rays;
	vec_ZZ startingRay = tmp->first;
	for (i = 1; i < numOfRays; i++)
	{
		tmp = tmp->rest;
		vec_ZZ endingRay = tmp->first - startingRay;

		for (int j = 1; j < numOfVars; ++j)
			mat[i - 1][j - 1] = endingRay[j];
	}

	//cout << "\nPRINTING MATRIX";
	//for (int i = 0; i < numOfRays - 1; i++)
	//{
	//	cout << '\n';
	//	for (int k = 0; k < numOfVars - 1; ++k)
	//		cout << mat[i][k] << " ";
	//}//for i.
	//cout << '\n';

	return to_RR(determinant(mat)) / to_RR(factorial(numOfRays - 1));

}

/**
 *
 */
RR findVolume(listCone * triangulatedPolytope, const int numOfVars)
{
	RR sum;

	sum = 0;

	for (listCone * ptr = triangulatedPolytope; ptr; ptr = ptr->rest)
	{
		sum += abs(findDetermiantForVolume(ptr, numOfVars));
	}

	return sum;

}



/* -------------------------------------------------------------- */
//dots 2 ve_ZZ's
ZZ dot(vec_ZZ& a, vec_ZZ& b) {
	ZZ ans = ZZ();
	assert(a.length() == b.length());
	for (int i = 0; i < a.length(); i++) {
		ans += a[i] * b[i];
	}
	return ans;
}

//overwrites the num and denom of the first number (the first 2 variables passed in)
//wih the value of adding it to the second number (the next 2 slots)
void add(ZZ& num1, ZZ& denom1, ZZ& num2, ZZ& denom2) {
	assert(denom1 != 0 && denom2 != 0);
	num2 *= denom1;
	num1 *= denom2;
	denom1 *= denom2;
	num1 += num2;
	ZZ tempNum = GCD(num1, denom1);
	num1 /= tempNum;
	denom1 /= tempNum;
}

void computeUniVolume(Polyhedron * poly) {
	srand((unsigned) time(0));
	vec_ZZ c = vec_ZZ();
	listCone * triangulatedCones;
	c.SetLength(poly->numOfVars);
	for (int i = 0; i < poly->numOfVars; i++)
		c[i] = rand() % 10000;
	ZZ valutation = ZZ();
	ZZ num = ZZ();
	ZZ denom = ZZ();
	denom = 1;
	ZZ tempNum = ZZ();
	ZZ tempDenom = ZZ();
	vec_ZZ vert = vec_ZZ();
	triangulatedCones = decomposeCones(poly->cones, false, parameters);
	cout << "c: ";
	for (int i = 0; i < poly->numOfVars; i++) {
		cout << c[i] << " ";
	}
	cout << endl;
	for (listCone * simplicialCone = triangulatedCones; simplicialCone; simplicialCone
			= simplicialCone->rest) {
		tempDenom = 1;
		vert = simplicialCone->vertex->vertex->numerators();
		for (int i = 0; i < poly->numOfVars; i++) {
			cout << vert[i] << " ";
		}
		cout << endl;
		tempNum = dot(vert, c);
		cout << "eta: " << tempNum << endl;
		tempNum = power(tempNum, poly->numOfVars);
		for (listVector * currentRay = simplicialCone->rays; currentRay; currentRay
				= currentRay->rest) {
			for (int i = 0; i < poly->numOfVars; i++) {
				cout << currentRay->first[i] << " ";
			}
			cout << endl;
			tempDenom *= -1 * dot(currentRay->first, c);
			cout << "xi: " << tempDenom << endl;
		}//for every ray

		tempNum = simplicialCone->coefficient * tempNum;
		add(num, denom, tempNum, tempDenom);
	}//for every simple cone.
	for (int i = 2; i <= poly->numOfVars; i++) {
		denom *= i;
	}
	tempNum = GCD(num, denom);
	num /= tempNum;
	denom /= tempNum;

	cout << "VOLUME: " << endl << num << endl << endl;
	if (denom != 1)
		cout << "/" << endl << endl << denom;
	cout << endl;
}

void computeTriangVolume(Polyhedron * poly) {
	srand((unsigned) time(0));
	vec_ZZ c = vec_ZZ();
	ZZ scale = ZZ();
	ZZ num = ZZ();
	ZZ denom = ZZ();
	denom = 1;
	ZZ tempNum = ZZ();
	ZZ tempDenom = ZZ();
	vec_ZZ vert = vec_ZZ();
	vec_ZZ ans = vec_ZZ();
	mat_ZZ mat;
	ZZ det = ZZ();
	mat.SetDims(poly->numOfVars, poly->numOfVars);
	listCone *triangulatedCones;

	c.SetLength(poly->numOfVars);
	for (int i = 0; i < poly->numOfVars; i++)
		c[i] = rand() % 10000;

	for (listCone *currentCone = poly->cones; currentCone; currentCone
			= currentCone->rest) {
		triangulatedCones = triangulateCone(currentCone, poly->numOfVars,
				&parameters);

		for (listCone * simplicialCone = triangulatedCones; simplicialCone; simplicialCone
				= simplicialCone->rest) {

			//find vertex
			vert = scaleRationalVectorToInteger(simplicialCone->vertex->vertex,
					poly->numOfVars, tempDenom);

			//raise f(vertex) to the power of the dimension
			tempNum = dot(vert, c);
			tempNum = power(tempNum, poly->numOfVars);
			tempDenom = power(tempDenom, poly->numOfVars);

			int col = 0;

			for (listVector * currentRay = simplicialCone->rays; currentRay; currentRay
					= currentRay->rest, col++) {

				//generate matrix
				for (int row = 0; row < poly->numOfVars; row++) {
					mat[row][col] = currentRay->first[row];
				}//for every component of the ray

			}//for every ray in the simple cone

			//solve for gammas, det of rays
			solve(det, ans, mat, c);

			//divide by the absolute value of the determinant
			tempDenom *= abs(det);

			//divide by the multiplication of the gammas
			tempNum *= power(det, poly->numOfVars);
			for (int i = 0; i < poly->numOfVars; i++) {
				tempDenom *= ans[i];
			}
			cout << "adding " << tempNum << " / " << tempDenom << " to " << num
					<< " / " << denom << endl;
			add(num, denom, tempNum, tempDenom);
		}//for every simple cone in the cone

	}//for every cone
	for (int i = 2; i <= poly->numOfVars; i++) {
		denom *= i;
	}
	tempNum = GCD(num, denom);
	num /= tempNum;
	denom /= tempNum;

	cout << "VOLUME: " << endl << num << endl;
	if (denom != 1)
		cout << "/" << endl << denom;
	cout << endl << to_RR(num) / to_RR(denom) << endl;
}

void printRationalFunction(Polyhedron * poly) {
	listCone * triangulatedCones;
	vec_ZZ vert = vec_ZZ();
	ZZ temp = ZZ();
	triangulatedCones = decomposeCones(poly->cones, false, parameters);
	cout << "( ";
	for (listCone * simplicialCone = triangulatedCones; simplicialCone; simplicialCone
			= simplicialCone->rest) {
		printConeToFile(cout, simplicialCone, poly->numOfVars);
		vert = scaleRationalVectorToInteger(simplicialCone->vertex->vertex,
							poly->numOfVars, temp);
		cout << "( ";
		for(int i = 0; i < poly->numOfVars - 1; i++){
			cout << vert[i];
			if(temp != 1)
				cout << " / " << temp;
			cout << " * c" << i << " + ";
		}
		cout <<  vert[poly->numOfVars - 1];
		if(temp != 1)
						cout << " / " << temp;
		cout << " * c" << poly->numOfVars - 1 << " ) ^ " << poly->numOfVars << " / ( ";

		if(poly->numOfVars % 2 == 1)
			cout << "-";
		for (listVector * currentRay = simplicialCone->rays; currentRay; currentRay
				= currentRay->rest) {
			cout << "( ";
			for(int i = 0; i < poly->numOfVars; i++){
				cout << currentRay->first[i] << " * c" << i;
				if(i != poly->numOfVars - 1){
						cout << " + ";
				}
			}
			cout << " )";
			if(currentRay->rest != NULL)
				cout << " * ";
		}//for every ray
		cout << " ) * ";
		cout << simplicialCone->coefficient;
		if(simplicialCone->rest != NULL)
			cout << " + ";
	}//for every simple cone.
	cout << ") / ( " << poly->numOfVars << "!";
	cout << " )" << endl;
}
/* ----------------------------------------------------------------- */

int main(int argc, char *argv[]) {
	BuildRandomPolytope *buildPolytope = 0;
	if ( argc > 2 )
	{
		buildPolytope = new BuildRandomPolytope(atoi(argv[1]));
		buildPolytope->setComments("Making random integer polytope for volume testing");
		buildPolytope->buildPolymakeFile(atoi(argv[2]));	//make the file
		buildPolytope->callPolymake();						//run polymake
		buildPolytope->findVolumeWithPolymake();			//run polymake for the volume
		buildPolytope->convertFacetEquations();				//fix facet equations
		buildPolytope->printFacetEquationsForLattE();		//make latte file.

	    read_polyhedron_data.filename = buildPolytope->getLatteFile();
	    read_polyhedron_data.expect_filename = false;

	}
	else
		set_program_name(argv[0]);

	int i;
	for (i = 1; i < argc; i++)
	{
		if (read_polyhedron_data.parse_option(argv[i]))
		{
		}
	}//for i.

	Polyhedron *poly = read_polyhedron_data.read_polyhedron(&parameters);
	//cout << "Finished reading poly. data." << endl;
	if (buildPolytope)
		delete buildPolytope;

	parameters.Number_of_Variables = poly->numOfVars;
#if 1

	//cout << "PRINT VERTICES" << endl;
	//for (listCone * c = poly->cones; c; c = c->rest)
	//	printRationalVectorToFile(cout, c->vertex->vertex, poly->numOfVars);
	//cout << "END PRINT VERTICES" << endl;


	PolytopeValuation polytopeValuation(poly, &parameters);
	RationalNTL rVolume =  polytopeValuation.findVolume();
	cout << "TOTAL VOLUME (from the class)(rational) = " << rVolume << " = " << rVolume.to_RR() << endl;
	cout << "TOTAL VOLUME (from the class)(float) = " << polytopeValuation.findVolume_old() << endl;
exit(0);

	listCone * polytopeAsOneCone = polytopeToCone(poly);
/*
	cout << "TESTING ENDING PTR" << endl;

	for (listVector * c = polytopeAsOneCone->rays; c; c = c->rest)
	{
		cout << c->first << endl;
	}
	cout << "END TESTING ENDING PTR" << endl;

	cout << "PRINT THE FINAL POLY" << endl;
	printConeToFile(cout, polytopeAsOneCone, poly->numOfVars + 1);
*/
	parameters.Number_of_Variables = poly->numOfVars + 1;
	listCone *triangulatedPolytope = triangulateCone(polytopeAsOneCone,
			poly->numOfVars + 1, &parameters);
/*
	cout << "PRINTING TRIANGULATED POLYTOPE" << endl;
	for (listCone * c = triangulatedPolytope; c; c = c->rest)
	printConeToFile(cout, c, poly->numOfVars + 1);
*/
	RR ans = findVolume(triangulatedPolytope, poly->numOfVars + 1);
	int asdf = 1 + 3;
	cout << "trangulated volume = " << ans << endl;

	freeListCone(triangulatedPolytope);

	exit(1);
#endif



	//computeUniVolume(poly);
	//computeTriangVolume(poly);
	printRationalFunction(poly);
	return 0;
}

