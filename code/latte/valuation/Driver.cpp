/*
 * Driver.cpp
 *
 *  Created on: Jun 24, 2010
 *      Author: bedutra
 */

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

BarvinokParameters parameters;
ReadPolyhedronData read_polyhedron_data;


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
	cout << "got here 50" << endl;
	//masterList->rest = NULL;

	cout << "BUILDING THE RAYS" << endl;
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

		cout << buildRay << endl;
		//listVector *appendVectorToListVector(const vec_ZZ &, listVector*);
		masterList->first = buildRay;

		masterList = appendVectorToListVector(buildRay, masterList);
		//newRay->first = buildRay;

		//newRay->rest = masterList->rest;


		//masterList = newRay;

		//buildRay.
		//(newray->first).
	}//for currentCone
	cout << "END  BUILDING THE RAYS" << endl;

	oneCone->rest = 0;
	oneCone->rays = masterList->rest;

	return oneCone;

}

ZZ factorial(const int n)
{
	ZZ product;
	product = 1;
	for(int i = n; i > 1; --i)
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

		for(int j = 1; j < numOfVars ; ++j)
			mat[i - 1][j - 1] = endingRay[j];
	}

	cout << "\nPRINTING MATRIX";
	for(int i = 0; i < numOfRays-1; i++)
	{
		cout << '\n';
		for(int k = 0; k < numOfVars -1; ++k)
			cout << mat[i][k] << " ";
	}//for i.
	cout << '\n';

	return to_RR(determinant(mat))/to_RR(factorial(numOfRays -1));


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

/* ----------------------------------------------------------------- */

int main(int argc, char *argv[])
{
	set_program_name(argv[0]);

	int i;
	for (i = 1; i < argc; i++)
	{
		if (read_polyhedron_data.parse_option(argv[i]))
		{
		} else
		{
			cerr << "Unknown argument: " << argv[i] << endl;
			exit(1);
		}
	}//for i.

	Polyhedron *poly = read_polyhedron_data.read_polyhedron(&parameters);
	cout << "Finished reading poly. data." << endl;

#if 1

	//cout << "PRINT VERTICES" << endl;
	//for (listCone * c = poly->cones; c; c = c->rest)
	//	printRationalVectorToFile(cout, c->vertex->vertex, poly->numOfVars);
	//cout << "END PRINT VERTICES" << endl;


	PolytopeValuation polytopeValuation(poly, &parameters);
	cout << "TOTAL VOLUME (from the class) = " << polytopeValuation.findVolume() << endl;
exit(0);


	listCone * polytopeAsOneCone = polytopeToCone(poly);

	cout << "TESTING ENDING PTR" << endl;

	for (listVector * c = polytopeAsOneCone->rays; c; c = c->rest)
	{
		cout << c->first << endl;
	}
	cout << "END TESTING ENDING PTR" << endl;

	cout << "PRINT THE FINAL POLY" << endl;
	printConeToFile(cout, polytopeAsOneCone, poly->numOfVars + 1);

	parameters.Number_of_Variables = poly->numOfVars + 1;
	listCone *triangulatedPolytope = triangulateCone(polytopeAsOneCone,
			poly->numOfVars + 1, &parameters);

	cout << "PRINTING TRIANGULATED POLYTOPE" << endl;
	for (listCone * c = triangulatedPolytope; c; c = c->rest)
		printConeToFile(cout, c, poly->numOfVars + 1);

	RR ans = findVolume(triangulatedPolytope, poly->numOfVars + 1);
	int asdf = 1 + 3;
	cout << "trangulated volume = " << ans << endl;

	freeListCone(triangulatedPolytope);

	exit(1);
#endif

#if 0
	parameters.Number_of_Variables = poly->numOfVars;

	//printListConeToFile(output_filename.c_str(), poly->cones, poly->numOfVars);
	cout << "poly.numofvars=" << poly->numOfVars << ", degree="
	<< read_polyhedron_data.degree << endl;

	listCone * triangulatedCones;
	ZZ valutation = ZZ();

	/*
	 cout << "starting decomposeCones_single call" << endl;

	 decomposeCones_Single(poly->cones, poly->numOfVars + 1,
	 read_polyhedron_data.degree, 8, argv[argc - 1], 1, false,
	 BarvinokParameters::DualDecomposition);

	 */

	/*
	 listCone* decomposedCones = decomposeCones(poly->cones, poly->numOfVars, 0,
	 argv[argc - 1], 1, true, BarvinokParameters::DualDecomposition);

	 printListConeToFile((output_filename + ".decomposed").c_str(),
	 decomposedCones, poly->numOfVars);

	 ZZ sum;
	 for (listCone * current = decomposedCones; current; current = current->rest)
	 sum += current->coefficient * abs(current->determinant);
	 cout << "SUM=" << sum << endl;
	 cout << "****************************************************" << endl;
	 */

	//decomposeCones_Single(cones, numOfVars, degree, flags, fileName, 1, false,
	//				      BarvinokParameters::DualDecomposition);
	//

	cout << "*************************\n";
	cout << "Printing the Poly's cones:" << endl;

	for(listCone *ptr = poly->cones; ptr; ptr = ptr->rest)
	printConeToFile(cout, ptr, poly->numOfVars);

	cout << "*******************************\n";
	cout << "Printing the Trangulated cones:" << endl;

	for (listCone *currentCone = poly->cones; currentCone; currentCone
			= currentCone->rest)
	{
		triangulatedCones = triangulateCone(currentCone, poly->numOfVars,
				&parameters);

		for (listCone * simplicalCone = triangulatedCones; simplicalCone; simplicalCone
				= simplicalCone->rest)
		{
			printConeToFile(cout, simplicalCone, poly->numOfVars);
			ZZ det = determinant((createConeDecMatrix(simplicalCone,
									poly->numOfVars, poly->numOfVars)));
			det = abs(det) * simplicalCone->coefficient;
			cout << "this det=" << det << endl;
			valutation += det;

		}//for every simply cone.
	}//for every cone

	cout << "TOTAL: " << valutation << endl;
#endif
	//printListConeToFile((output_filename + ".triangulated").c_str(), triangulatedCones, 4); //poly->numOfVars);

	return 0;
}

