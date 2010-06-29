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
#include <cstdlib>
#include <ctime>
#include <NTL/vec_ZZ.h>
#include <NTL/ZZ.h>

using namespace std;

BarvinokParameters parameters;
ReadPolyhedronData read_polyhedron_data;
string output_filename;

//dots 2 ve_ZZ's
ZZ dot(vec_ZZ& a, vec_ZZ& b){
	ZZ ans = ZZ();
	assert(a.length() == b.length());
	for(int i = 0; i < a.length(); i++){
		ans += a[i] * b[i];
	}
	return ans;
}

//overwrites the num and denom of the first number (the first 2 variables passed in)
//wih the value of adding it to the second number (the next 2 slots)
void add(ZZ& num1, ZZ& denom1, ZZ& num2, ZZ& denom2){
	assert(denom1 != 0 && denom2 != 0);
	num2 *= denom1;
	num1 *= denom2;
	denom1 *= denom2;
	num1 += num2;
	ZZ tempNum = GCD(num1, denom1);
	num1 /= tempNum;
	denom1 /= tempNum;
}

/* ----------------------------------------------------------------- */

int main(int argc, char *argv[])
{
	set_program_name(argv[0]);
	srand((unsigned)time(0));
	vec_ZZ c = vec_ZZ();

	int i;
	for (i = 1; i < argc; i++)
	{
		if (read_polyhedron_data.parse_option(argv[i]))
		{
		} else if (strncmp(argv[i], "--output-cones=", 15) == 0)
		{
			output_filename = argv[i] + 15;
		} else
		{
			cerr << "Unknown argument: " << argv[i] << endl;
			exit(1);
		}
	}
	Polyhedron *poly = read_polyhedron_data.read_polyhedron(&parameters);
	//cout << "Finished reading poly. data." << endl;

	parameters.Number_of_Variables = poly->numOfVars;

	//cout << "poly.numofvars=" << poly->numOfVars << ", degree="
	//		<< read_polyhedron_data.degree << endl;

	listCone * triangulatedCones;
	c.SetLength(poly->numOfVars);
	for(int i = 0; i < poly->numOfVars; i++)
		c[i] = rand() % 10000;

	ZZ valutation = ZZ();
	ZZ num = ZZ();
	ZZ denom = ZZ();
	denom = 1;
	ZZ tempNum = ZZ();
	ZZ tempDenom = ZZ();
	vec_ZZ vert = vec_ZZ();




	//for (listCone *currentCone = poly->cones; currentCone; currentCone
	//		= currentCone->rest)
	//{
		triangulatedCones =	decomposeCones(poly->cones, false, parameters);

		//triangulatedCones = triangulateCone(currentCone, poly->numOfVars,
		//		&parameters);

		cout << "c: ";
		for(int i = 0; i < poly->numOfVars; i++){
			cout << c[i] << " ";
		}
		cout << endl;

		for (listCone * simplicialCone = triangulatedCones; simplicialCone; simplicialCone
				= simplicialCone->rest)
		{
			//printConeToFile(cout, simplicialCone , poly->numOfVars);

			tempDenom = 1;
			vert = simplicialCone->vertex->vertex->numerators();

			//cout << vert[i] << endl;
			//vert = constructRay(simplicialCone->vertex->vertex, new rationalVector(poly->numOfVars), poly->numOfVars);
			for(int i = 0; i < poly->numOfVars; i++){
				cout << vert[i] << " ";
			}
			cout << endl;

			tempNum = dot(vert, c);
			cout << "eta: " << tempNum << endl;
			tempNum = power(tempNum, poly->numOfVars);

			//tempDenom *= abs(determinant(createConeDecMatrix(simplicialCone, poly->numOfVars, poly->numOfVars)));

			for(listVector * currentRay = simplicialCone->rays; currentRay; currentRay = currentRay->rest){
				for(int i = 0; i < poly->numOfVars; i++){
					cout << currentRay->first[i] << " ";
				}
				cout << endl;
				tempDenom *= -1 * dot(currentRay->first, c);
				cout << "xi: " << tempDenom << endl;
			}//for every ray

			num = simplicialCone->coefficient * num;
			add(num, denom, tempNum, tempDenom);
		}//for every simple cone.
	//}//for every cone
	for(int i = 2; i <= poly->numOfVars; i++){
		denom *= i;
	}
	tempNum = GCD(num, denom);
	num /= tempNum;
	denom /= tempNum;

	cout << "VOLUME: "<< endl << -num << endl << endl;
	if(denom != 1)
		cout << "/" << endl << endl << denom;
	cout << endl;

	//printListConeToFile((output_filename + ".triangulated").c_str(), triangulatedCones, 4); //poly->numOfVars);

	return 0;
}

