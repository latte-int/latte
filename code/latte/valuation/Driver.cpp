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
#include <NTL/vec_ZZ.h>
#include <NTL/ZZ.h>
#include <NTL/mat_ZZ.h>

using namespace std;

BarvinokParameters parameters;
ReadPolyhedronData read_polyhedron_data;
string output_filename;

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

		num = simplicialCone->coefficient * num;
		add(num, denom, tempNum, tempDenom);
	}//for every simple cone.
	for (int i = 2; i <= poly->numOfVars; i++) {
		denom *= i;
	}
	tempNum = GCD(num, denom);
	num /= tempNum;
	denom /= tempNum;

	cout << "VOLUME: " << endl << -num << endl << endl;
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
	set_program_name(argv[0]);

	int i;
	for (i = 1; i < argc; i++) {
		if (read_polyhedron_data.parse_option(argv[i])) {
		} else if (strncmp(argv[i], "--output-cones=", 15) == 0) {
			output_filename = argv[i] + 15;
		} else {
			cerr << "Unknown argument: " << argv[i] << endl;
			exit(1);
		}
	}
	Polyhedron *poly = read_polyhedron_data.read_polyhedron(&parameters);

	parameters.Number_of_Variables = poly->numOfVars;

	//computeUniVolume(poly);
	//computeTriangVolume(poly);
	printRationalFunction(poly);
	return 0;
}

