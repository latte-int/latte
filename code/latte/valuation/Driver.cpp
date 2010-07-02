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
#include <NTL/RR.h>
#include <NTL/mat_ZZ.h>
#include "testEhrhart/BuildRandomPolytope.h"

/* EHRHART INCLUDES */
#include <string.h>
#include <stdio.h>

#include "config.h"
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
//  #include "jesus.h"
#include "preprocess.h"
#include "print.h"
#include "ramon.h"
#include "rational.h"
#include "timing.h"
#include "flags.h"
//#include "testing.h"
#include "IntegralHull.h"
#include "ReadingFile.h"
#include "binarySearchIP.h"
#include "CheckEmpty.h"
#include "ProjectUp.h"

#include "banner.h"
#include "convert.h"
#include "latte_system.h"

#include "gnulib/progname.h"
/* END EHRHART INCLUDES */
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
ZZ dot(vec_ZZ& a, vec_ZZ& b)
{
	ZZ ans = ZZ();
	assert(a.length() == b.length());
	for (int i = 0; i < a.length(); i++)
	{
		ans += a[i] * b[i];
	}
	return ans;
}

//overwrites the num and denom of the first number (the first 2 variables passed in)
//wih the value of adding it to the second number (the next 2 slots)
void add(ZZ& num1, ZZ& denom1, ZZ& num2, ZZ& denom2)
{
	assert(denom1 != 0 && denom2 != 0);
	num2 *= denom1;
	num1 *= denom2;
	denom1 *= denom2;
	num1 += num2;
	ZZ tempNum = GCD(num1, denom1);
	num1 /= tempNum;
	denom1 /= tempNum;
	if (denom1 < 0)
	{
		denom1 *= -1;
		num1 *= -1;
	}//if negative.
}

void computeUniVolume(listCone * triangulatedCones, int numOfVars)
{
	srand((unsigned) time(0));
	vec_ZZ c = vec_ZZ();
	//listCone * triangulatedCones;
	c.SetLength(numOfVars);
	for (int i = 0; i < numOfVars; i++)
		c[i] = rand() % 10000;
	ZZ valutation = ZZ();
	ZZ num = ZZ();
	ZZ denom = ZZ();
	denom = 1;
	ZZ tempNum = ZZ();
	ZZ tempDenom = ZZ();
	vec_ZZ vert = vec_ZZ();
	//	triangulatedCones = decomposeCones(poly->cones, false, parameters);
	//	cout << "c: ";
	//	for (int i = 0; i < poly->numOfVars; i++) {
	//		out << c[i] << " ";
	//	}
	//	cout << endl;
	for (listCone * simplicialCone = triangulatedCones; simplicialCone; simplicialCone
			= simplicialCone->rest)
	{
		tempDenom = 1;
		vert = scaleRationalVectorToInteger(simplicialCone->vertex->vertex,
							numOfVars, tempDenom);
		//		for (int i = 0; i < poly->numOfVars; i++) {
		//			cout << vert[i] << " ";
		//		}
		//		cout << endl;

		//		cout << "eta: " << tempNum << endl;
		tempNum = dot(vert, c);
		tempNum = power(tempNum, numOfVars);
		tempDenom = power(tempDenom, numOfVars);

		for (listVector * currentRay = simplicialCone->rays; currentRay; currentRay
				= currentRay->rest)
		{
			//			for (int i = 0; i < poly->numOfVars; i++) {
			//				cout << currentRay->first[i] << " ";
			//			}
			//			cout << endl;
			tempDenom *= -1 * dot(currentRay->first, c);
			//			cout << "xi: " << tempDenom << endl;
		}//for every ray

		tempNum = simplicialCone->coefficient * tempNum;
		add(num, denom, tempNum, tempDenom);
	}//for every simple cone.
	for (int i = 2; i <= numOfVars; i++)
	{
		denom *= i;
	}
	tempNum = GCD(num, denom);
	num /= tempNum;
	denom /= tempNum;

	cout << "computeUniVolume(): VOLUME: " << endl << num << endl << endl;
	if (denom != 1)
		cout << "/" << endl << endl << denom;
	cout << endl;
}

void computeTriangVolume(listCone * inputCones, int numOfVars)
{
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
	mat.SetDims(numOfVars, numOfVars);
//	listCone *triangulatedCones;

	c.SetLength(numOfVars);
	for (int i = 0; i < numOfVars; i++)
		c[i] = rand() % 10000;

//	for (listCone *currentCone = inputCones; currentCone; currentCone
//			= currentCone->rest)
//	{
//		triangulatedCones
//				= triangulateCone(currentCone, numOfVars, &parameters);

		for (listCone * simplicialCone = inputCones /*triangulatedCones*/; simplicialCone; simplicialCone
				= simplicialCone->rest)
		{

			//find vertex
			vert = scaleRationalVectorToInteger(simplicialCone->vertex->vertex,
					numOfVars, tempDenom);

			//raise f(vertex) to the power of the dimension
			tempNum = dot(vert, c);
			tempNum = power(tempNum, numOfVars);
			tempDenom = power(tempDenom, numOfVars);

			int col = 0;

			for (listVector * currentRay = simplicialCone->rays; currentRay; currentRay
					= currentRay->rest, col++) {
				//divide by the dot product of c and the ray
				tempDenom *= -1 * dot(c, currentRay->first);

				//generate matrix
				for (int row = 0; row < numOfVars; row++) {
					mat[row][col] = currentRay->first[row];
				}//for every component of the ray

			}//for every ray in the simple cone

			//get the determinant
			determinant(det, mat);

			//multiply by the absolute value of the determinant
			tempNum *= abs(det);

			//add current term to the running total
			//cout << "adding " << tempNum << " / " << tempDenom << " to " << num
			//		<< " / " << denom << endl;
			add(num, denom, tempNum, tempDenom);
		}//for every simple cone in the cone

//	}//for every cone
	for (int i = 2; i <= numOfVars; i++)
	{
		denom *= i;
	}
	tempNum = GCD(num, denom);
	num /= tempNum;
	denom /= tempNum;

	cout << "computeTriangVolume(): VOLUME: " << endl << num << endl;
	if (denom != 1)
		cout << "/" << endl << denom;
	cout << endl << to_RR(num) / to_RR(denom) << endl;
}

void printRationalFunction(Polyhedron * poly)
{
	listCone * triangulatedCones;
	vec_ZZ vert = vec_ZZ();
	ZZ temp = ZZ();
	triangulatedCones = decomposeCones(poly->cones, false, parameters);
	cout << "( ";
	for (listCone * simplicialCone = triangulatedCones; simplicialCone; simplicialCone
			= simplicialCone->rest)
	{
		printConeToFile(cout, simplicialCone, poly->numOfVars);
		vert = scaleRationalVectorToInteger(simplicialCone->vertex->vertex,
				poly->numOfVars, temp);
		cout << "( ";
		for (int i = 0; i < poly->numOfVars - 1; i++)
		{
			cout << vert[i];
			if (temp != 1)
				cout << " / " << temp;
			cout << " * c" << i << " + ";
		}
		cout << vert[poly->numOfVars - 1];
		if (temp != 1)
			cout << " / " << temp;
		cout << " * c" << poly->numOfVars - 1 << " ) ^ " << poly->numOfVars
				<< " / ( ";

		if (poly->numOfVars % 2 == 1)
			cout << "-";
		for (listVector * currentRay = simplicialCone->rays; currentRay; currentRay
				= currentRay->rest)
		{
			cout << "( ";
			for (int i = 0; i < poly->numOfVars; i++)
			{
				cout << currentRay->first[i] << " * c" << i;
				if (i != poly->numOfVars - 1)
				{
					cout << " + ";
				}
			}
			cout << " )";
			if (currentRay->rest != NULL)
				cout << " * ";
		}//for every ray
		cout << " ) * ";
		cout << simplicialCone->coefficient;
		if (simplicialCone->rest != NULL)
			cout << " + ";
	}//for every simple cone.
	cout << ") / ( " << poly->numOfVars << "!";
	cout << " )" << endl;
}
/* ----------------------------------------------------------------- */

int main(int argc, char *argv[])
{
#if 1
	set_program_name(argv[0]);

	BarvinokParameters params;
#ifdef SUN
	struct tms tms_buf;
#endif
	float z;
	int i, numOfVars, numOfAllVars, degree = 1;
	unsigned int flags = 0, print_flag = 0, output_cone = 0;
	vec_ZZ dim, v, w;
	int oldnumofvars;
	vec_ZZ *generators;
	char fileName[127], invocation[127], decompose[10], equationsPresent[10],
			assumeUnimodularCones[127], dualApproach[127], taylor[127],
			printfile[127], rationalCone[127], nonneg[127], Memory_Save[127],
			Save_Tri[127], Load_Tri[127], Print[127], inthull[127],
			cddstyle[127], grobner[127], removeFiles[127], command[127],
			maximum[127], Singlecone[127], LRS[127], Vrepresentation[127],
			dilation[127], minimize[127], binary[127], interior[127];
	listVector *matrix, *equations, *inequalities, *rays, *endRays, *tmpRays,
			*matrixTmp;
	vec_ZZ cost;
	listVector *templistVec;
	listCone *cones, *tmp, *tmpcones;

	latte_banner(cerr);

	z = 0;
	//setbuf(stdout,0);

	strcpy(invocation, "Invocation: ");
	strcat(invocation, argv[0]);
	strcat(invocation, " ");
	/*    strcat(invocation,argv[argc-1]); */
	/*    strcat(invocation,"\n\n"); */
	/*    printf(invocation); */

	strcpy(Vrepresentation, "no");
	strcpy(interior, "no");
	strcpy(dilation, "no");
	strcpy(binary, "no");
	strcpy(Singlecone, "no");
	strcpy(removeFiles, "no");
	strcpy(grobner, "no");
	strcpy(maximum, "no");
	strcpy(minimize, "no");
	strcpy(decompose, "yes");
//	strcpy(dualApproach, "yes");
	strcpy(dualApproach, "no");
	strcpy(equationsPresent, "no");
	strcpy(assumeUnimodularCones, "no");
	strcpy(printfile, "no");
	strcpy(taylor, "no");
	strcpy(rationalCone, "no");
	strcpy(nonneg, "no");
	strcpy(Memory_Save, "yes");
	strcpy(Save_Tri, "no");
	strcpy(Load_Tri, "no");
	strcpy(Print, "no");
	strcpy(inthull, "no");
	strcpy(cddstyle, "no");
	strcpy(LRS, "no");

	flags |= DUAL_APPROACH;

	for (i = 1; i < argc - 1; i++)
	{
		strcat(invocation, argv[i]);
		strcat(invocation, " ");

		if (strncmp(argv[i], "simp", 3) == 0)
			strcpy(rationalCone, "yes");

		if (strncmp(argv[i], "cdd", 3) == 0)
			strcpy(cddstyle, "yes");

		if (strncmp(argv[i], "lrs", 3) == 0)
			strcpy(LRS, "yes");

		if (strncmp(argv[i], "trisave", 7) == 0)
		{
			strcpy(Save_Tri, "yes");
			flags |= SAVE;
		}
		if (strncmp(argv[i], "triload", 7) == 0)
		{
			strcpy(Load_Tri, "yes");
			flags |= LOAD;
		}
	}

	/* if cdd and lrs are NOT called. */
	if ((cddstyle[0] == 'n') && (LRS[0] == 'n'))
	{
		if (argc == 2)
		{
			strcpy(printfile, "yes");
			flags |= PRINT;
		} else if (argc == 3)
		{
			if (strncmp(argv[1], "sim", 3) == 0)
				strcpy(rationalCone, "yes");
			else
				strcpy(taylor, "yes");
		} else
		{
			cerr << "Too many arguments.  Check the manual for command line."
					<< endl;
			exit(1);
		}
	}

	/* if cdd is called but and lrs is NOT called. */
	else if ((cddstyle[0] == 'y') && (LRS[0] == 'n'))
	{
		if (argc == 3)
		{
			strcpy(printfile, "yes");//hit
			flags |= PRINT; //hit
		} else if (argc == 4)
		{
			for (i = 1; i < 3; i++)
				if (strncmp(argv[i], "sim", 3) == 0)
					strcpy(rationalCone, "yes");
			if (rationalCone[0] == 'n')
				strcpy(taylor, "yes");
		} else
		{
			cerr << "Too many arguments.  Check the manual for command line."
					<< endl;
			exit(1);
		}
	}

	/* if cdd is NOT called and lrs is called. */
	else if ((cddstyle[0] == 'n') && (LRS[0] == 'y'))
	{
		if (argc == 3)
		{
			strcpy(printfile, "yes");
			flags |= PRINT;
		} else if (argc == 4)
		{
			for (i = 1; i < 3; i++)
				if (strncmp(argv[i], "sim", 3) == 0)
					strcpy(rationalCone, "yes");
			if (rationalCone[0] == 'n')
				strcpy(taylor, "yes");
		} else
		{
			cerr << "Too many arguments.  Check the manual for command line."
					<< endl;
			exit(1);
		}
	}

	/* if cdd and lrs are called. */
	else if ((cddstyle[0] == 'y') && (LRS[0] == 'y'))
	{
		if (argc == 4)
		{
			strcpy(printfile, "yes");
			flags |= PRINT;
		} else if (argc == 5)
		{
			for (i = 1; i < 4; i++)
				if (strncmp(argv[i], "sim", 3) == 0)
					strcpy(rationalCone, "yes");
			// else if(strncmp(argv[3],"sim",3)==0) strcpy(rationalCone,"yes");
			if (rationalCone[0] == 'n')
				strcpy(taylor, "yes");
		} else
		{
			cerr << "Too many arguments.  Check the manual for command line."
					<< endl;
			exit(1);
		}
	}

	flags |= DUAL_APPROACH;//hit
	if (minimize[0] == 'y')
		strcpy(maximum, "yes");
	if (grobner[0] == 'y')
		strcpy(equationsPresent, "yes");
	if (binary[0] == 'y')
	{
		strcpy(maximum, "yes");
		strcpy(Memory_Save, "no");
	}
	if (maximum[0] == 'y')
		strcpy(Memory_Save, "no");
	if (printfile[0] == 'y')
		strcpy(Memory_Save, "no");
	if (rationalCone[0] == 'y') //hit
		strcpy(Memory_Save, "no");
	if (printfile[0] == 'y')
		print_flag = 1;
	if (taylor[0] == 'y')
	{
		degree = atoi(argv[argc - 2]);
	}
	if (rationalCone[0] == 'y')
	{

		//HugInt digit(argv[1]);
		//conv(output_cone, digit.BigInt);
		// User can use only Mode one
		output_cone = 3;
	}
	int dilation_const = 1;

	if (dilation[0] == 'y')
		dilation_const = atoi(argv[argc - 2]);

	if (output_cone > 3)
		output_cone = 0;
	flags |= (output_cone << 1);//hit

	if ((dualApproach[0] == 'y') && (nonneg[0] == 'y') && (equationsPresent[0]
			== 'n'))
	{
		cerr << "You cannot use + and dua at the same time." << endl;
		exit(2);
	}

	if ((cddstyle[0] == 'y') && (Vrepresentation[0] == 'y'))
	{
		cerr << "Use not cdd style and v-representation." << endl;
		exit(2);
	}

	if ((Memory_Save[0] == 'y') && (inthull[0] == 'y'))
	{
		cerr << "You cannot use int and memsave at the same time." << endl;
		exit(3);
	}

	strcat(invocation, argv[argc - 1]); //hit
	strcat(invocation, "\n\n");
	cerr << invocation;
	char costFile[127];
	if (maximum[0] == 'y')
	{
		strcpy(fileName, argv[argc - 1]);
		strcpy(costFile, argv[argc - 1]);
		strcat(costFile, ".cost");
	} else
		strcpy(fileName, argv[argc - 1]); //hit
	//  cerr << fileName << " " << costFile << endl;
	if (maximum[0] == 'y')
	{
		ifstream ReadTest(fileName);
		if (!ReadTest)
		{
			cerr << "Need a polytope input file." << endl;
			exit(2);
		}
		//    cerr << fileName << " " << costFile << endl;
		ifstream INCost(costFile);
		if (!INCost)
		{
			cerr << "Need a cost input file." << endl;
			exit(3);
		}
		int costDim, dummy;
		INCost >> dummy >> costDim;
		cost.SetLength(costDim);
		for (i = 0; i < costDim; i++)
			INCost >> cost[i];
	}
	//strcpy (fileName,"stdin");

	/* Check input file. */
	if (Vrepresentation[0] == 'n') //hit
	{
		if ((cddstyle[0] == 'n') && (grobner[0] == 'n') && (maximum[0] == 'n')
				&& (minimize[0] == 'n'))
		{ //hit
			CheckInputFile(fileName);
			CheckLength(fileName, equationsPresent);
		}
		if (minimize[0] == 'y')
			strcpy(maximum, "yes");

		if ((cddstyle[0] == 'n') && (grobner[0] == 'n') && (maximum[0] == 'y'))
		{
			CheckInputFile(fileName);
			CheckLength(fileName, equationsPresent);
		}

		if (cddstyle[0] == 'y')
		{
			CheckInputFileCDDRep(argv[argc - 1]);
			CheckInputFileCDDRep1(argv[argc - 1]);
			CheckInputFileCDDRep3(argv[argc - 1]);
			CheckInputFileCDDRep4(argv[argc - 1]);
		}
	} else
		CheckInputFileVrep(fileName);
	CheckEmpty(fileName);//hit
	//vector cost;
	/* Read problem data. */
	if ((cddstyle[0] == 'n') && (Vrepresentation[0] == 'n'))
		CheckRed(fileName, equationsPresent, maximum, nonneg, interior,
				dilation, dilation_const);//hit
	//file, yes, no, no, no, no, 1

	dilation_const = 1;
	if ((cddstyle[0] == 'n') && (grobner[0] == 'n'))//hit ?
		readLatteProblem(fileName, &equations, &inequalities, equationsPresent,
				&numOfVars, nonneg, dualApproach, grobner, Vrepresentation);
	// if(minimize[0] == 'y') cost = -cost;
	if (cddstyle[0] == 'y')
	{
		int tmpoutput;
		CDDstylereadLatteProblem(fileName, &equations, &inequalities,
				equationsPresent, &numOfVars, nonneg, dualApproach, taylor,
				degree, rationalCone, tmpoutput, Memory_Save,
				assumeUnimodularCones, inthull, grobner);
		output_cone = 3;
		if (dualApproach[0] == 'y')
		{
			flags |= DUAL_APPROACH;
		}
	}

	// cerr << grobner << endl;
	vec_ZZ holdCost;
	if (minimize[0] == 'y')
		cost = -cost;
	holdCost = cost;//hit
	//cerr <<"Cost is: " << cost << endl;
	vec_RR holdcost_RR;
	holdcost_RR.SetLength(holdCost.length());
	for (i = 0; i < holdCost.length(); i++)
		conv(holdcost_RR[i], holdCost[i]);

	if (minimize[0] == 'y')
		holdcost_RR = -holdcost_RR;
	if (grobner[0] == 'y')
	{

		CheckGrobner(fileName, cddstyle);
		SolveGrobner(fileName, nonneg, dualApproach, grobner, equationsPresent,
				cddstyle);
	} else
	{//hit
		if ((dualApproach[0] == 'y') && (nonneg[0] == 'y')
				&& (equationsPresent[0] == 'n'))
		{
			cerr << "You cannot use + and dua at the same time." << endl;
			exit(2);
		}

		if ((Memory_Save[0] == 'y') && (inthull[0] == 'y'))
		{
			cerr << "You cannot use int and memsave at the same time." << endl;
			exit(3);
		}

		if ((Vrepresentation[0] == 'y') && (equationsPresent[0] == 'y'))
		{
			cerr << "You cannot use vrep and equ at the same time." << endl;
			exit(4);
		}

		numOfVars--;//hit twice?
		/* Binary seach IP*/

		if (binary[0] == 'y')
		{
			cerr << "The number of optimal solutions: " << binarySearch(
					equations, inequalities, cost, numOfVars, minimize) << endl;
			cerr << "Time: " << GetTime() << endl;
			exit(0);
		}

		numOfAllVars = numOfVars;
		mat_ZZ ProjU;
		ProjU.SetDims(numOfVars, numOfVars);
		oldnumofvars = numOfVars;
		generators = createArrayVector(numOfVars);
		if (equationsPresent[0] == 'y')//hit
		{
			/*    if(grobner[0] == 'y')
			 {
			 matrixTmp=Grobner(equations,inequalities,&generators,&numOfVars, &templistVec, oldnumofvars);

			 }*/
			matrixTmp = preprocessProblem(equations, inequalities, &generators,
					&numOfVars, cost, ProjU, interior, dilation_const);
			templistVec = transformArrayBigVectorToListVector(ProjU,
					ProjU.NumCols(), ProjU.NumRows());
		} else
		{
			dilateListVector(inequalities, numOfVars, dilation_const);
			matrixTmp = inequalities;
		}
		if ((dualApproach[0] == 'y') && (equationsPresent[0] == 'y'))
		{
			matrix = TransformToDualCone(matrixTmp, numOfVars);//hit
		} else
		{
			matrix = matrixTmp;
		}
		/* Now matrix contains the new inequalities. */
		RR LP_OPT;
		cerr << "\nTime: " << GetTime() << " sec\n\n";
		//   cerr << "Project down cost function: " << cost << endl;
		vec_RR Rat_solution, tmp_den, tmp_num;
		mat_RR ProjU_RR;
		ProjU_RR.SetDims(ProjU.NumRows(), ProjU.NumCols());
		for (i = 0; i < ProjU.NumRows(); i++) //hit
			for (int j = 0; j < ProjU.NumCols(); j++)
				conv(ProjU_RR[i][j], ProjU[i][j]);
		//cerr << ProjU << ProjU_RR << endl;
		Rat_solution.SetLength(numOfVars);
		tmp_den.SetLength(numOfVars);
		tmp_num.SetLength(numOfVars);
		/* Compute vertices and edges. */
		rationalVector* LP_vertex;
		if ((dualApproach[0] == 'n') && (Vrepresentation[0] == 'n'))
		{
			if (LRS[0] == 'n')
				tmpcones = computeVertexCones(fileName, matrix, numOfVars);
			else
				tmpcones
						= computeVertexConesViaLrs(fileName, matrix, numOfVars);
			if (maximum[0] == 'y')
			{
				LP_vertex = LP(matrix, cost, numOfVars);
				vec_RR Rat_cost;
				Rat_cost.SetLength(numOfVars);
				for (i = 0; i < numOfVars; i++)
				{
					conv(tmp_num[i], LP_vertex->numerators()[i]);
					conv(tmp_den[i], LP_vertex->denominators()[i]);
					Rat_solution[i] = tmp_num[i] / tmp_den[i];
					conv(Rat_cost[i], cost[i]);
				}
				if (Singlecone[0] == 'y')
					cones = CopyListCones(tmpcones, numOfVars, LP_vertex);
				else
					cones = tmpcones;
				if (lengthListCone(cones) == 1)
					cerr << "\nWe found a single vertex cone for IP.\n" << endl;
				cerr << "A vertex which we found via LP is: "
						<< ProjectingUpRR(ProjU_RR, Rat_solution, numOfVars)
						<< endl;
				//printRationalVector(LP_vertex, numOfVars);
				LP_OPT = Rat_cost * Rat_solution; //cerr << cost << endl;
				cerr << "The LP optimal value is: " << holdcost_RR
						* ProjectingUpRR(ProjU_RR, Rat_solution, numOfVars)
						<< endl;
			} else
			{
				cones = tmpcones;
				cerr << "\nThe polytope has " << lengthListCone(cones)
						<< " vertices.\n";
				//system_with_error_check("rm -f numOfLatticePoints");
				cerr << endl;
			}
		}

		/* Reading from the vertex representation. */

		if (Vrepresentation[0] == 'y')
			cones = computeVertexConesFromVrep(fileName, numOfVars);

		/* Compute triangulation or decomposition of each vertex cone. */

		//CONES SHOULD BE FULLY DIM. --not true?


		if (dualApproach[0] == 'y')
		{ //hit
			cones = createListCone();
			cones->vertex = new Vertex(createRationalVector(numOfVars));
			rays = createListVector(createVector(numOfVars));
			endRays = rays;
			tmpRays = matrix;
			while (tmpRays)
			{
				v = createVector(numOfVars);
				for (i = 0; i < numOfVars; i++)
					v[i] = -(tmpRays->first)[i + 1];
				endRays->rest = createListVector(v);
				endRays = endRays->rest;
				tmpRays = tmpRays->rest;
			}
			cones->rays = rays->rest;

			if (Memory_Save[0] == 'n')
			{//hit prints Decomposing all cones. Triangulating cone... done. All cones have been decomposed.80 cones in total.600


				BarvinokParameters brandonParameters;

				brandonParameters.Flags = flags;
				brandonParameters.Number_of_Variables = numOfVars;
				brandonParameters.max_determinant = 1;
				brandonParameters.File_Name = fileName;
				brandonParameters.decomposition
						= BarvinokParameters::DualDecomposition;

				ofstream file;
				file.open("printAllTheConesDual.txt");
				file << "cones that are passed into decomposedCones()" << endl;
				for (listCone * c = cones; c; c = c->rest)
					printConeToFile(file, c, numOfVars);
				//PolytopeValuation polytopeValuation(cones, numOfVars, &brandonParameters);
				//RationalNTL rVolume =  polytopeValuation.findVolume();
				//cout << "polytopeValuation.findVolume(): VOLUME" << rVolume << " = " << rVolume.to_RR() << endl;

				cout << "CALLING DECOMPOSED CONES\n" << endl;
				cones = decomposeCones(cones, numOfVars, flags, fileName, 1,
						false, BarvinokParameters::DualDecomposition);

				file << "after call to decompsost cones " << endl;
				for (listCone * c = cones; c; c = c->rest)
					printConeToFile(file, c, numOfVars);
				computeUniVolume(cones, numOfVars);
				computeTriangVolume(cones, numOfVars);

				//cout << "PRINTING DECOMPOSED CONES THAT SHOULD BE UNIMODULAR" << endl;


			}

			else
			{
				decomposeCones_Single(cones, numOfVars, degree, flags,
						fileName, 1, false,
						BarvinokParameters::DualDecomposition);
			}

		}

		else
		{

			if (assumeUnimodularCones[0] == 'n')
			{
				if (decompose[0] == 'y')
				{
					if (Memory_Save[0] == 'n')
					{
						ofstream file;
						file.open("printAllTheConesNotDual.txt");
						file << "cones that are passed into decomposedCones()"
								<< endl;
						for (listCone * c = cones; c; c = c->rest)
							printConeToFile(file, c, numOfVars);

						BarvinokParameters brandonParameters;
						brandonParameters.Flags = flags;
						brandonParameters.Number_of_Variables = numOfVars;
						brandonParameters.max_determinant = 1;
						brandonParameters.File_Name = fileName;
						brandonParameters.decomposition
								= BarvinokParameters::DualDecomposition;

						PolytopeValuation polytopeValuation(cones, numOfVars,
								&brandonParameters);
						RationalNTL rVolume = polytopeValuation.findVolume();
						cout << "VOLUME (ehrhart)" << rVolume << " = "
								<< rVolume.to_RR() << endl;

						cones = decomposeCones(cones, numOfVars, flags,
								fileName, 1, true,
								BarvinokParameters::DualDecomposition);
						computeUniVolume(cones, numOfVars);
						computeTriangVolume(cones, numOfVars);

						//int k = 0;
						file << "after call to decompsost cones " << endl;
						for (listCone * c = cones; c; c = c->rest)
						{
							//cout << ++k << endl;
							printConeToFile(file, c, numOfVars);
						}
						cout << "printed ok" << endl;

						// Iterator through simplicial cones, DFS
					} else
					{
						decomposeCones_Single(cones, numOfVars, degree, flags,
								fileName, 1, true,
								BarvinokParameters::DualDecomposition);
					}//else
				}
			}
		}

		/* Compute points in parallelepipeds, unless we already did using memsave version!  */

		if (Memory_Save[0] == 'n')
		{ //hit
			cerr
					<< "Computing the points in the Parallelepiped of the unimodular Cones."
					<< endl;
			tmp = cones;
			int Cones_Processed_Count = 0;
			while (tmp)
			{
				tmp->latticePoints = pointsInParallelepiped(tmp, numOfVars,
						&params);
				tmp = tmp->rest;

				Cones_Processed_Count++;

				if ((Cones_Processed_Count % 1000) == 0)
					cerr << Cones_Processed_Count << " cones processed."
							<< endl;
			}
		}

		if (grobner[0] == 'y')
		{

			cones = ProjectUp(cones, oldnumofvars, numOfVars, templistVec);
			numOfVars = oldnumofvars;

		}
		if (Print[0] == 'y')
			printListCone(cones, numOfVars);

		if (inthull[0] == 'y')
			;
		//    printListVector(IntegralHull(cones, inequalities, equations,  numOfVars),numOfVars);

		if (maximum[0] == 'y')
		{
			listCone * Opt_cones;
			if (Singlecone[0] == 'n')
			{
				Opt_cones = CopyListCones(cones, numOfVars);
				ZZ NumOfLatticePoints; //printListCone(Opt_cones, numOfVars);
				NumOfLatticePoints = Residue(Opt_cones, numOfVars);
				cerr << "Finished computing a rational function. " << endl;
				cerr << "Time: " << GetTime() << " sec." << endl;
				if (IsZero(NumOfLatticePoints) == 1)
				{
					cerr
							<< "Integrally empty polytope.  Check the right hand side."
							<< endl;
					exit(0);
				} else
				{
					int singleCone = 0;
					if (Singlecone[0] == 'y')
						singleCone = 1;
					vec_ZZ Opt_solution;
					if (minimize[0] == 'y')
						holdCost = -holdCost;
					Opt_solution = SolveIP(cones, matrix, cost, numOfVars,
							singleCone);
					if (minimize[0] == 'y')
						cost = -cost;
					cerr << "An optimal solution for " << holdCost << " is: "
							<< ProjectingUp(ProjU, Opt_solution, numOfVars)
							<< "." << endl;
					cerr << "The projected down opt value is: " << cost
							* Opt_solution << endl;
					cerr << "The optimal value is: " << holdCost
							* ProjectingUp(ProjU, Opt_solution, numOfVars)
							<< "." << endl;
					ZZ IP_OPT;
					IP_OPT = cost * Opt_solution;
					RR tmp_RR;

					conv(tmp_RR, cost * Opt_solution);
					// cerr << tmp_RR << " " << LP_OPT << endl;
					if (minimize[0] == 'y')
						LP_OPT = -LP_OPT;
					cerr << "The gap is: " << abs(tmp_RR - LP_OPT) << endl;
					cerr << "Computation done." << endl;
					cerr << "Time: " << GetTime() << " sec." << endl;
					strcpy(command, "rm -f ");
					strcat(command, fileName);
					strcat(command, ".ext");
					system_with_error_check(command);

					strcpy(command, "rm -f ");
					strcat(command, fileName);
					strcat(command, ".cdd");
					system_with_error_check(command);

					strcpy(command, "rm -f ");
					strcat(command, fileName);
					strcat(command, ".ead");
					system_with_error_check(command);

					if (cddstyle[0] == 'n')
					{
						strcpy(command, "rm -f ");
						strcat(command, fileName);
						system_with_error_check(command);
					}

					exit(0);
				}
			} else
			{
				int singleCone = 0;
				if (Singlecone[0] == 'y')
					singleCone = 1;
				vec_ZZ Opt_solution;
				if (minimize[0] == 'y')
					holdCost = -holdCost;
				Opt_solution = SolveIP(cones, matrix, cost, numOfVars,
						singleCone);
				cerr << "An optimal solution for " << holdCost << " is: "
						<< ProjectingUp(ProjU, Opt_solution, numOfVars) << "."
						<< endl;
				if (minimize[0] == 'y')
					cost = -cost;
				cerr << "The projected down opt value is: " << cost
						* Opt_solution << endl;
				cerr << "The optimal value is: " << holdCost * ProjectingUp(
						ProjU, Opt_solution, numOfVars) << "." << endl;
				ZZ IP_OPT;
				IP_OPT = cost * Opt_solution;
				RR tmp_RR;
				conv(tmp_RR, IP_OPT);
				// cerr << cost * Opt_solution << endl;
				if (minimize[0] == 'y')
					LP_OPT = -LP_OPT;
				cerr << "The gap is: " << abs(tmp_RR - LP_OPT) << endl;
				cerr << "Computation done." << endl;
				cerr << "Time: " << GetTime() << " sec." << endl;
				strcpy(command, "rm -f ");
				strcat(command, fileName);
				strcat(command, ".ext");
				system_with_error_check(command);

				strcpy(command, "rm -f ");
				strcat(command, fileName);
				strcat(command, ".cdd");
				system_with_error_check(command);

				strcpy(command, "rm -f ");
				strcat(command, fileName);
				strcat(command, ".ead");
				system_with_error_check(command);

				if (cddstyle[0] == 'n')
				{
					strcpy(command, "rm -f ");
					strcat(command, fileName);
					system_with_error_check(command);
				}

				exit(0);
			}
		} else
		{
			/* if(maximum[0] == 'y') {
			 listCone * Opt_cones;
			 if(Singlecone[0] == 'n'){
			 Opt_cones = CopyListCones(cones, numOfVars);
			 ZZ NumOfLatticePoints; //printListCone(Opt_cones, numOfVars);
			 NumOfLatticePoints = Residue(Opt_cones, numOfVars);
			 cerr <<"Finished computing a rational function. " << endl;
			 cerr <<"Time: " << GetTime() << " sec." << endl;
			 if(IsZero(NumOfLatticePoints) == 1){
			 cerr<<"Integrally empty polytope.  Check the right hand side."<< endl;
			 exit(0);}
			 else{
			 int singleCone = 0;
			 if(Singlecone[0] == 'y') singleCone = 1;
			 vec_ZZ Opt_solution;
			 if(minimize[0] == 'y') holdCost = -holdCost;
			 Opt_solution = SolveIP(cones, inequalities, equations,  cost, numOfVars, singleCone);
			 cerr << "An optimal solution for " <<  holdCost << " is: " << ProjectingUp(ProjU, Opt_solution, numOfVars) << "." << endl;
			 cerr << "The projected down opt value is: " << cost * Opt_solution << endl;
			 cerr <<"The optimal value is: " << holdCost * ProjectingUp(ProjU, Opt_solution, numOfVars) << "." << endl;
			 ZZ IP_OPT; IP_OPT = cost*Opt_solution;
			 RR tmp_RR;
			 conv(tmp_RR, IP_OPT);
			 // cerr << cost * Opt_solution << endl;
			 cerr <<"The gap is: "<< abs(tmp_RR - LP_OPT) << endl;
			 cerr << "Computation done." << endl;
			 cerr <<"Time: " << GetTime() << " sec." << endl;
			 exit(0);
			 }
			 }
			 else{
			 int singleCone = 0;
			 if(Singlecone[0] == 'y') singleCone = 1;
			 vec_ZZ Opt_solution;
			 if(minimize[0] == 'y') holdCost = -holdCost;
			 Opt_solution = SolveIP(cones, inequalities, equations,  cost, numOfVars, singleCone);
			 cerr << "An optimal solution for " <<  holdCost << " is: " << ProjectingUp(ProjU, Opt_solution, numOfVars) << "." << endl;
			 cerr << "The projected down opt value is: " << cost * Opt_solution << endl;
			 cerr <<"The optimal value is: " << holdCost * ProjectingUp(ProjU, Opt_solution, numOfVars) << "." << endl;
			 ZZ IP_OPT; IP_OPT = cost*Opt_solution;
			 RR tmp_RR;
			 conv(tmp_RR, IP_OPT);
			 // cerr << cost * Opt_solution << endl;
			 cerr <<"The gap is: "<< abs(tmp_RR - LP_OPT) << endl;
			 cerr << "Computation done." << endl;
			 cerr <<"Time: " << GetTime() << " sec." << endl;
			 exit(0);
			 }
			 }else{*/
			if (Memory_Save[0] == 'n')
			{//hit

				if (dualApproach[0] == 'n')
				{ //not hit ?
					cerr << "Creating generating function.\n";
					//printListVector(templistVec, oldnumofvars); cerr << ProjU << endl;
					if (equationsPresent[0] == 'y')
					{
						cones = ProjectUp(cones, oldnumofvars, numOfVars,
								templistVec);
						numOfVars = oldnumofvars;
					}

					createGeneratingFunctionAsMapleInput(fileName, cones,
							numOfVars);
				}
				//printListCone(cones, numOfVars);

				cerr << "Printing decomposed cones to decomposed_cones."
						<< endl;
				printListConeToFile("decomposed_cones", cones, numOfVars);

				if (dualApproach[0] == 'n')
				{ //hit
					cerr << "Starting final computation.\n";
					cerr << endl << "****  The number of lattice points is: "
							<< Residue(cones, numOfVars) << "  ****" << endl
							<< endl;
				}

				if (dualApproach[0] == 'y')
				{
					cerr << "Starting final computation.\n";
					//cerr << "output_cone: " << output_cone;
					ResidueFunction(cones, numOfVars, print_flag, degree,
							output_cone, &params);
					//  Else we have already computed the residue.
				}

			}
		}

		if (rationalCone[0] == 'y')
		{
			cerr << endl << "Rational function written to " << argv[argc - 1]
					<< ".rat" << endl << endl;
			strcpy(command, "mv ");
			strcat(command, "simplify.sum ");
			strcat(command, argv[argc - 1]);
			strcat(command, ".rat");
			system_with_error_check(command);
		}

		if (printfile[0] == 'y')
		{ //hit
			cerr << endl << "Rational function written to " << argv[argc - 1]
					<< ".rat" << endl << endl;
			strcpy(command, "mv ");
			strcat(command, "func.rat ");
			strcat(command, argv[argc - 1]);
			strcat(command, ".rat");
			system_with_error_check(command);
		}
		if (removeFiles[0] == 'y')
		{

			strcpy(command, "rm -f ");
			strcat(command, fileName);
			strcat(command, ".ext");
			system_with_error_check(command);

			strcpy(command, "rm -f ");
			strcat(command, fileName);
			strcat(command, ".cdd");
			system_with_error_check(command);

			if (Memory_Save[0] == 'n')
			{
				strcpy(command, "rm -f ");
				strcat(command, fileName);
				strcat(command, ".maple");
				system_with_error_check(command);
			}

			if (cddstyle[0] == 'n')
			{
				strcpy(command, "rm -f ");
				strcat(command, fileName);
				strcat(command, ".ead");
				system_with_error_check(command);
			}
		}
	}

	if (cddstyle[0] == 'n')
	{ //hit
		strcpy(command, "rm -f ");
		strcat(command, fileName);
		system_with_error_check(command);
	}

	cerr << "Computation done. " << endl;
	cerr << "Time: " << GetTime() << " sec\n\n";

	return (0);
#endif
	BuildRandomPolytope *buildPolytope = 0;
#if 0

	if ( argc > 2 )
	{
		buildPolytope = new BuildRandomPolytope(atoi(argv[1]));
		buildPolytope->setComments("Making random integer polytope for volume testing");
		buildPolytope->buildPolymakeFile(atoi(argv[2])); //make the file
		buildPolytope->callPolymake(); //run polymake
		buildPolytope->findVolumeWithPolymake(); //run polymake for the volume
		buildPolytope->convertFacetEquations(); //fix facet equations
		buildPolytope->printFacetEquationsForLattE(); //make latte file.

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

#endif
	Polyhedron *poly = read_polyhedron_data.read_polyhedron(&parameters);
	//cout << "Finished reading poly. data." << endl;
	if (buildPolytope)
		delete buildPolytope;

	parameters.Number_of_Variables = poly->numOfVars;

	//cout << "PRINT VERTICES" << endl;
	//for (listCone * c = poly->cones; c; c = c->rest)
	//	printRationalVectorToFile(cout, c->vertex->vertex, poly->numOfVars);
	//cout << "END PRINT VERTICES" << endl;

#if 0
	PolytopeValuation polytopeValuation(poly, &parameters);
	RationalNTL rVolume = polytopeValuation.findVolume();
	cout << "TOTAL VOLUME (from the class)(rational) = " << rVolume << " = " << rVolume.to_RR() << endl;
	cout << "TOTAL VOLUME (from the class)(float) = " << polytopeValuation.findVolume_old() << endl;
	exit(0);
#endif

#if 0
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
	//printRationalFunction(poly);
	return 0;
}

