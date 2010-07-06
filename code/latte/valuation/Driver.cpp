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

void computeVolume(listCone * cones, BarvinokParameters &myParameters,
		const char *valuationType);

BarvinokParameters parameters;
ReadPolyhedronData read_polyhedron_data;
string output_filename;


/* ----------------------------------------------------------------- */

/**
 * Computes the volume of the polytope.
 */
void computeVolume(listCone * cones, BarvinokParameters &myParameters,
		const char *valuationAlg, const char * print)
{
	RationalNTL ans1, ans2;
	if (strcmp(valuationAlg, "triangulate") == 0 || strcmp(valuationAlg, "all")
			== 0)
	{
		PolytopeValuation polytopeValuation(cones,
				PolytopeValuation::VertexRayCones,
				myParameters.Number_of_Variables, myParameters);
		ans1 = polytopeValuation.findVolume(PolytopeValuation::DeterminantVolume);

		cout << "(using triangulation-determinant method) VOLUME " << ans1 << endl;
	}//if triangulate

	if (strncmp(valuationAlg, "lawrence", 8) == 0
			|| strcmp(valuationAlg, "all") == 0)
	{

		PolytopeValuation polytopeValuation(cones,
				PolytopeValuation::VertexRayCones,
				myParameters.Number_of_Variables, myParameters);
		ans2 = polytopeValuation.findVolume(PolytopeValuation::LawrenceVolume);

		cout << "(using Lawrence method) VOLUME " << ans2 << endl;

		if (*print == 'y')
			polytopeValuation.printLawrenceVolumeFunction(); //print the lawrence rational function.

	}//if lawrence

	if ( strcmp(valuationAlg, "all") == 0 && ans1 != ans2)
	{
		cout << "Driver.cpp: the two methods are different." << ans1 << "!=" << ans2 << endl;
		exit(1);
	}//if error.

}//computeVolume

int mainValuationDriver(char *argv[], int argc)
{

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
	vec_ZZ *generators = 0;
	char fileName[127], invocation[127], decompose[10], equationsPresent[10],
			assumeUnimodularCones[127], dualApproach[127], taylor[127],
			printfile[127], rationalCone[127], nonneg[127], Memory_Save[127],
			Save_Tri[127], Load_Tri[127], Print[127], inthull[127],
			cddstyle[127], grobner[127], removeFiles[127], command[127],
			maximum[127], Singlecone[127], LRS[127], Vrepresentation[127],
			dilation[127], minimize[127], binary[127], interior[127];
	char valuationType[127], valuationAlg[127], printLawrenceFunction[127];
	listVector *matrix = 0, *equations = 0, *inequalities = 0, *rays = 0, *endRays = 0, *tmpRays = 0,
			*matrixTmp = 0;
	vec_ZZ cost;
	listVector *templistVec = 0;
	listCone *cones = 0, *tmp = 0, *tmpcones = 0;

	latte_banner(cerr);

	z = 0;

	strcpy(invocation, "Invocation: ");
	strcat(invocation, argv[0]);
	strcat(invocation, " ");

	strcpy(Vrepresentation, "no");
	strcpy(interior, "no");
	strcpy(dilation, "no");
	strcpy(binary, "no");
	strcpy(Singlecone, "no");
	strcpy(removeFiles, "yes");
	strcpy(grobner, "no");
	strcpy(maximum, "no");
	strcpy(minimize, "no");
	strcpy(decompose, "yes");
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
	strcpy(valuationType, "volume");
	strcpy(valuationAlg, "all");
	strcpy(printLawrenceFunction, "no");

	flags |= DUAL_APPROACH;



	if ( argc == 2 && !strcmp(argv[1], "help"))
	{
		cout << "usage: valuation-exe [valuation  type] [valuation algorithm] <latte file>\n";
		cout << "valuation types: volume\n";
		cout << "  volume algorithm: [lawrence  <printLawrenceFunction> | triangulate | all]\n";
		exit(0);
	}//if need help.

	for (i = 1; i < argc - 1; i++) //i = 1...argc-1 because we do not want to process the file name. The file name could be "simp"
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
		if ( strcmp(argv[i], "keepFiles") == 0)
			strcpy(removeFiles, "no");
		if (strncmp(argv[i], "lawrence", 8) == 0)
			strcpy(valuationAlg, "lawrence");
		if (strcmp(argv[i], "triangulate") == 0)
			strcpy(valuationAlg, "triangulate");
		if (strcmp(argv[i], "all") == 0)
			strcpy(valuationAlg, "all");
		if ( strcmp(argv[i], "printLawrenceFunction") == 0)
			strcpy(printLawrenceFunction, "yes");

	}//for i.
#if 0
	/* if cdd and lrs are NOT called. */
	if ((cddstyle[0] == 'n') && (LRS[0] == 'n'))
	{//hit
		if (argc == 2)
		{
			strcpy(printfile, "yes");
			flags |= PRINT;//hit
		} else if (argc == 3)
		{
			if (strncmp(argv[1], "sim", 3) == 0)
				strcpy(rationalCone, "yes");
			else
				strcpy(taylor, "yes");
		}
		//
		//else
		//{
		//	cerr << "Too many arguments.  Check the manual for command line."
		//			<< endl;
		//	exit(1);
		//}
	}


	/* if cdd is called but and lrs is NOT called. */
	else if ((cddstyle[0] == 'y') && (LRS[0] == 'n'))
	{
		if (argc == 3)
		{
			strcpy(printfile, "yes");
			//flags |= PRINT; //hit
			flags &= ~PRINT; //don't print.
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
			flags |= PRINT;//hit
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
#endif

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
	if (rationalCone[0] == 'y')
		strcpy(Memory_Save, "no");
	if (printfile[0] == 'y')
		print_flag = 1;

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
	flags |= (output_cone << 1);

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

	strcat(invocation, argv[argc - 1]);
	strcat(invocation, "\n\n");
	cerr << invocation;
	char costFile[127];

	strcpy(fileName, argv[argc - 1]);

	/* Check input file. */
	if (Vrepresentation[0] == 'n')
	{//hit
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

	/* Read problem data. */
	if ((cddstyle[0] == 'n') && (Vrepresentation[0] == 'n'))//hit
		CheckRed(fileName, equationsPresent, maximum, nonneg, interior,
				dilation, dilation_const);
	//file, yes, no, no, no, no, 1

	dilation_const = 1;
	if ((cddstyle[0] == 'n') && (grobner[0] == 'n'))//hit
		readLatteProblem(fileName, &equations, &inequalities, equationsPresent,
				&numOfVars, nonneg, dualApproach, grobner, Vrepresentation);

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

	if ((dualApproach[0] == 'y') && (nonneg[0] == 'y') && (equationsPresent[0]
			== 'n'))
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

	numOfVars--;
	numOfAllVars = numOfVars;
	mat_ZZ ProjU;
	ProjU.SetDims(numOfVars, numOfVars);
	oldnumofvars = numOfVars;
	generators = createArrayVector(numOfVars);
	if (equationsPresent[0] == 'y') //hit
	{
		matrixTmp = preprocessProblem(equations, inequalities, &generators,
				&numOfVars, cost, ProjU, interior, dilation_const);
		templistVec = transformArrayBigVectorToListVector(ProjU,
				ProjU.NumCols(), ProjU.NumRows());
	} else
	{
		dilateListVector(inequalities, numOfVars, dilation_const);
		matrixTmp = inequalities;
	}

	matrix = matrixTmp;

	/* Now matrix contains the new inequalities. */

	RR LP_OPT;
	cerr << "\nTime: " << GetTime() << " sec\n\n";
	vec_RR Rat_solution, tmp_den, tmp_num;
	mat_RR ProjU_RR;
	ProjU_RR.SetDims(ProjU.NumRows(), ProjU.NumCols());
	for (i = 0; i < ProjU.NumRows(); i++)
		for (int j = 0; j < ProjU.NumCols(); j++)
			conv(ProjU_RR[i][j], ProjU[i][j]);
	//cerr << ProjU << ProjU_RR << endl;
	Rat_solution.SetLength(numOfVars);
	tmp_den.SetLength(numOfVars);
	tmp_num.SetLength(numOfVars);

	/* Compute vertices and edges. */

	rationalVector* LP_vertex;
	if ((dualApproach[0] == 'n') && (Vrepresentation[0] == 'n'))
	{//hit
		if (LRS[0] == 'n')//hit
			tmpcones = computeVertexCones(fileName, matrix, numOfVars);
		else
			tmpcones = computeVertexConesViaLrs(fileName, matrix, numOfVars);

		cones = tmpcones;
		cerr << "\nThe polytope has " << lengthListCone(cones)
				<< " vertices.\n";
		//system_with_error_check("rm -f numOfLatticePoints");
		cerr << endl;
	}

	/* Reading from the vertex representation. */

	if (Vrepresentation[0] == 'y') //not hit
		cones = computeVertexConesFromVrep(fileName, numOfVars);

	/* Compute triangulation or decomposition of each vertex cone. */

	if (dualApproach[0] == 'n') // hit
	{
		if (assumeUnimodularCones[0] == 'n')
		{
			if (decompose[0] == 'y')
			{

				BarvinokParameters myParameters;
				myParameters.Flags = flags;
				myParameters.Number_of_Variables = numOfVars;
				myParameters.max_determinant = 1;
				myParameters.File_Name = fileName;
				myParameters.decomposition
						= BarvinokParameters::DualDecomposition;

				if (strcmp(valuationType, "volume") == 0)
				{
					computeVolume(cones, myParameters, valuationAlg, printLawrenceFunction);
				} else
				{
					cout << "Ops, " << valuationType << " is not supported" << endl;
					exit(1);
				}//if not finding volume.

			}//if decompose
		}//assumeUnimodularCones
	}//dualApproach


//	if (rationalCone[0] == 'y')
//	{
//		cerr << endl << "Rational function written to " << argv[argc - 1]
//				<< ".rat" << endl << endl;
//		strcpy(command, "mv ");
//		strcat(command, "simplify.sum ");
//		strcat(command, argv[argc - 1]);
//		strcat(command, ".rat");
//		system_with_error_check(command);
//	}

//	if (printfile[0] == 'y')
//	{
//		cerr << endl << "Rational function written to " << argv[argc - 1]
//				<< ".rat" << endl << endl;
//		strcpy(command, "mv ");
//		strcat(command, "func.rat ");
//		strcat(command, argv[argc - 1]);
//		strcat(command, ".rat");
//		system_with_error_check(command);
//	}

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

		if (cddstyle[0] == 'n')
		{
			strcpy(command, "rm -f ");
			strcat(command, fileName);
			strcat(command, ".ead");
			system_with_error_check(command);

			strcpy(command, "rm -f ");
			strcat(command, fileName);
			system_with_error_check(command);
		}
	}// remove the latte_nonredundant_input files.


	//free up memory.
	//if ( matrix) 		freeListVector(matrix);
	if ( equations) 	freeListVector(equations);
	if (inequalities) 	freeListVector(inequalities);
	//if (rays) 			freeListVector(rays);
	//if (endRays) 		freeListVector(endRays);
	//if (tmpRays)		freeListVector(tmpRays);
	//if (matrixTmp)		freeListVector(matrixTmp);
	//if (templistVec)	freeListVector(templistVec);

	//if (generators)		delete generators;


	if (cones) 			freeListCone(cones);
	//if (tmp) 			freeListCone(tmp);
	//if (tmpcones) 		freeListCone(tmpcones);

	return (0);
}//mainValuationDrivver()




void runOneTest(int ambientDim, int numPoints)
{
	char * argv[] = {"runTests()", "all", 0};
	stringstream comments;
	comments << "Making random integer polytope with " << numPoints << " points in R^" << ambientDim<< " for volume testing";

	BuildRandomPolytope buildPolytope(ambientDim);
	buildPolytope.setComments(comments.str().c_str());
	buildPolytope.setIntegerPoints(false); //make random rational points.
	buildPolytope.buildPolymakeFile(numPoints); //make the file
	buildPolytope.callPolymake(); //run polymake
	buildPolytope.findVolumeWithPolymake(); //run polymake for the volume
	buildPolytope.convertFacetEquations(); //fix facet equations
	buildPolytope.printFacetEquationsForLattE(); //make latte file.
	cout << comments.str().c_str();
	//buildPolytope.findEhrhardPolynomial();
	//buildPolytope.findVolumeWithPolymake();

	string file = buildPolytope.getLatteFile();


	char * sFile = new char[file.size() + 1];
	strcpy(sFile, file.c_str());
	argv[2] = sFile;
	mainValuationDriver(argv, 3);
	delete [] sFile;
}//RunOneTest

void runTests()
{
	int startAmbientDim = 6, endAmbientDim = 50;
	int pointStepSize = 5;


	for(int ambientDim = startAmbientDim; ambientDim < endAmbientDim; ambientDim = ambientDim + 3)
	{
		for(int numberPoints = startAmbientDim / 2; numberPoints < startAmbientDim/2 + startAmbientDim + 10; numberPoints = numberPoints + pointStepSize)
			runOneTest(ambientDim, numberPoints);

	}//for ambientDim

}//runTests


int main(int argc, char *argv[])
{
	//mainValuationDriver(argv, argc);
	runTests();
	//runOneTest(atoi(argv[1]), atoi(argv[2]));

	return 0;
}//main()
