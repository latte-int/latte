/*
 * Driver.cpp
 *
 *  Created on: Jun 24, 2010
 *      Author: Brandon Dutra and Gregory Pinto
 *
 *
 *
 *  type --help to print the help menu.
 */




#include <cstdlib>
#include <iostream>
#include <string>
#include <cassert>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <vector>







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
#include "testEhrhart/BuildHypersimplexEdgePolytope.h"

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
#include "preprocess.h"
#include "print.h"
#include "ramon.h"
#include "rational.h"
#include "timing.h"
#include "flags.h"
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




typedef struct
{
	RationalNTL triangulate;
	RationalNTL lawrence;

} VolumesContainer;


VolumesContainer computeVolume(listCone * cones, BarvinokParameters &myParameters,
		const char *valuationType);
VolumesContainer mainValuationDriver(char *argv[], int argc);
void runOneTest(int ambientDim, int numPoints);
void runTests();

//BarvinokParameters parameters;
//ReadPolyhedronData read_polyhedron_data;
//string output_filename;


/* ----------------------------------------------------------------- */

/**
 * Computes the volume of the polytope.
 */
VolumesContainer computeVolume(listCone * cones, BarvinokParameters &myParameters,
		const char *valuationAlg, const char * print)
{
	VolumesContainer ans;
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


	ans.lawrence = ans2;
	ans.triangulate = ans1;

	return ans;

}//computeVolume

VolumesContainer mainValuationDriver(char *argv[], int argc)
{

	set_program_name(argv[0]);

	BarvinokParameters params;
	VolumesContainer volumeAns;
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


	if ( argc == 1)
	{
		cout << "type --help to print a help menu." << endl;
		exit(1);
	}//if too few parameters.
	if ( argc == 2 && !strcmp(argv[1], "--help"))
	{
		cout << "usage: valuation-exe [valuation  type] [valuation algorithm] <latte file>\n";
		cout << "valuation types: volume\n";
		cout << "  volume algorithm: [--lawrence  <--printLawrenceFunction> | --triangulate | --all]\n";
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
		if (strncmp(argv[i], "--lawrence", 8) == 0)
			strcpy(valuationAlg, "--lawrence");
		if (strcmp(argv[i], "--triangulate") == 0)
			strcpy(valuationAlg, "triangulate");
		if (strcmp(argv[i], "--all") == 0)
			strcpy(valuationAlg, "all");
		if ( strcmp(argv[i], "--printLawrenceFunction") == 0)
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
					volumeAns = computeVolume(cones, myParameters, valuationAlg, printLawrenceFunction);
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

	return (volumeAns);
}//mainValuationDrivver()




void runOneTest(int ambientDim, int numPoints)
{
	char * argv[] = {"runTests()", "--all", 0};
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


void runHyperSimplexTests()
{
	char * argv[] = {"runHyperSimplexTests()", "--all", 0};
						//   n  k  num/denom
	int  hyperSimplexData[][4] = { /*{4, 1, 1, 6},
							{4, 2, 2, 3},
							{5, 1, 1, 24},
							{5, 2, 11, 24},
							{6, 1, 1, 120},
							{6, 2, 13, 60},
							{6, 3, 11, 20},
							{7, 1, 1, 720},
							{7, 2, 19, 240},
							{7, 3, 151, 360},
							{8, 1, 1, 5040},
							{8, 2, 1, 42},
							{8, 3, 397, 1680},
							{8, 4, 151, 315},
							{9, 1, 1, 40320},
							{9, 2, 247, 40320},
							{9, 3, 477, 4480},
							{9, 4, 15619, 40320},
							{10, 1, 1, 362880},
							{10, 2, 251, 181440},
							{10, 3, 913, 22680},*/
							{10, 4, 44117, 18440},
							{10, 5, 15619, 36288},
							{11, 1, 1, 3628800},
							{11, 2, 1013, 3628800},
							{11, 3, 299, 22680},
							{11, 4, 56899, 45300},
							{11, 5, 655177, 1814400},
							{12, 1, 1, 39916800},
							{12, 2, 509, 9979200},
							{12, 3, 50879, 13305600},
							{12, 4, 1093, 19800},
							{12, 5, 1623019, 6652800},
							{12, 6, 655177, 1663200}
	};//hyperSimplexData

	int numberTestCases = 34;

	VolumesContainer volumeAnswer;
	for(int i = 0; i < numberTestCases; ++i)
	{
		stringstream comments;
		BuildHypersimplexEdgePolytope hyperSimplex(hyperSimplexData[i][0], hyperSimplexData[i][1]);

		comments << "finding volume of Hypersimplex(" << hyperSimplexData[i][0] << ", " << hyperSimplexData[i][1] << ")";
		hyperSimplex.buildPolymakeFile();



		hyperSimplex.setComments(comments.str().c_str());

		hyperSimplex.buildPolymakeFile(); //make the file
		hyperSimplex.callPolymake(); //run polymake
		//hyperSimplex.findVolumeWithPolymake(); //run polymake for the volume
		hyperSimplex.convertFacetEquations(); //fix facet equations
		hyperSimplex.printFacetEquationsForLattE(); //make latte file.

		string file = hyperSimplex.getLatteFile();

		cout << comments.str().c_str() << endl;

		char * sFile = new char[file.size() + 1];
		strcpy(sFile, file.c_str());
		argv[2] = sFile;
		volumeAnswer = mainValuationDriver(argv, 3);
		delete [] sFile;

		RationalNTL correctVolumeAnswer(hyperSimplexData[i][2], hyperSimplexData[i][3]);
		if ( correctVolumeAnswer != volumeAnswer.lawrence
				|| correctVolumeAnswer != volumeAnswer.triangulate)
		{
			cout << "******* ERROR ******" << endl;
			cout << "correct answer: " << correctVolumeAnswer << endl;
			cout << "lawrence: " << volumeAnswer.lawrence << endl;
			cout << "triangulate: " << volumeAnswer.triangulate << endl;
			cout << "see file " << file.c_str() << endl;
			exit(1); //dont' delete the latte file.
		}//if error
		else
			cout << comments.str().c_str() << " CORRECT!" << endl;


	}//for i.



}//runHyperSimplexTests


/**
 * Get list of files in directory dir.
 */
int getdir (string dir, vector<string> &files)
{
    DIR *dp;
    struct dirent *dirp;
    if((dp  = opendir(dir.c_str())) == NULL) {
        cout << "Error(" << errno << ") opening " << dir << endl;
        return errno;
    }

    while ((dirp = readdir(dp)) != NULL) {
        files.push_back(string(dirp->d_name));
    }
    closedir(dp);
    return 0;
}


void runBirkhoffTests()
{
    string dir = string("../../EXAMPLES/birkhoff/");
    vector<string> files = vector<string>();
	VolumesContainer volumeAnswer;
	char * argv[] = {"runHyperSimplexTests()", "--all", 0};

    getdir(dir,files);

    for(int i = 0; i < file.size(); ++i)
    {
    	if ( strrncmp(file[i], ".latte", 6) == 0)
    	{


    		char * sFile = new char[file[i].length() + 1];
    		strcpy(sFile, file[i].c_str());
    		argv[2] = sFile;
    		volumeAnswer = mainValuationDriver(argv, 3);
    		delete [] sFile;


    	}//if a latte file.
    }//for ever file in the directory



}//runBirkhoffTests


int main(int argc, char *argv[])
{
	//mainValuationDriver(argv, argc);
	runHyperSimplexTests();
	//runBirkhoffTests();
	//runTests();
	//runOneTest(atoi(argv[1]), atoi(argv[2]));

	return 0;
}//main()
