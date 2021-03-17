/* ReadPolyhedron.cpp -- Handle command-line args to read a polyhedron

 Copyright 2007 Matthias Koeppe
 Derived from count.cpp, which is:
 Copyright 2002, 2003 Raymond Hemmecke, Ruriko Yoshida
 Copyright 2006 Matthias Koeppe

 This file is part of LattE.

 LattE is free software; you can redistribute it and/or modify it
 under the terms of the version 2 of the GNU General Public License
 as published by the Free Software Foundation.

 LattE is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with LattE; if not, write to the Free Software Foundation,
 Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
 */

#include <cstring>

#include "CheckEmpty.h"
#include "ReadLatteStyle.h"
#include "ReadPolyhedron.h"
#include "vertices/cdd.h"
#include "convert.h"
#include "preprocess.h"
#include "ramon.h"
#include "ReadSubcones.h"
#include "ProjectUp.h"
#ifdef HAVE_FORTYTWO_LIB
#  include "VertexConesWith4ti2.h"
#endif
#include "print.h"
#include "dual.h"

ReadPolyhedronData::ReadPolyhedronData() {
	strcpy(Vrepresentation, "no");
	strcpy(interior, "no");
	strcpy(dilation, "no");
	strcpy(dualApproach, "no");
	strcpy(nonneg, "no");
	strcpy(cddstyle, "no");
	strcpy(equationsPresent, "no");
	strcpy(Singlecone, "no");
	strcpy(grobner, "no");
	strcpy(maximum, "no");
	strcpy(minimize, "no");
	strcpy(assumeUnimodularCones, "no");
	strcpy(taylor, "no");
	strcpy(rationalCone, "no");
	strcpy(Memory_Save, "yes");

	vertexcones =
#ifdef HAVE_FORTYTWO_LIB
			VertexConesWith4ti2
#else
	VertexConesWithCdd
#endif
	;
	redundancycheck = FullRedundancyCheckWithCddlib;

	expect_dilation_factor = false;
	dilation_const = 1;
	expect_filename = true;
	degree = 1;

	input_homog_cone = false;
	input_vertex_cones = false;
	input_dualized = false;
	have_subcones = false;
	input_listcone_format = false;

	matrix = NULL;
	templistVec = NULL;
}

void ReadPolyhedronData::show_options(ostream &stream) {
	stream << "Standard input specifications:" << endl
			<< "  FILENAME                                 Inequalities in LattE format"
			<< endl
			<< "  --vrep FILENAME                          Vertices in LattE format"
			<< endl
			<< "  --cdd FILENAME.{ext,ine}                 Inequalities or vertices in CDD format"
			<< endl << "Input modifications:" << endl
			<< "  --dilation=DILATION-FACTOR               Dilate by DILATION-FACTOR"
			<< endl
			<< "  --interior                               Handle the interior of the polyhedron"
			<< endl
	//	 << "    +                                      - Add non-negativity constraints" << endl
			<< "Intermediate input specifications:" << endl
			<< "  --input-primal-homog-cone=CONE.ext       The homogenized polyhedron given by a "
			<< endl
			<< "                                           full-dimensional cone in CDD format"
			<< endl
			<< "  --input-dual-homog-cone=CONE.ext         The dual of the homogenized polyhedron given by a "
			<< endl
			<< "                                           full-dimensional cone in CDD format"
			<< endl
			<< "  --subcones=FILENAME                      Use a subdivision of the above specified"
			<< endl
			<< "                                           cone (up to lower-dimensional cones), given by "
			<< endl
			<< "                                           ray indicator vectors"
			<< endl
			<< "  --input-primal-homog-cones=CONES         The homogenized polyhedron given by a "
			<< endl
			<< "                                           union of cones (up to lower-dimensional cones) "
			<< endl
			<< "                                           in LattE's internal format"
			<< endl
			<< "  --input-dual-homog-cones=CONES           The dual of the homogenized polyhedron given by a "
			<< endl
			<< "                                           union of cones (up to lower-dimensional cones) "
			<< endl
			<< "                                           in LattE's internal format"
			<< endl
			<< "  --input-vertex-cones=CONES               The collection of vertex cones "
			<< endl
			<< "                                           in LattE's internal format"
			<< endl << "Input handling options:" << endl
			<< "  --compute-vertex-cones={cdd,lrs,4ti2}    Use this method for computing vertex cones"
			<< endl
			<< "  --redundancy-check={none,cddlib,full-cddlib}   Use this method for computing an irredundant "
			<< endl
			<< "                                           representation"
			<< endl << "Algorithmic option:" << endl
			<< "  --homog                                  Compute in homomogenized mode (by coning over the polytope) "
			<< endl
			<< "                                           rather than using the vertex cones"
			<< endl;
}

bool ReadPolyhedronData::parse_option(const char *arg) {
	/* Parse traditional LattE options. */
	if (strncmp(arg, "vrep", 3) == 0)
		strcpy(Vrepresentation, "yes");
	else if (strncmp(arg, "int", 3) == 0) {
		strcpy(interior, "yes");
		cerr
				<< "WARNING: Options `--interior' and `int' are broken for most methods."
				<< endl;
		THROW_LATTE(LattException::bug_Unknown,0);
	} else if (strncmp(arg, "homog", 3) == 0)
		strcpy(dualApproach, "yes");
	else if (strncmp(arg, "equ", 3) == 0) {
		cerr << "Warning: Ignoring the old-style LattE option `equ', "
				<< "since we detect the presence of equations ourselves."
				<< endl;
	} else if (strncmp(arg, "+", 1) == 0) {
		cerr << "Note: Recommend specifying nonnegativity constraints in the "
				<< endl
				<< "      input file rather than using the old-style LattE option `+'."
				<< endl;
		strcpy(nonneg, "yes");
	} else if (strncmp(arg, "cdd", 3) == 0)
		strcpy(cddstyle, "yes");
	else if (strncmp(arg, "dil", 3) == 0) {
		cerr << "Note: Old-style LattE option `dil FACTOR' corresponds to "
				<< endl << "      new option `--dilation=FACTOR'." << endl;
		strcpy(dilation, "yes");
		expect_dilation_factor = true;
	} else if (strncmp(arg, "lrs", 3) == 0) {
		cerr << "Note: Old-style LattE option `lrs' corresponds to " << endl
				<< "      new option `--compute-vertex-cones=lrs'." << endl;
		vertexcones = ReadPolyhedronData::VertexConesWithLrs;
	}
	/* Parse new options. */
	else if (strncmp(arg, "--dilation=", 11) == 0) {
		strcpy(dilation, "yes");
		dilation_const = atoi(arg + 11);
	} else if (strcmp(arg, "--interior") == 0) {
		cerr
				<< "WARNING: Options `--interior' and `int' are broken for most methods."
				<< endl;
		THROW_LATTE(LattException::bug_Unknown,0); //If we cannot stand by our computation, we should not let the user run it ~Brandon 2010. (I added the exit statement).
		strcpy(interior, "yes");
	} else if (strcmp(arg, "--vrep") == 0) {
		strcpy(Vrepresentation, "yes");
	} else if (strcmp(arg, "--homog") == 0) {
		strcpy(dualApproach, "yes");
	} else if (strcmp(arg, "--cdd") == 0) {
		strcpy(cddstyle, "yes");
	} else if (strncmp(arg, "--input-primal-homog-cone=", 26) == 0) {
		filename = arg + 26;
		expect_filename = false;
		input_homog_cone = true;
		input_dualized = false;
		strcpy(dualApproach, "yes");
	} else if (strncmp(arg, "--input-dual-homog-cone=", 24) == 0) {
		filename = arg + 24;
		expect_filename = false;
		input_homog_cone = true;
		input_dualized = true;
		strcpy(dualApproach, "yes");
	} else if (strncmp(arg, "--subcones=", 11) == 0) {
		subcones_filename = string(arg + 11);
		have_subcones = true;
	} else if (strncmp(arg, "--input-primal-homog-cones=", 27) == 0) {
		filename = arg + 27;
		expect_filename = false;
		input_homog_cone = true;
		input_dualized = false;
		input_listcone_format = true;
		strcpy(dualApproach, "yes");
	} else if (strncmp(arg, "--input-dual-homog-cones=", 25) == 0) {
		filename = arg + 25;
		expect_filename = false;
		input_homog_cone = true;
		input_dualized = true;
		input_listcone_format = true;
		strcpy(dualApproach, "yes");
	} else if (strncmp(arg, "--input-vertex-cones=", 21) == 0) {
		filename = arg + 21;
		expect_filename = false;
		input_vertex_cones = true;
		input_dualized = false;
		input_listcone_format = true;
	} else if (strncmp(arg, "--compute-vertex-cones=", 23) == 0) {
		if (strcmp(arg + 23, "cdd") == 0)
			vertexcones = VertexConesWithCdd;
		else if (strcmp(arg + 23, "lrs") == 0)
			vertexcones = VertexConesWithLrs;
		else if (strcmp(arg + 23, "4ti2") == 0)
			vertexcones = VertexConesWith4ti2;
		else {
			cerr << "Unknown vertex cone method: " << arg + 23 << endl;
			THROW_LATTE(LattException::ue_BadCommandLineOption, 0);
		}
	} else if (strncmp(arg, "--redundancy-check=", 19) == 0) {
		if (strcmp(arg + 19, "none") == 0)
			redundancycheck = NoRedundancyCheck;
		else if (strcmp(arg + 19, "cddlib") == 0)
			redundancycheck = RedundancyCheckWithCddlib;
		else if (strcmp(arg + 19, "full-cddlib") == 0)
			redundancycheck = FullRedundancyCheckWithCddlib;
		else {
			cerr << "Unknown redundancy check method: " << arg + 19 << endl;
			THROW_LATTE(LattException::ue_BadCommandLineOption, 0);
		}
	} else if (strncmp(arg, "--", 2) != 0) {
		// Regular argument, see if we expect one
		if (expect_dilation_factor) {
			dilation_const = atoi(arg);
			expect_dilation_factor = false;
		} else if (expect_filename) {
			filename = arg;
			expect_filename = false;
		} else
			return false;
	} else
		return false;
	return true;
}

/**
 * Converts a matrix to its v-rep or finds the tangent cones. The result is saved in the Polyhedron
 * @parm theMatirx: Full-dimensional inequality matrix in cdd format
 * @parm numberOfVars (not including RHS)
 * @parm Poly: the answer is returned here.
 * @parm parms: only used to record time.
 */
void  ReadPolyhedronData::matrixToVerticesOrCones(listVector * theMatrix, int numOfVars, Polyhedron *& Poly, BarvinokParameters *&params)
{
	if (dualApproach[0] == 'y') {
		Poly->numOfVars = numOfVars + 1;
		listVector *rays = NULL, *endRays, *tmpRays;
		Poly->cones = createListCone();
		Poly->cones->vertex = new Vertex(createRationalVector(numOfVars + 1));
		rays = createListVector(createVector(numOfVars + 1));
		endRays = rays;
		tmpRays = theMatrix;
		vec_ZZ v;
		v.SetLength(numOfVars + 1);
		while (tmpRays) {
			/* Change from CDD format ( b | -A ) to LattE's homogenized format ( A | -b ). */
			int i;
			for (i = 0; i < numOfVars; i++)
				v[i] = -(tmpRays->first)[i + 1];
			v[numOfVars] = -(tmpRays->first)[0];
			endRays->rest = createListVector(v);
			endRays = endRays->rest;
			tmpRays = tmpRays->rest;
		}
		Poly->cones->rays = rays->rest;
		delete rays; // deletes dummy head
		Poly->dualized = true;
		Poly->homogenized = true;
	} else {
		Poly->numOfVars = numOfVars;
		/* Compute vertices and edges. */
		listCone *tmpcones;

		params->vertices_time.start();
		switch (vertexcones) {
		case VertexConesWithCdd:
			tmpcones = computeVertexCones(filename.c_str(), theMatrix, numOfVars);
			break;
		case VertexConesWithLrs:
			tmpcones = computeVertexConesViaLrs(filename.c_str(), theMatrix,
					numOfVars);
			break;
		case VertexConesWith4ti2:
#ifdef HAVE_FORTYTWO_LIB
			tmpcones = computeVertexConesWith4ti2(theMatrix, numOfVars,
					Poly->unbounded);
#else
			cerr << "VertexConesWith4ti2 not compiled in, sorry" << endl;
			THROW_LATTE(LattException::ue_BadCommandLineOption,0);
#endif
			break;
		default:
			cerr << "Bad VertexConesType" << endl;
			abort();
		};

		Poly->cones = tmpcones;
		cerr << "The polytope has " << lengthListCone(Poly->cones)
				<< " vertices." << endl;
		//system_with_error_check("rm -f numOfLatticePoints");
		params->vertices_time.stop();
		cerr << params->vertices_time;
		Poly->homogenized = false;
	}


}//matrixToVerticesOrCones

static listCone *
read_cone_cdd_format(const string &filename) {
	FILE *in = fopen(filename.c_str(), "r");
	if (in == NULL) {
		cerr << "Unable to open CDD-style input file " << filename << endl;
		THROW_LATTE(LattException::fe_Open, 0);
	}
	dd_MatrixPtr M;
	dd_ErrorType err = dd_NoError;
	M = dd_PolyFile2Matrix(in, &err);
	if (err != dd_NoError) {
		cerr << "Parse error in CDD-style input file " << filename << endl;
		THROW_LATTE(LattException::fe_Parse, 0);
	}
	listCone *cone = cddlib_matrix_to_cone(M);
	dd_FreeMatrix(M);
	return cone;
}




/**
 * Returns the system Ax <= b where the polytope is full dimensional, starting from the file name. We read in the file and check the polytope is not empty for latte-files.
 * @parm params: I don't know if I will use this.
 * @return something. Caller is in charge of freeing the memory.
 */
listVector *ReadPolyhedronData::read_full_rank_inequality_matrix(BarvinokParameters *params)
{

	cout << "I think it is save to delete this function::ReadPolyhedronData::read_full_rank_inequality_matrix" << endl;
	THROW_LATTE(LattException::bug_Unknown);

	if (expect_filename) {
		cerr << "The input file name is missing." << endl;
		THROW_LATTE(LattException::ue_FileNameMissing,0);
	}

	dd_MatrixPtr M;

	if (cddstyle[0] == 'y') {
		/* Read an input file in CDD input format. */
		if (Vrepresentation[0] == 'y') {
			cerr
					<< "ReadPolyhedronData::read_full_rank_inequality_matrix:: Sorry, cannot compute projected H-rep starting from a V-rep.";
			THROW_LATTE(LattException::bug_NotImplementedHere, 0);
		}
		cerr << "Warning: Not performing check for empty polytope, "
				<< "because it is unimplemented for the CDD-style input format. "
				<< endl;
		M = ReadCddStyleMatrix(filename);
	} else {
		/* Read an input file in LattE format. */
		if (Vrepresentation[0] == 'y') {
			/* The polyhedron is given by its V-representation in a
			 LattE-style input format. */
			cerr << "ReadPolyhedronData::read_full_rank_inequality_matrix:: Sorry, cannot compute projected H-rep starting from a V-rep.";
			THROW_LATTE(LattException::bug_NotImplementedHere, 0);

		}

		/* Not VREP. */
		CheckEmpty(filename.c_str());
		M = ReadLatteStyleMatrix(filename.c_str(), /* vrep: */false,
		/* homogenize: */false,
		/* nonnegative: */nonneg[0] == 'y');
	}//if cdd or latte file format.

	/* Now we have in M the H-representation*/

	Polyhedron * Poly = new Polyhedron;
	int numOfVars = M->colsize - 1; /* Number of variables, not including RHS. */
	polyhedronRedundancyCheck(redundancycheck, M);
	matrix = projectOutVariables(M, numOfVars, Poly);
	dd_FreeMatrix(M);
	delete Poly; //not used right now.
	Poly = NULL;

	return matrix;
}//read_full_rank_inequality_matrix

Polyhedron *
ReadPolyhedronData::read_polyhedron(BarvinokParameters *params) {
	if (expect_filename) {
		cerr << "The input file name is missing." << endl;
		THROW_LATTE(LattException::ue_FileNameMissing, 0);
	}

	if (input_homog_cone)
		return read_polyhedron_from_homog_cone_input(params);
	else if (input_vertex_cones)
		return read_polyhedron_from_vertex_cone_input(params);
	else
		return read_polyhedron_hairy(params);
}


Polyhedron *ReadPolyhedronData::read_polyhedron(dd_MatrixPtr M, BarvinokParameters *params, const ReadPolyhedronOutput readPolyhedronOutput)
{
	Polyhedron * Poly;



	if ( set_card(M->linset) )
		strcpy(equationsPresent, "yes");
	else
		strcpy(equationsPresent, "no");


	if ( readPolyhedronOutput == ReadPolyhedronData::computeVertices)
	{
		strcpy(dualApproach, "yes");
	}
	else if ( readPolyhedronOutput == ReadPolyhedronData::computePrimalCones)
	{
		strcpy(dualApproach, "no");
	}

	switch (M->representation) {
	case dd_Generator: //given v-rep, want v-rep. done.
		if ( dualApproach[0] == 'y')
			Poly = PolyhedronFromVrepMatrix(M, true);
		else
			Poly = PolyhedronFromVrepMatrix(M, false);
		break;
	case dd_Inequality: //given h-rep, want v-rep.
		Poly = PolyhedronFromHrepMatrix(M, params);
		break;
	default:
		cerr << "Unknown representation" << endl;
		THROW_LATTE(LattException::bug_Unknown);
	}

	params->Number_of_Variables = Poly->numOfVars;

	/*
	cout << "poly before anything" << endl;
	Poly->printPolyhedron();
	cout << "poly after anything" << endl;
	cout << "dd h=" << ( M->representation == dd_Inequality ? 1 : 0) << endl;
	cout << "find cone" << (readPolyhedronOutput == ReadPolyhedronData::computePrimalCones ? 1 : 0) << endl;
	//find the primal rep.
	*/

	if ( Poly->dualized)
	{
		dualizeCones(Poly->cones, Poly->numOfVars, params);
		Poly->dualized = false;
	}
	if ( Poly->cones->rays == NULL)
	{
		dualizeCones(Poly->cones, Poly->numOfVars, params);
		dualizeCones(Poly->cones, Poly->numOfVars, params);
	}

	/*
	//if hrep to vertices
	if ( M->representation == dd_Inequality && readPolyhedronOutput == ReadPolyhedronData::computeVertices)
	{
		dualizeCones(Poly->cones, Poly->numOfVars, params);
		dualizeCones(Poly->cones, Poly->numOfVars, params);
		dualizeCones(Poly->cones, Poly->numOfVars, params);
		Poly->dualized = false;
	}
	else if (M->representation == dd_Inequality && readPolyhedronOutput == ReadPolyhedronData::computePrimalCones)
	{
		dualizeCones(Poly->cones, Poly->numOfVars, params);
		dualizeCones(Poly->cones, Poly->numOfVars, params);
	}
	*/

	return Poly;
}//read_polyhedron


Polyhedron *
ReadPolyhedronData::read_polyhedron_from_homog_cone_input(
		BarvinokParameters *params) {
	/* We are already given a full-dimensional, homogenized cone
	 or a list of those. */
	ConeProducer *producer = NULL;
	if (input_listcone_format) {
		if (have_subcones) {
			listCone *cones = readListConeFromFile(filename.c_str());
			if (lengthListCone(cones) != 1) {
				cerr
						<< "A subcones file can only be given for a single-cone file."
						<< endl;
				THROW_LATTE(LattException::bug_Unknown);
			}
			producer = new SubconeReadingConeProducer(cones, subcones_filename);
		} else {
			producer = new ListConeReadingConeProducer(filename);
		}
	} else {
		listCone *cone = read_cone_cdd_format(filename);
		if (have_subcones) {
			// Also a subcones file given.
			producer = new SubconeReadingConeProducer(cone, subcones_filename);
		} else {
			producer = new SingletonConeProducer(copyCone(cone));
		}
	}
	/* Use the producer to create the polyhedron. */
	CollectingConeConsumer ccc;
	producer->Produce(ccc);
	delete producer;
	Polyhedron *Poly = new Polyhedron;
	Poly->cones = ccc.Collected_Cones;
	int numOfVars;
	if (Poly->cones == NULL || Poly->cones->rays == NULL)
		numOfVars = 0;
	else
		numOfVars = Poly->cones->rays->first.length();
	Poly->numOfVars = numOfVars;
	Poly->homogenized = true;
	Poly->dualized = input_dualized;
	return Poly;
}

Polyhedron *
ReadPolyhedronData::read_polyhedron_from_vertex_cone_input(
		BarvinokParameters *params) {
	ConeProducer *producer;
	producer = new ListConeReadingConeProducer(filename);
	CollectingConeConsumer ccc;
	producer->Produce(ccc);
	delete producer;
	Polyhedron *Poly = new Polyhedron;
	Poly->cones = ccc.Collected_Cones;
	int numOfVars;
	if (Poly->cones == NULL)
		numOfVars = 0;
	else
		numOfVars = ambient_cone_dimension(Poly->cones);
	//printListCone(Poly->cones, numOfVars);
	Poly->numOfVars = numOfVars;
	Poly->homogenized = false;
	Poly->dualized = input_dualized;
	return Poly;
}

static dd_MatrixPtr ReadCddStyleMatrix(const string &filename) {
	FILE *in = fopen(filename.c_str(), "r");
	if (in == NULL) {
		cerr << "Unable to open CDD-style input file " << filename << endl;
		THROW_LATTE(LattException::fe_Open, 0);
	}
	dd_MatrixPtr M;
	dd_ErrorType err = dd_NoError;
	M = dd_PolyFile2Matrix(in, &err);
	if (err != dd_NoError) {
		cerr << "Parse error in CDD-style input file " << filename << endl;
		THROW_LATTE(LattException::fe_Parse, 0);
	}
	return M;
}

Polyhedron *
ReadPolyhedronData::read_polyhedron_hairy(BarvinokParameters *params) {
	Polyhedron *Poly = NULL;

	if (expect_filename) {
		cerr << "The input file name is missing." << endl;
		THROW_LATTE(LattException::ue_FileNameMissing, 0);
	}

	dd_MatrixPtr M;

	if (cddstyle[0] == 'y') {
		/* Read an input file in CDD input format. */
		if (Vrepresentation[0] == 'y') {
			cerr
					<< "The command-line keyword `vrep' denotes the use of a LattE-style "
					<< endl
					<< "input format giving the V-representation.  If you want to give "
					<< endl
					<< "the a V-representation in CDD format, just do that, but don't use "
					<< endl << "the `vrep' keyword." << endl;
			THROW_LATTE(LattException::ue_BadCommandLineOption, 0);
		}
		cerr << "Warning: Not performing check for empty polytope, "
				<< "because it is unimplemented for the CDD-style input format. "
				<< endl;
		M = ReadCddStyleMatrix(filename);
	} else {
		/* Read an input file in LattE format. */
		if (Vrepresentation[0] == 'y') {
			/* The polyhedron is given by its V-representation in a
			 LattE-style input format. */
			if (dilation_const != 1) {
				cerr << "Dilation unimplemented for `vrep' input" << endl;
				THROW_LATTE(LattException::ue_BadCommandLineOption, 0);
			}
			if (dualApproach[0] != 'y') {
				/* FIXME: Special case that ought to be handled uniformly.
				 *
				 * Brandon: I think I fixed this by changing PolyhedronFromVrepMatrix.
				 *          I think we can just delete this if-statement, and it will work.
				 */
				/* Don't homogenize. */
				Polyhedron *P = new Polyhedron;
				P->cones = computeVertexConesFromVrep(filename.c_str(),
						P->numOfVars);

				P->dualized = false;
				P->homogenized = false;
				return P; /* Directly deliver the polyhedron. */
			}
			M = ReadLatteStyleMatrix(filename.c_str(), /* vrep: */true,
			/* homogenize: */false);
		} else {
			/* Not VREP. */
			CheckEmpty(filename.c_str());
			M = ReadLatteStyleMatrix(filename.c_str(), /* vrep: */false,
			/* homogenize: */false,
			/* nonnegative: */nonneg[0] == 'y');
		}
	}

	/* Now we have in M the H-representation or the V-representation. */

	switch (M->representation) {
	case dd_Generator: /* V-representation */
		return PolyhedronFromVrepMatrix(M, /* homogenize: */dualApproach[0]
				== 'y');
	case dd_Inequality: /* H-representation */
		return PolyhedronFromHrepMatrix(M, params);
	default:
		cerr << "Unknown representation" << endl;
		THROW_LATTE(LattException::bug_Unknown, 0);
	}
}

Polyhedron *
ReadPolyhedronData::PolyhedronFromHrepMatrix(dd_MatrixPtr M,
		BarvinokParameters *params) {
	Polyhedron *Poly = new Polyhedron;
	int numOfVars = M->colsize - 1; /* Number of variables, not
	 including RHS. */

	params->read_time.start();


	polyhedronRedundancyCheck(redundancycheck, M);

	matrix = projectOutVariables(M, numOfVars, Poly);
	dd_FreeMatrix(M);


	//matrix = matrixTmp;
	params->read_time.stop();
	cerr << params->read_time;
	/* Now matrix contains the new inequalities. */


	matrixToVerticesOrCones(matrix, numOfVars, Poly, params);

	return Poly;
}

/**
 * Finds hidden equalities and inequalities.
 * @parm redunType: which method to use
 * @parm M: Note, we take M by reference, and it is updated.
 */
void ReadPolyhedronData::polyhedronRedundancyCheck(RedundancyCheckType redunType, dd_MatrixPtr &M)
{

	switch (redunType) {
	case NoRedundancyCheck:
		break;
	case RedundancyCheckWithCddlib: {
		cerr << "Finding hidden equalities using cddlib...";
		cerr.flush();
		dd_rowset impl_lin;
		dd_rowindex newpos;
		dd_ErrorType err;
		dd_MatrixCanonicalizeLinearity(&M, &impl_lin, &newpos, &err);
		check_cddlib_error(err, "PolyhedronFromHrepMatrix");
		cerr << "done. " << endl;
		break;
	}
	case FullRedundancyCheckWithCddlib: {
		cerr
				<< "Removing redundant inequalities and finding hidden equalities using cddlib...";
		cerr.flush();
		dd_rowset impl_lin, redset;
		dd_rowindex newpos;
		dd_ErrorType err;
		dd_MatrixCanonicalize(&M, &impl_lin, &redset, &newpos, &err);
		check_cddlib_error(err, "polyhedronRedundancyCheck");
		cerr << "done. " << endl;
		break;
	}
	default:
		cerr << "Unknown redundancy check" << endl;
		abort();
	}
}//polyhedronRedundancyCheck

/**
 * reduces the input matrix to a full-dimensional matrix.
 * @parm numOfVars: Note that preprocessProblem() does change this value.
 * @parm Poly: the projecting_up_transducer is saved here.
 * @parm M: the H-rep matrix after removing redundant rows and adding hiden equalities.
 */
listVector * ReadPolyhedronData::projectOutVariables(dd_MatrixPtr &M, int &numOfVars, Polyhedron *& Poly)
{
	listVector *equations, *inequalities;
	cddlib_matrix_to_equations_and_inequalities(M, &equations, &inequalities);


	cerr << "Ax <= b, given as (b|-A):\n";
	cerr << "=========================\n";
	printListVectorToFile(cerr, inequalities, numOfVars + 1);
	cerr << endl;
	cerr << "Ax = b, given as (b|-A):\n";
	cerr << "========================\n";
	printListVectorToFile(cerr, equations, numOfVars + 1);
	cerr << endl;

	if (equations != NULL)
		strcpy(equationsPresent, "yes");
	else
		strcpy(equationsPresent, "no");

	/* Project out variables using equations. */
	mat_ZZ ProjU, ProjU2;
	ProjU.SetDims(numOfVars, numOfVars);
	ProjU2.SetDims(numOfVars, numOfVars);
	oldnumofvars = numOfVars;

	listVector *matrixTmp;
	if (equationsPresent[0] == 'y') {
		{
			vec_ZZ *generators = NULL;
			matrixTmp = preprocessProblem(equations, inequalities, &generators,
					&numOfVars, cost, ProjU, interior, dilation_const);
			if (generators)
				delete[] generators;
		}
		freeListVector(equations);
		freeListVector(inequalities);
		ProjU2 = transpose(ProjU);
		bb = ProjU2[0];
		mat_ZZ AAA;
		AAA.SetDims(ProjU2.NumRows() - 1, ProjU2.NumCols());
		int i;
		for (i = 1; i <= numOfVars; i++) {
			AAA[i - 1] = ProjU2[i];
		}
		AA = transpose(AAA);
		// cerr << ProjU << determinant(transpose(AAA)*AAA) <<  endl;
		templistVec = transformArrayBigVectorToListVector(ProjU,
				ProjU.NumCols(), ProjU.NumRows()); //what is this for?
		Poly->projecting_up_transducer = new ProjectingUpConeTransducer(
				oldnumofvars, numOfVars, AA, bb);
	} else {
		/* No equations. */
		dilateListVector(inequalities, numOfVars, dilation_const);
		matrixTmp = inequalities;
	}
	return matrixTmp;
}//project_out_variables


/**
 * Finds the transformation that reduces the input matrix to a full-dimensional matrix,
 *   but DOES NOT perform the transformation.
 * @parm numOfVars: Note that preprocessProblem() does NOT change this value.
 * @parm Poly: the projecting_up_transducer is saved here. I don't know what this is.
 * @parm M: the H-rep matrix after removing redundant rows and adding hidden equalities.
 * @return mattrix of the lattice basis. Each row is a basis element.
 */
 mat_ZZ ReadPolyhedronData::findLatticeBasis(dd_MatrixPtr &M, int &numOfVars)
{
	listVector *equations, *inequalities;
	cddlib_matrix_to_equations_and_inequalities(M, &equations, &inequalities);


	cerr << "Ax <= b, given as (b|-A):\n";
	cerr << "=========================\n";
	printListVectorToFile(cerr, inequalities, numOfVars + 1);
	cerr << endl;
	cerr << "Ax = b, given as (b|-A):\n";
	cerr << "========================\n";
	printListVectorToFile(cerr, equations, numOfVars + 1);
	cerr << endl;

	if (equations != NULL)
		strcpy(equationsPresent, "yes");
	else
		strcpy(equationsPresent, "no");

	/* Project out variables using equations. */
	mat_ZZ ProjU, ProjU2;
	ProjU.SetDims(numOfVars, numOfVars);
	ProjU2.SetDims(numOfVars, numOfVars);
	oldnumofvars = numOfVars;

	listVector *matrixTmp;
	vec_ZZ *generators = NULL;

	int numGenerators = numOfVars;
	if (equationsPresent[0] == 'y') {
		//reuse the preprocessProblem function, but only keep the lattice basis.

		matrixTmp = preprocessProblem_hack(equations, inequalities, &generators,
				&numGenerators, cost, ProjU, interior, dilation_const);
		//freeListVector(matrixTmp);

		/*
		ProjU2 = transpose(ProjU);
		bb = ProjU2[0];
		mat_ZZ AAA;
		AAA.SetDims(ProjU2.NumRows() - 1, ProjU2.NumCols());
		int i;

		for (i = 1; i < ProjU2.NumRows(); i++) {
			AAA[i - 1] = ProjU2[i];
		}

		AA = transpose(AAA);
		// cerr << ProjU << determinant(transpose(AAA)*AAA) <<  endl;
		templistVec = transformArrayBigVectorToListVector(ProjU,
				ProjU.NumCols(), ProjU.NumRows()); //what is this for?

		Poly->projecting_up_transducer = new ProjectingUpConeTransducer(
				oldnumofvars, numOfVars, AA, bb);
		*/
	}else
	{
		cout << "ReadPolyhedronData::findLatticeBasis: should only be called when the polytope has equations, error." << endl;
		THROW_LATTE(LattException::pe_UnexpectedRepresentation, 0);
	}//else {
		/* No equations. */
		//dilateListVector(inequalities, numOfVars, dilation_const);
		//matrixTmp = inequalities;
	//}
	//freeListVector(equations);
	//freeListVector(inequalities);

	assert(generators[0].length() == numOfVars);


	mat_ZZ basis;
	basis.kill();
	basis.SetDims(numOfVars, numGenerators);

	cout << "print the generators" << numGenerators << endl;
	for ( int i = 0; i < numGenerators; ++i)
	{
		cout << "i=" << i << " ";
		for (int j = 0; j < numOfVars; ++j)
		{
			cout << generators[i][j] << ", ";
			basis[j][i] = generators[i][j];
		}
		cout << endl;
	}

	delete [] generators;
	return basis;
}//findLatticeBasis

 /**
  * @parm matrix: vertex matrix.
  * @parm homogenize: true = we want to compute the vertices of the polytope.
  * 				  false = we want to compute the tangent cones.
  */
Polyhedron *PolyhedronFromVrepMatrix(dd_MatrixPtr matrix, bool homogenize) {
	Polyhedron *P = new Polyhedron;

	if (homogenize) {
		/* Homogenize. */
		dd_ErrorType error;
		dd_rowset redundant = dd_RedundantRows(matrix, &error);
		check_cddlib_error(error, "ReadLatteStyleVrep");
		/* The non-redundant rows are the rays of the homogenization. */
		int i;
		listCone *cone = createListCone();
		P->numOfVars = matrix->colsize;
		vec_ZZ ray;
		ray.SetLength(matrix->colsize);
		for (i = 1; i <= matrix->rowsize; i++) {
			if (!set_member(i, redundant)) {
				int j;
				/* CDD has homogenization in the 0-th,
				 LattE expects it in the last coordinate. */
				for (j = 0; j < matrix->colsize - 1; j++)
					ray[j] = convert_mpq_to_ZZ(matrix->matrix[i - 1][j + 1]);
				ray[matrix->colsize - 1] = convert_mpq_to_ZZ(matrix->matrix[i - 1][0]);
				cone->rays = appendVectorToListVector(ray, cone->rays);
				cone->vertex = new Vertex(createRationalVector(P->numOfVars));
			}
		}
		dd_FreeMatrix(matrix);
		P->cones = cone;
		P->dualized = false;
		P->homogenized = true;
	} else {
		/* Don't homogenize: compute the tangetn cones. */

		P->cones = computeVertexConesFromVrep(matrix, P->numOfVars);
		P->dualized = false;
		P->homogenized = false;
	}
	return P;
}

