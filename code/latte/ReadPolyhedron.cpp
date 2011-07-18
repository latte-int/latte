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
		exit(1); //If we cannot stand by our computation, we should not let the user run it ~Brandon 2010. (I added the exit statement).
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
			exit(1);
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
			exit(1);
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
			exit(1);
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
		exit(1);
	}
	dd_MatrixPtr M;
	dd_ErrorType err = dd_NoError;
	M = dd_PolyFile2Matrix(in, &err);
	if (err != dd_NoError) {
		cerr << "Parse error in CDD-style input file " << filename << endl;
		exit(1);
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
	exit(1);

	if (expect_filename) {
		cerr << "The input file name is missing." << endl;
		exit(2);
	}

	dd_MatrixPtr M;

	if (cddstyle[0] == 'y') {
		/* Read an input file in CDD input format. */
		if (Vrepresentation[0] == 'y') {
			cerr
					<< "ReadPolyhedronData::read_full_rank_inequality_matrix:: Sorry, cannot compute projected H-rep starting from a V-rep.";
			exit(2);
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
			exit(2);

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
		exit(2);
	}

	if (input_homog_cone)
		return read_polyhedron_from_homog_cone_input(params);
	else if (input_vertex_cones)
		return read_polyhedron_from_vertex_cone_input(params);
	else
		return read_polyhedron_hairy(params);
}

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
				exit(1);
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
		exit(1);
	}
	dd_MatrixPtr M;
	dd_ErrorType err = dd_NoError;
	M = dd_PolyFile2Matrix(in, &err);
	if (err != dd_NoError) {
		cerr << "Parse error in CDD-style input file " << filename << endl;
		exit(1);
	}
	return M;
}

Polyhedron *
ReadPolyhedronData::read_polyhedron_hairy(BarvinokParameters *params) {
	Polyhedron *Poly = NULL;

	if (expect_filename) {
		cerr << "The input file name is missing." << endl;
		exit(2);
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
			exit(2);
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
				exit(1);
			}
			if (dualApproach[0] != 'y') {
				/* FIXME: Special case that ought to be handled uniformly. */
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
		abort();
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
		exit(1);
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
				ray[matrix->colsize - 1] = convert_mpq_to_ZZ(matrix->matrix[i
						- 1][0]);
				cone->rays = appendVectorToListVector(ray, cone->rays);
				cone->vertex = new Vertex(createRationalVector(P->numOfVars));
			}
		}
		dd_FreeMatrix(matrix);
		P->cones = cone;
		P->dualized = false;
		P->homogenized = true;
	} else {
		/* Don't homogenize. */
		cerr << "PolyhedronFromVrepMatrix: Unimplemented for homogenize=false."
				<< endl;
		abort();
	}
	return P;
}

// *********************************************************************
// ********* Start of the  ReadPolyhedronDataRecursive functions *******
// *********************************************************************

/**
 * Constructor.
 * @parm rpd: We should copy the input options. The vertex-cone and redundancy
 *   check are the two most important options. I'm not sure if we need the others,
 *   but lets just copy them too. After we get everyhing working, we then can comment things out
 *   and see what we really need.
 */
ReadPolyhedronDataRecursive::ReadPolyhedronDataRecursive(const ReadPolyhedronData & rpd)
{

	// copy the maze of twisty parameters.
	strcpy( equationsPresent, rpd.equationsPresent);
	strcpy( nonneg, rpd.nonneg);
	strcpy( cddstyle, rpd.cddstyle);
	strcpy( Vrepresentation, rpd.Vrepresentation);
	strcpy( dilation, rpd.dilation);
	strcpy( interior, rpd.interior);
	dilation_const = rpd.dilation_const;
	strcpy( dualApproach, rpd.dualApproach);
	filename = rpd.filename;
	strcpy( Memory_Save, rpd.Memory_Save);
	strcpy( grobner, rpd.grobner);
	strcpy( maximum, rpd.maximum);
	strcpy( minimize, rpd.minimize);
	strcpy( taylor, rpd.taylor);
	strcpy( rationalCone, rpd.rationalCone);
	strcpy( assumeUnimodularCones, rpd.assumeUnimodularCones);
	strcpy( Singlecone, rpd.Singlecone);
	degree = rpd.degree;

	// copy data for input of cones
	input_homog_cone = rpd.input_homog_cone;
	input_vertex_cones = rpd.input_vertex_cones;
	input_dualized = rpd.input_dualized;
	have_subcones = rpd.have_subcones;
	input_listcone_format = rpd.input_listcone_format;
	subcones_filename = rpd.subcones_filename;


	//These two are the most important to copy over.
	vertexcones = rpd.vertexcones;
	redundancycheck = rpd.redundancycheck;

	// I don't think we need to copy these.
	cost = rpd.cost;
	matrix = NULL;  // don't copy.
	//AA = rpd.AA;
	//bb = rpd.bb;
	oldnumofvars = rpd.oldnumofvars;
	templistVec = NULL; //don't copy.

	//I don't think we need these.
	expect_dilation_factor = rpd.expect_dilation_factor;
	expect_filename = rpd.expect_filename;

}//constructor.


ReadPolyhedronDataRecursive::~ReadPolyhedronDataRecursive()
{
	cout << "ReadPoly data rec deconstructor start" << endl;
	freeListVector(templistVec);
	freeListVector(matrix);
}


void ReadPolyhedronDataRecursive::dilatePolytope()
{
	Polyhedron *Poly = findTangentCones();
/*
	cout << "rpdr::dd matrix before dilations" << endl;
	for(int i = 0; i < ddHrep->rowsize; ++i)
	{
		for(int j = 0; j < ddHrep->colsize; ++j)
			cout << ddHrep->matrix[i][j] << " ";
		cout << endl;
	}
	cout << "end rpdr::dd matrix before dilation" << endl;
*/
	//now find dilation factor.
	dilationNum = 1;
	for(listCone * cone = Poly->cones; cone; cone = cone->rest)
	{
		for(int i = 0; i < Poly->numOfVars; ++i) //find gcd of b in Ax <= b
			dilationNum = (dilationNum *  (cone->vertex->vertex->denominators())[i])/ GCD(dilationNum, (cone->vertex->vertex->denominators())[i]);
	}

	if (dilationNum == 1)
		return;

	//update the hrep matrix and the list cones.

	//dd_mul(x, y, z)   Set x to be the multiplication of y and z.
	mpq_class x;
	x = convert_ZZ_to_mpz(dilationNum);
	for(int i = 0; i < ddHrep->rowsize; ++i)
		dd_mul(ddHrep->matrix[i][0], ddHrep->matrix[i][0], x.get_mpq_t() );

	//update the tangent cones
	for (listCone * cone = Poly->cones; cone; cone = cone->rest)
	{
		cone->vertex->vertex->scalarMultiplication(dilationNum, to_ZZ(1));
	}//for every vertex.

	//cout << "dilation factor for rpdr dilate is " << dilationNum << endl;
	//cout << "rpdr::dilate" << endl;
	//printListCone(Poly->cones, Poly->numOfVars);
	//cout << "end rpdr::dilate" << endl;

	/*
	cout << "rpdr::dd matrix" << endl;
	for(int i = 0; i < ddHrep->rowsize; ++i)
	{
		for(int j = 0; j < ddHrep->colsize; ++j)
			cout << ddHrep->matrix[i][j] << " ";
		cout << endl;
	}
	cout << "end rpdr::dd matrix" << endl;
	*/
	delete Poly;
}//dilatePolytope();

int ReadPolyhedronDataRecursive::dimension()
{
	return (matrix->first).length();
}

int ReadPolyhedronDataRecursive::getFullDimensionCount() const
{
	//cout << "ddhrep colsie=" <<ddHrep->colsize << endl;
	//cout << "total size=" << set_card(ddHrep->linset) << endl;
	//for(int i = 1; i <= ddHrep->rowsize; ++i)
	//	cout << "i=" << i << ":" << set_member(i,ddHrep->linset) << endl;

	return (ddHrep->colsize -1- getNumberEqualities());
}

Polyhedron * ReadPolyhedronDataRecursive::findTangentCones()
{
	assert(dualApproach[0] != 'y');
	Polyhedron *Poly = new Polyhedron;

	Poly->numOfVars = ddHrep->colsize -1;

	/* Compute vertices and edges. */
	listCone *tmpcones;
	tmpcones = computeVertexCones(filename.c_str(), ddHrep);
	Poly->cones = tmpcones;
	cerr << "The polytope has " << lengthListCone(Poly->cones)
		<< " vertices." << endl;
	Poly->homogenized = false;

	//cout << " ReadPolyhedronDataRecursive::findTangentCones printing the tangent cones" << endl;
	//printListCone(Poly->cones, Poly->numOfVars);
	//cout << "end printing tangent cones";

	return Poly;
}//findTangentCones

void ReadPolyhedronDataRecursive::getFacetPolytope(int row, ReadPolyhedronDataRecursive &newMatrix, vec_ZZ & l, RationalNTL &lDotNormal)
{
	if ( set_member(row, ddHrep->linset) == true)
	{
		lDotNormal = 0;
		return;
	}//if this row is already an equality.

	//find <l, n_row>
	//ddHrep is in the form [b -A].
	//convert l to gmp vector.
	vector<mpq_class> L;
	L.resize(l.length());
	for (int i = 0; i < l.length(); ++i) L[i] = convert_ZZ_to_mpq(l[i]);

	mpq_class LDotNormal;
	LDotNormal = 0;
	assert(l.length() == ddHrep->colsize -1);
	mytype dotProdSum;
	dd_init(dotProdSum);
	cout << "get facet:: going to do dot prod, row " << row << " m.row=" << ddHrep->rowsize << endl;
	/*
	for(int i = 0; i < ddHrep->rowsize; ++i)
	{
		for(int j = 0; j < ddHrep->colsize; ++j)
			cout << ddHrep->matrix[i][j] << ",";
		cout << endl;
	}
	*/
	for(int i = 0; i < l.length(); ++i)
	{
		cout << ddHrep->matrix[row-1][i+1] << " : " << l[i] << endl;
	}
	for(int i = 0; i < l.length(); ++i)
	{
		mytype x;
		dd_init(x);
		dd_mul(x, L[i].get_mpq_t(), (ddHrep->matrix[row-1][i+1]));
		dd_add(dotProdSum, dotProdSum, x);
		//LDotNormal += L[i] * (ddHrep->matrix[row][i+1]);
	}
	LDotNormal = mpq_class(dotProdSum);
	LDotNormal *= (-1);
	lDotNormal = RationalNTL(convert_mpz_to_ZZ(LDotNormal.get_num()), convert_mpz_to_ZZ(LDotNormal.get_den())); //convert to rationalNTL.
	cout << "dot produce worked ok =" << lDotNormal << endl;
	if (LDotNormal == 0)
	{
		assert(lDotNormal.getNumerator() == 0);
		//return;
	}//<y, n> = 0, so stop.

	//now find the new polytope.
	newMatrix.setMatrix(dd_CopyMatrix(ddHrep));
	newMatrix.setInequality(row);
	newMatrix.readHrepMatrix();

	//now do redundancy on the new matrix.


}//getFacetPolytope

void ReadPolyhedronDataRecursive::readHrepMatrix()
{
	/*start of read PolyhedronFromHrepMatrix */
	int numOfVars = ddHrep->colsize - 1; /* Number of variables, not including RHS. */


	polyhedronRedundancyCheck(redundancycheck, ddHrep);

	if ( set_card(ddHrep->linset) > 0)
	{
		latticeBasis = findLatticeBasis(ddHrep, numOfVars);
	}


	/* Now matrix contains the new inequalities. */


	//ok, ddHrep now contains the polytpoe, and latticeBasis contains the basis for the lattice affine space of the polytope if there are equations.


	//dilatePolytope(); not needed?
}


//assumes latticeInverse() has already been called.
const mat_ZZ * ReadPolyhedronDataRecursive::getLatticeInverse() const
{
	return &latticeLeftInverse;
}

const ZZ * ReadPolyhedronDataRecursive::getLatticeInverseDilation() const
{
	return &latticeLeftInverseDilation;
}


RationalNTL ReadPolyhedronDataRecursive::getNormalFactor() const
{
	mat_ZZ mat;
	mat.SetDims(ddHrep->colsize -1, ddHrep->colsize -1);

	for(int i = 0; i < mat.NumRows(); ++i)
	{
		for(int j = 1; j < mat.NumCols(); ++i)
			mat[i][j] = latticeBasis[i][j-1];

cout << "ReadPolyhedronDataRecursive::getNormalFactor(). start here, what is matrix col." << endl;
exit(1);
//		assert(ddHrep->matrix[i][])
	}
//start here --uncomment
//	mpz_class nDilation
//	dilationNum = 1;
//	for(listCone * cone = Poly->cones; cone; cone = cone->rest)
//	{
//		for(int i = 0; i < Poly->numOfVars; ++i) //find gcd of b in Ax <= b
//			dilationNum = (dilationNum *  (cone->vertex->vertex->denominators())[i])/ GCD(dilationNum, (cone->vertex->vertex->denominators())[i]);
//	}


}//ReadPolyhedronDataRecursive

int ReadPolyhedronDataRecursive::getNumberEqualities() const
{
	return set_card(ddHrep->linset);
}

int ReadPolyhedronDataRecursive::getNumberRows() const
{
	return ddHrep->rowsize;
}


void ReadPolyhedronDataRecursive::latticeInverse()
{

	//cout << "got to lattice Inverse()" << endl;

	if ( latticeBasis.NumCols() == latticeBasis.NumRows())
		return; //this is not needed, the polytope is full dimentional.

	assert(latticeBasis.NumCols() < latticeBasis.NumRows()); //we assume the columns are independent.

	RationalNTL **A; //adjoint matrix A=[L | I ]
	const int row = latticeBasis.NumRows();
	const int col = latticeBasis.NumCols();
	int pivot;
	A = new RationalNTL*[row];
	for(int i = 0; i  < row; ++i)
	{
		A[i] = new RationalNTL[col + row];
		for(int j = 0; j < col; ++j)
			A[i][j] = latticeBasis[i][j];
		A[i][col + i] = 1;
	}

	//do GE on the first col columns of A.
	//assumes A has full column rank.
	for(int curCol = 0; curCol < col; ++curCol)
	{
		//find smallst pivot on col curCol and curRow >= curCol.
		pivot = -1;
		for(int i = curCol; i < row && pivot == -1; ++i)
		{
			if( A[i][curCol].getNumerator() > 0)
				pivot = i;
		}
		assert(pivot > -1);


		if ( pivot != curCol)
		{
			RationalNTL *t;
			t = A[curCol];
			A[curCol] = A[pivot];
			A[pivot] = t;
		}//swap rows (pointers)

		assert(A[curCol][curCol].getNumerator() > 0);

		//divide by leading element.
		for(int r = curCol+1; r < col + row; ++r)
		{
//			cout << "r" << r << "cur =" << curCol << "row " << row << "col " << col << endl;
			A[curCol][r] = A[curCol][r]/A[curCol][curCol];
		}
		A[curCol][curCol] = 1;

//		cout << "A[]" << curCol << " = " << A[curCol][curCol] << endl;
		//Eliminate the rows.
		for(int i = 0; i < row; ++i)
		{
			if ( i == curCol)
				continue;
			else if( A[i][curCol].getNumerator() != 0)
			{
				RationalNTL mult = A[i][curCol];
//				cout << "mult " << mult << endl;
				for(int r = curCol; r < col + row; ++r)
				{
//					cout << "r=" << r << endl;
					A[i][r] = A[i][r] - mult*A[curCol][r];
				}
			}//add rows.
		}//for each row other than the pivor row.

	}//GE.



/*	cout << "print left iverse matrix" << endl;
	for(int i = 0; i < row; ++i)
	{
		for(int j = 0; j < row; ++j)
			cout << A[i][col+j] << ",";
		cout << endl;
	}
	cout << "end of printing matix" << endl;
*/

	//check we really have the left-inverse.
	for(int i = 0; i < latticeBasis.NumCols(); ++i)
	{
		for(int j = 0; j < row; ++j)
		{
			RationalNTL sum;
			for(int k = 0; k < row; ++k)
				sum += A[j][col+k]*latticeBasis[k][i];
			if (j == i)
				assert(sum == 1);
			else
				assert(sum == 0);
		}// j th row of A times the i th col of lattceBasis.
	}

	//find the lcm of every denominator.
	ZZ theLCM;
	theLCM = 1;
	for(int i = 0; i < /*row*/ col; ++i)
		for(int j = col; j < col + row; ++j)
			theLCM = (theLCM*A[i][j].getDenominator())/GCD(theLCM, A[i][j].getDenominator());

	//same the factor and the integer left "inverse"
	latticeLeftInverseDilation = theLCM;
	latticeLeftInverse.kill();
	latticeLeftInverse.SetDims(/*row*/ col, row);
	for(int i = 0; i < /*row*/ col; ++i)
		for(int j = 0; j < row; ++j)
		{
			RationalNTL prod = (A[i][col+j]*latticeLeftInverseDilation);
			latticeLeftInverse[i][j] = prod.getNumerator();
			assert(prod.getDenominator() == 1);
		}

	//free memory.
	for(int i = 0; i < row; ++i)
		delete [] A[i];
	delete A;
}//latticeInverse()

void ReadPolyhedronDataRecursive::readHrepMatrixFromFile(string filename, BarvinokParameters *params)
{
	dd_MatrixPtr M;



	if (cddstyle[0] == 'y') {
		cout << "readHrepMatrixFromFile:: we can only work with latte h-reps currently, sorry." << endl;
		exit(1);
	} else {
		/* Read an input file in LattE format. */
		if (Vrepresentation[0] == 'y') {
			cout << "readHrepMatrixFromFile:: we can only work with latte h-reps currently, sorry." << endl;
			exit(1);
		} else {
			/* Not VREP. */
			CheckEmpty(filename.c_str());
			M = ReadLatteStyleMatrix(filename.c_str(), /* vrep: */false,
			/* homogenize: */false,
			/* nonnegative: */nonneg[0] == 'y');
		}
	}

	/* Now we have in M the H-representation*/
	if ( M->representation != dd_Inequality)
	{
		cout << "readHrepMatrixFromFile:: M is not an h-rep, error" << endl;
		exit(1);
	}

	/*start of read PolyhedronFromHrepMatrix */
	int numOfVars = M->colsize - 1; /* Number of variables, not including RHS. */

	params->read_time.start();
	polyhedronRedundancyCheck(redundancycheck, M);

	if ( set_card(M->linset) > 0)
	{
		latticeBasis = findLatticeBasis(M, numOfVars);
	}
	//matrix = projectOutVariables(M, numOfVars, Poly);
	//dd_FreeMatrix(M);


	//matrix = matrixTmp;
	params->read_time.stop();
	params->Number_of_Variables = M->colsize -1;
	cerr << params->read_time;

	ddHrep = M;
	/* Now matrix contains the new inequalities. */


	//ok, ddHrep now contains the polytpoe, and latticeBasis contains the basis for the lattice affine space of the polytope if there are equations.


	dilatePolytope();
}

//we are in charge of freeing the memory of matrix.
//void ReadPolyhedronDataRecursive::setInequalityMatrix(listVector *newMatrix)
//{
//	matrix = newMatrix;
//}

/**
 * Given the current matrix, sets the ith inequality to an equality and recomputes the new polytope.
 * The outline is:
 * 1) convert the listVector to dd_matrix
 * 2) set a facet to equality
 * 3) run redundancy check
 * 4) run projection.
 */
void ReadPolyhedronDataRecursive::setInequalityToEquality(int i, listVector * &newMatrix, BarvinokParameters &newParm)
{
		cout << "ReadPolyhedronDataRecursive::setInequalityToEquality is not finished, sorry" << endl;
	exit(1);

}//setInequalityToEquality

void ReadPolyhedronDataRecursive::setInequality(int row)
{
	assert( set_member(row, ddHrep->linset) == false);
	set_addelem(ddHrep->linset, row);
}

void ReadPolyhedronDataRecursive::setMatrix(dd_MatrixPtr m)
{
	//assert(Poly == NULL); //make sure we never computed the lattice, polyhedron, or other objects.
	//if(Poly)
	//	printListCone(Poly->cones, ddHrep->colsize -1);
	//cout << "setMatrix::end of print poly" << endl;
	ddHrep = m;
	latticeBasis.kill();
}



RationalNTL ReadPolyhedronDataRecursive::volumeCorrection(const RationalNTL & a) const
{
	return a / power(dilationNum, getFullDimensionCount());
}

