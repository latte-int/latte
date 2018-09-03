/*
 * testLibraryExamples.cpp
 *
 *  Created on: Jul 12, 2011
 *      Author: bedutra
 *  This files contains examples on how to use latte functions in a
 *  library setting.
 *
 *  We will always be working with the 3-polytope that is a
 *  square pyramid sitting on top of a cube for the examples.
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cassert>

#include "valuation.h"
#include "count.h"
#include "latte_cddlib.h"


//global option on how to build the polytopes.
int USE_FILES = 0;
bool PRINT_POLYHEDRA = false ;


/*****************************************************************************/
/********** start of function headers. Nothing interesting here **************/
/*****************************************************************************/


using namespace std;

typedef enum
{
	vrepToCones, vrepToVertices, hrepToCones, hrepToVertices
} ReadProblemType;

//write h or v files
void building_h_file();
void building_v_file();

//Reading the polyhedron
Polyhedron * building_polyhedron(const ReadProblemType readProblemType,
		BarvinokParameters * parms);
Polyhedron * building_polyhedron_with_files(const ReadProblemType readProblemType,
		BarvinokParameters * parms);
Polyhedron * building_polyhedron_without_files(const ReadProblemType readProblemType,
		BarvinokParameters * parms);


//finding volumes and integrals.
void finding_integral_polynomial_coneDecomp();
void finding_integral_linearForm_triangulation();
void finding_volume_coneDecomp();
void finding_volume_triangulation();

//counting
void finding_lattice_count();
void finding_ehrhart_polynomial();
void finding_ehrhart_taylor();

//helper functions to make the example polytope and integrand.
dd_MatrixPtr getHrepMatrix();
dd_MatrixPtr getVrepMatrix();
linFormSum getLinearForms();
monomialSum getPolynomial();

//A B C D E F G H I J K L M N O P Q R S T U V W X Y Z.


string FILE_BASE;

int main(int argc, char *argv[])
{
	/*	//feel free to comment this back in if you want.
	cout << "Hello, this file is designed to help developers use LattE as\n"
			<< " a library. I recommend commenting out the examples you do\n"
			<< " not care to see and un-commenting the print statements\n"
			<< " to the examples you want to see\n" << endl;

	cout << "Press enter. Lots of things will be printed to the screen\n"
			<< " depending on how many examples you have compiled." << endl;
	cin.get(); //stop, so user can read the message above.

	*/

	//global variable. Used to make temp. files.
	FILE_BASE = argv[0];


	/* //I should delete this block
	int file, num;
	for(cin >> file >> num; file != -1; cin >> file >> num)
	{
		BarvinokParameters *params = new BarvinokParameters;
		USE_FILES = file;
		switch(num)
		{
		case 1:
			cout << "doing h to cones" << endl;
			building_polyhedron(hrepToCones, params);
			break;
		case 2:
			cout << "doing h to v" << endl;
			building_polyhedron(hrepToVertices, params);
			break;
		case 3:
			cout << "doing v to cones" << endl;
			building_polyhedron(vrepToCones, params);
			break;
		case 4:
			cout << "doing v to v" << endl;
			building_polyhedron(vrepToVertices, params);
			break;
		}
		delete params;
	}
	*/



	//Reading in a polytope. The "easy" way.
	//Note that these "building" functions are leaking memory because I
	//	am not deleting the returned polytopes.
	BarvinokParameters *params = new BarvinokParameters;
	building_polyhedron(hrepToCones, params);
	building_polyhedron(hrepToVertices, params);
	building_polyhedron(vrepToCones, params);
	building_polyhedron(vrepToVertices, params);




	//volume examples
	finding_volume_coneDecomp();
	finding_volume_triangulation();

	//integral examples.
	finding_integral_polynomial_coneDecomp();
	finding_integral_linearForm_triangulation();

	//counting examples. Note that the rational functions are
	//  written to the .rat file.
	finding_lattice_count();
	finding_ehrhart_polynomial();
	finding_ehrhart_taylor();

	cout << "ALL TESTS PASSED" << endl;
	return 0;
}//main()


void building_h_file()
{
	string hFile = FILE_BASE + ".hrep";
	dd_MatrixPtr M;

	M = getHrepMatrix();
	fstream hOut;
	hOut.open(hFile.c_str(), fstream::out);
	hOut << M->rowsize << " " << M->colsize << endl;
	for (int i = 0; i < M->rowsize; ++i)
	{
		for (int j = 0; j < M->colsize; ++j)
		{
			hOut << M->matrix[i][j] << " ";
		}
		hOut << endl;
	}
	hOut.close();

	dd_FreeMatrix(M);

}///building_h_file

void building_v_file()
{
	string vFile = FILE_BASE + ".vrep";
	dd_MatrixPtr M;

	//make the v file
	M = getVrepMatrix();
	fstream vOut;
	vOut.open(vFile.c_str(), fstream::out);
	vOut << M->rowsize << " " << M->colsize << endl;
	for (int i = 0; i < M->rowsize; ++i)
	{
		for (int j = 0; j < M->colsize; ++j)
			vOut << M->matrix[i][j] << " ";
		vOut << endl;
	}
	vOut.close();

	dd_FreeMatrix(M);

}//building_v_file()

/**
 * The easy way to have latte construct a polyhedron from your v or h
 * representation is to write it to a file.
 *
 *
 */
Polyhedron * building_polyhedron(const ReadProblemType readProblemType,
		BarvinokParameters * params)
{

	if( USE_FILES == 1 )
		return building_polyhedron_with_files(readProblemType, params);
	else
		return building_polyhedron_without_files(readProblemType, params);

}

Polyhedron * building_polyhedron_with_files(const ReadProblemType readProblemType,
		BarvinokParameters * params)
{
	string hFile = FILE_BASE + ".hrep";
	string vFile = FILE_BASE + ".vrep";

	ReadPolyhedronData readPolyData;
	Polyhedron *Poly;

	if (readProblemType == hrepToCones || readProblemType == hrepToVertices)
	{
		building_h_file();
	} else
	{
		building_v_file();
	}


	//start from the h file, find the vertex representation.
	if (readProblemType == hrepToVertices)
	{

		strcpy(readPolyData.dualApproach, "yes"); //compute the vertices

		//note that at this point, you can set additional options.
		//readPolyData.parse_option("--cdd"); //if the input is in a cdd-style file.
		//readPolyData.parse_option("--compute-vertex-cones=cdd"); //find the tangent cones with cdd
		//readPolyData.parse_option("--compute-vertex-cones=4ti2");// or 4ti2
		//readPolyData.parse_option("--redundancy-check=none"); //check h-rep for redundancies or not
		//readPolyData.parse_option("--redundancy-check=cddlib");
		readPolyData.parse_option("--redundancy-check=full-cddlib");
		//readPolyData.parse_option("--vrep"); //if the input file is a latte style v-rep file.

		readPolyData.parse_option(hFile.c_str()); //set the file name last.


		//params->max_determinant = 1; //if you want unimodular cones.

#ifdef HAVE_FORTYTWO_LIB
		parse_standard_triangulation_option("--triangulation=4ti2", params);
		parse_standard_dualization_option("--dualization=4ti2", params);
#endif

		//ok that's all the set up we need. Read in the polyhedron.
		Poly = readPolyData.read_polyhedron(params);

		params->Number_of_Variables = Poly->numOfVars;
		params->File_Name = (char *) readPolyData.filename.c_str();

		//hey, the next print statement shows we don't have vertices. We have the dual h-rep.
		if ( PRINT_POLYHEDRA )
			Poly->printPolyhedron();
		// homogenized: true: cones represent the homogenization of the polyhedron.
		//				false:cones represent the supporting cones of all vertices of the polyhedron.
		// dualized: if the cone is dualized.
		// NOTE: In the Polyhedron class, vertices are saved with the 'leanding 1' at the end (v, 1)
		//		 For hyperplanes, they are saved as 'ax - b < 0'. This is opposite for how they are
		//			saved in the v- and h-rep files.

		//to fix this, dualize the cone again.
		dualizeCones(Poly->cones, Poly->numOfVars, params);
		Poly->dualized = false;
		if ( PRINT_POLYHEDRA )
			Poly->printPolyhedron();

	}
	//start from the h file, find the tangent cones
	else if (readProblemType == hrepToCones)
	{
		strcpy(readPolyData.dualApproach, "no"); //compute the cones

		readPolyData.parse_option("--redundancy-check=full-cddlib");
		readPolyData.parse_option(hFile.c_str()); //set the file name last.

#ifdef HAVE_FORTYTWO_LIB
		parse_standard_triangulation_option("--triangulation=4ti2", params);
		parse_standard_dualization_option("--dualization=4ti2", params);
#endif

		//ok that's all the set up we need. Read in the polyhedron.
		Poly = readPolyData.read_polyhedron(params);

		params->Number_of_Variables = Poly->numOfVars;
		params->File_Name = (char *) readPolyData.filename.c_str();

		//hey, the next print statement shows we don't have the rays of the tangent cones,
		//  we have the facets of the tangent cones.
		if ( PRINT_POLYHEDRA )
			Poly->printPolyhedron();
		//to fix this, dualize the cone again to compute the rays.
		//but then the cone is dualized. So dualize back (it just swaps the ray/facet lists).
		dualizeCones(Poly->cones, Poly->numOfVars, params);
		dualizeCones(Poly->cones, Poly->numOfVars, params);
		if ( PRINT_POLYHEDRA )
			Poly->printPolyhedron();
	}
	//start from the v-rep file, find the vertices.
	else if (readProblemType == vrepToVertices)
	{
		strcpy(readPolyData.dualApproach, "yes"); //compute the vertices

		readPolyData.parse_option(vFile.c_str()); //set the file name last.
		readPolyData.parse_option("--vrep"); //the input file is a latte style v-rep file.

#ifdef HAVE_FORTYTWO_LIB
		parse_standard_triangulation_option("--triangulation=4ti2", params);
		parse_standard_dualization_option("--dualization=4ti2", params);
#endif

		//ok that's all the set up we need. Read in the polyhedron.
		Poly = readPolyData.read_polyhedron(params);

		params->Number_of_Variables = Poly->numOfVars;
		params->File_Name = (char *) readPolyData.filename.c_str();

		//the polyhedron is in the form you think it is!
		if ( PRINT_POLYHEDRA )
			Poly->printPolyhedron();

	} else if (readProblemType == vrepToCones)
	{

		strcpy(readPolyData.dualApproach, "no"); //compute the cones

		readPolyData.parse_option(vFile.c_str()); //set the file name last.
		readPolyData.parse_option("--vrep"); //the input file is a latte style v-rep file.

#ifdef HAVE_FORTYTWO_LIB
		parse_standard_triangulation_option("--triangulation=4ti2", params);
		parse_standard_dualization_option("--dualization=4ti2", params);
#endif

		//ok that's all the set up we need. Read in the polyhedron.
		Poly = readPolyData.read_polyhedron(params);

		params->Number_of_Variables = Poly->numOfVars;
		params->File_Name = (char *) readPolyData.filename.c_str();

		//the polyhedron is in the form you think it is!
		if ( PRINT_POLYHEDRA )
			Poly->printPolyhedron();

	} else
	{
		cerr << "unknown type" << endl;
		exit(1);
	}

	freeListVector(readPolyData.templistVec);
	freeListVector(readPolyData.matrix);
	return Poly;
}//building_polyhedron()


Polyhedron * building_polyhedron_without_files(const ReadProblemType readProblemType,
		BarvinokParameters * params)
{
	string hFile = FILE_BASE + ".hrep";
	string vFile = FILE_BASE + ".vrep";

	ReadPolyhedronData readPolyData;
	Polyhedron *Poly;

	//start from the h file, find the vertex representation.
	if (readProblemType == hrepToVertices)
	{
		dd_MatrixPtr M;
		M = getHrepMatrix();
		//note that at this point, you can set additional options.
		//readPolyData.parse_option("--compute-vertex-cones=cdd"); //find the tangent cones with cdd
		//readPolyData.parse_option("--compute-vertex-cones=4ti2");// or 4ti2
		//readPolyData.parse_option("--redundancy-check=none"); //check h-rep for redundancies or not
		//readPolyData.parse_option("--redundancy-check=cddlib");
		readPolyData.parse_option("--redundancy-check=full-cddlib");

#ifdef HAVE_FORTYTWO_LIB
		parse_standard_triangulation_option("--triangulation=4ti2", params);
		parse_standard_dualization_option("--dualization=4ti2", params);
#endif

		//ok that's all the set up we need. Read in the polyhedron.
		Poly = readPolyData.read_polyhedron(M, params, ReadPolyhedronData::computeVertices);

		if ( PRINT_POLYHEDRA )
			Poly->printPolyhedron();
	}
	//start from the h file, find the tangent cones
	else if (readProblemType == hrepToCones)
	{
		dd_MatrixPtr M;
		M = getHrepMatrix();
		readPolyData.parse_option("--redundancy-check=full-cddlib");

#ifdef HAVE_FORTYTWO_LIB
		parse_standard_triangulation_option("--triangulation=4ti2", params);
		parse_standard_dualization_option("--dualization=4ti2", params);
#endif

		//ok that's all the set up we need. Read in the polyhedron.
		Poly = readPolyData.read_polyhedron(M, params, ReadPolyhedronData::computePrimalCones);

		if ( PRINT_POLYHEDRA )
			Poly->printPolyhedron();
	}
	//start from the v-rep file, find the vertices.
	else if (readProblemType == vrepToVertices)
	{
		dd_MatrixPtr M;
		M = getVrepMatrix();
		readPolyData.parse_option("--redundancy-check=full-cddlib");

#ifdef HAVE_FORTYTWO_LIB
		parse_standard_triangulation_option("--triangulation=4ti2", params);
		parse_standard_dualization_option("--dualization=4ti2", params);
#endif

		//ok that's all the set up we need. Read in the polyhedron.
		Poly = readPolyData.read_polyhedron(M, params, ReadPolyhedronData::computeVertices);

		if ( PRINT_POLYHEDRA )
			Poly->printPolyhedron();
	} else if (readProblemType == vrepToCones)
	{
		dd_MatrixPtr M;
		M = getVrepMatrix();
		readPolyData.parse_option("--redundancy-check=full-cddlib");

#ifdef HAVE_FORTYTWO_LIB
		parse_standard_triangulation_option("--triangulation=4ti2", params);
		parse_standard_dualization_option("--dualization=4ti2", params);
#endif

		//ok that's all the set up we need. Read in the polyhedron.
		Poly = readPolyData.read_polyhedron(M, params, ReadPolyhedronData::computePrimalCones);

		if ( PRINT_POLYHEDRA )
			Poly->printPolyhedron();
	} else
	{
		cerr << "unknown type" << endl;
		exit(1);
	}

	freeListVector(readPolyData.templistVec);
	freeListVector(readPolyData.matrix);


	return Poly;
}//building_polyhedron()





/**
 * Integrate a polynomial over a polytope using the cone decomposition method.
 *
 * Integrating via the triangulation method is very similar. I will not give
 * an explicit example using the triangulation method on polynomials.
 *
 */
void finding_integral_polynomial_coneDecomp()
{
	RationalNTL ans;
	BarvinokParameters *params = new BarvinokParameters;
	Polyhedron * poly = building_polyhedron(hrepToCones, params);
	PolytopeValuation polytopeValuation(poly, *params);
	monomialSum polynomial = getPolynomial();

	//the cone-decomposition method needs to start from tangent-cones
	ans = polytopeValuation.findIntegral(polynomial,
			PolytopeValuation::integratePolynomialAsLinearFormCone);
	//If a dilation is needed, the polyhedron will now be dilated.
	//If you wanted to use the triangulation method, you just need to preplace
	// LawrenceIntegration with TriangulationIntegration. Also, you could
	// ask for a v-rep by replacing hrepToCones with vrepToVertices because
	// the triangulation method can start from both polyhedron descriptions.

	cout << "Integral using the cone decomposition method is \n"
			<< "Rational: " << ans << '\n' << "Real    : " << ans.to_RR()
			<< endl;

	assert( ans == RationalNTL("734124064/2079"));

	delete poly;
	delete params;
}

/**
 * Computes the integral of powers of linear forms using the triangulation
 * method.
 */
void finding_integral_linearForm_triangulation()
{
	RationalNTL ans;
	BarvinokParameters * params = new BarvinokParameters;
	Polyhedron * poly = building_polyhedron(hrepToVertices, params);
	PolytopeValuation polytopeValuation(poly, *params);
	linFormSum lForms = getLinearForms();

	ans = polytopeValuation.findIntegral(lForms,
			PolytopeValuation::integrateLinearFormTriangulation);

	cout << "Integral using the triangulation method is \n" << "Rational: "
			<< ans << '\n' << "Real    : " << ans.to_RR() << endl;

	//note the string, number too big for an int.
	assert( ans == RationalNTL("2448974416/15"));

	delete poly;
	delete params;
}//finding_integral_polynomial_coneDecomp


/**
 * Computes the volume using the cone decomposition method.
 */
void finding_volume_coneDecomp()
{
	RationalNTL ans;
	BarvinokParameters *params = new BarvinokParameters;

	Polyhedron * poly = building_polyhedron(vrepToCones, params);

	PolytopeValuation polytopeValuation(poly, *params);


	//the cone-decomposition method needs to start from tangent-cones
	ans = polytopeValuation.findVolume(PolytopeValuation::volumeCone);
	//if a dilation is needed, the polyhedron will now be dilated.

	cout << "Volume using the cone decomposition method is \n" << "Rational: "
			<< ans << '\n' << "Real    : " << ans.to_RR() << endl;

	assert( ans == RationalNTL(208,3));

	delete poly;
	delete params;
}//finding_volume_coneDecomp


/**
 * Computes the volume using the triangulation method.
 */
void finding_volume_triangulation()
{
	RationalNTL ans;
	BarvinokParameters *params = new BarvinokParameters;
	Polyhedron * poly = building_polyhedron(vrepToVertices, params);
	PolytopeValuation polytopeValuation(poly, *params);

	//note, the triangulation method could start from the tangent-cones or a v-rep.
	ans = polytopeValuation.findVolume(PolytopeValuation::volumeTriangulation);
	//the original polytope is not lost if we started from a v-rep.

	cout << "Volume using the triangulation method is \n" << "Rational: "
			<< ans << '\n' << "Real    : " << ans.to_RR() << endl;

	assert( ans == RationalNTL(208,3));

	delete poly;
	delete params;
}//finding_volume_triangulation


/**
 * Computes the number of lattice points in the polytope.
 * Unlike in the other examples, we cannot start from the Polhedron, we have
 * to start from the command line!
 *
 * The improvement to using this method to find the count instead of simply
 * calling system with the LattE count command is that the output is directly
 * returned to you here. If you where to call a system command, you would have
 * to parse the output which could case errors if the LattE output every changes.
 */
void finding_lattice_count()
{
	string fileName = FILE_BASE + ".hrep";
	building_h_file();

	char * argv[] = { "./dummyName", "--count-lattice-points", "" };
	argv[2] = (char *) fileName.c_str();

	CountAnswerContainer ans;
	ans = mainCountDriver(3, argv);

	cout << "The number of lattice points is " << ans.numLaticePoints << endl;

	assert( ans.numLaticePoints == to_ZZ(126));
}//finding_lattice_count

/**
 * Finds the ehrhart polynomial
 */
void finding_ehrhart_polynomial()
{
	string fileName = FILE_BASE + ".hrep";
	building_h_file();

	char * argv[] = { "./dummyName", "--ehrhart-polynomial", "" };
	argv[2] = (char *) fileName.c_str();

	CountAnswerContainer ans;
	ans = mainCountDriver(3, argv);

	const char * ehrhartPoly[4] = { "1", "35/3", "44", "208/3" };

	for (int i = 0; i < 4; ++i)
		assert( ans.ehrhart_coefficients[i] == mpq_class(ehrhartPoly[i]));

}//finding_ehrhart_polynomial


/*
 * Finds the first 5 terms in the taylor expansion.
 */
void finding_ehrhart_taylor()
{
	string fileName = FILE_BASE + ".hrep";
	building_h_file();

	char * argv[] = { "./dummyName", "--ehrhart-taylor=5", "" };
	argv[2] = (char *) fileName.c_str();

	CountAnswerContainer ans;
	ans = mainCountDriver(3, argv);

	const char * taylorPoly[6] = { "1", "126", "755", "2304", "5189t", "9826" };

	for (int i = 0; i < 6; ++i)
		assert( ans.seriesExpansion[i] == to_ZZ(taylorPoly[i]));

}//finding_ehrhart_taylor

/**
 * Simple function to set up a dd_matrix of facets.
 */
dd_MatrixPtr getHrepMatrix()
{
	int h[][9] = { { 0, 1, 0, 0 }, { 0, 0, 1, 0 }, { 0, 0, 0, 1 }, { 4, -1, 0,
			0 }, { 4, 0, -1, 0 }, { 32, 4, 0, -8 }, { 48, 0, -4, -8 }, { 48,
			-4, 0, -8 }, { 32, 0, 4, -8 } }; //latte-style h rep (0 <= b - Ax)


	dd_MatrixPtr matrix =
			dd_CreateMatrix(9 /*num rows*/, 4 /*num of vars + 1*/);

	matrix->representation = dd_Inequality; //h-rep.

	for (int i = 0; i < 9; ++i)
		for (int j = 0; j < 4; ++j)
			dd_set(matrix->matrix[i][j], mpq_class(h[i][j]).get_mpq_t());

	return matrix;
}//getHrepMatrix

/**
 * Create a v-rep matrix for the polytope with vertices
 * 0 0 0
 * 4 0 0
 * 0 4 0
 * 4 4 0
 * 0 0 4
 * 4 0 4
 * 0 4 4
 * 4 4 4
 * (the above is a 3-cube)
 * and
 * (2,2,5)
 * (the above is a square pyramid on top of a cube).
 */
dd_MatrixPtr getVrepMatrix()
{
	int v[][4] = { { 1, 0, 0 }, { 1, 4, 0 }, { 1, 0, 4 }, { 1, 4, 4 } }; //I'm to lazy to fill up the matrix one element at a time, so save a part of the vertices here.

	dd_MatrixPtr matrix =
			dd_CreateMatrix(9 /*num vertices*/, 4 /*num of vars + 1*/);

	matrix->representation = dd_Generator; //v-rep.

	for (int z = 0; z <= 1; ++z)
	{
		for (int i = 0; i < 4; ++i)
		{
			for (int j = 0; j < 3; ++j)
				dd_set(matrix->matrix[4 * z + i][j],
						mpq_class(v[i][j]).get_mpq_t());
			dd_set(matrix->matrix[4 * z + i][3], mpq_class(z * 4).get_mpq_t());
		}
	}
	dd_set(matrix->matrix[8][0], mpq_class(1).get_mpq_t());
	dd_set(matrix->matrix[8][1], mpq_class(2).get_mpq_t());
	dd_set(matrix->matrix[8][2], mpq_class(2).get_mpq_t());
	dd_set(matrix->matrix[8][3], mpq_class(5).get_mpq_t());

	return matrix;
}//getVrepMatrix


/**
 * Returns the powers of linear forms
 * 1/2(1x +2y + 3z)^2 + (30x +20y + 9z)^3 + 1
 */
linFormSum getLinearForms()
{
	linFormSum lform;
	lform.varCount = 3;
	vec_ZZ l;
	l.SetLength(3);

	FormLoadConsumer<RationalNTL> loader;
	loader.setFormSum(lform);

	//the coefficient must be adjusted by the (degree!)

	l[0] = 1;
	l[1] = 2;
	l[2] = 3;
	loader.ConsumeLinForm(RationalNTL(2, 2), 2, l);// add 1/2(1x +2y + 3z)^2
	//2!=2

	l[0] = 30;
	l[1] = 20;
	l[2] = 9;
	loader.ConsumeLinForm(RationalNTL(6, 1), 3, l);// add (30x +20y + 9z)^3
	//3!=6.

	l[0] = 8;
	l[1] = 11;
	l[2] = 123456789;
	loader.ConsumeLinForm(RationalNTL(1, 1), 0, l);// add 1 = 1*( anything)^0


	return lform;
}//getLinearForms

/**
 * Returns the polynomial
 * 1/4 x^2y^3z^5 +  3/4 x^2y^2z^2 + 1/9
 */
monomialSum getPolynomial()
{
	monomialSum polynomial;
	polynomial.varCount = 3;
	int exp[3];

	MonomialLoadConsumer<RationalNTL> loader;
	loader.setMonomialSum(polynomial);

	exp[0] = 2;
	exp[1] = 3;
	exp[2] = 5;
	loader.ConsumeMonomial(RationalNTL(1, 4), exp); //add 1/4 x^2y^3z^5

	exp[0] = 2;
	exp[1] = 2;
	exp[2] = 2;
	loader.ConsumeMonomial(RationalNTL(3, 4), exp); //add 3/4 x^2y^2z^2

	exp[0] = 0;
	exp[1] = 0;
	exp[2] = 0;
	loader.ConsumeMonomial(RationalNTL(1, 9), exp); //add 1/9

	return polynomial;
}//getPolynomial
