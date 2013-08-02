/* cdd.cpp -- Computation of all vertex cones via CDD

 Copyright 2002 Raymond Hemmecke, Ruriko Yoshida
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

#include "config.h"
#include "cone.h"
#include "print.h"
#include "ramon.h"
#include "rational.h"
#include <list>
#include <cassert>
#include "latte_system.h"
#include "latte_relocatable.h"
#include "latte_cddlib.h"
#include "gnulib/pathmax.h"
#include "../LattException.h"

using namespace std;
/* ----------------------------------------------------------------- */

listVector* CopyListVector(listVector* A, int numOfVars)
{

	vec_ZZ v;
	listVector *L, *endL;

	v = createVector(numOfVars);
	L = createListVector(v);
	endL = L;

	while (A)
	{
		v = A -> first;
		endL->rest = createListVector(v);
		endL = endL->rest;
		A = A->rest;
	}

	return (L->rest);
}

/* ----------------------------------------------------------------- */
listCone* CopyListCones(listCone* RudyCones, int numOfVars,
		rationalVector* Opt_vertex)
{
	int s;
	listCone *cones, *endCones, *newCone, *tmp;

	cones = createListCone();
	endCones = cones;
	tmp = RudyCones;
	while (tmp)
	{

		if ((tmp->vertex->vertex->numerators() == Opt_vertex->numerators())
				&& (tmp->vertex->vertex->denominators()
						== Opt_vertex->denominators()))
		{
			newCone = createListCone();
			s = tmp->coefficient;
			newCone->coefficient = s;
			newCone->rays = CopyListVector(tmp->rays, numOfVars);
			newCone->facets = CopyListVector(tmp->facets, numOfVars);
			newCone->vertex = tmp->vertex;
			endCones->rest = newCone;
			endCones = endCones->rest;
		}
		tmp = tmp ->rest;
	}
	//  printListCone(cones->rest,numOfVars);
	return (cones->rest);
}

/* ----------------------------------------------------------------- */
listCone* CopyListCones(listCone* RudyCones, int numOfVars)
{
	int s;
	listCone *cones, *endCones, *newCone, *tmp;
	listVector *latticepoint;
	cones = createListCone();
	endCones = cones;
	tmp = RudyCones;
	while (tmp)
	{

		//    if((tmp->vertex->numerators()==Opt_vertex->numerators()) &&
		//(tmp->vertex->denominators()==Opt_vertex->denominators())){
		newCone = createListCone();
		latticepoint = createListVector(tmp->latticePoints->first);
		newCone->latticePoints = latticepoint;
		s = tmp->coefficient;
		newCone->coefficient = s;
		newCone->rays = CopyListVector(tmp->rays, numOfVars);
		newCone->facets = CopyListVector(tmp->facets, numOfVars);
		newCone->vertex = tmp->vertex;
		endCones->rest = newCone;
		endCones = endCones->rest;

		tmp = tmp ->rest;
	}
	//  printListCone(cones->rest,numOfVars);
	return (cones->rest);
}

/* ----------------------------------------------------------------- */

void createCddIneFile(listVector* matrix, int numOfVars)
{
	int i;
	listVector* tmp;

	ofstream OUT;
	OUT.open("latte_cdd.ine");
	OUT << "H-representation" << endl;
	OUT << "begin " << endl;
	OUT << lengthListVector(matrix) << " " << numOfVars << " integer" << endl;
	tmp = matrix;
	while (tmp)
	{
		for (i = 0; i < (numOfVars); i++)
			OUT << (tmp->first)[i] << " ";
		OUT << endl;
		tmp = tmp->rest;
	}
	OUT << "end" << endl;
	OUT << "adjacency" << endl;
	OUT << "incidence" << endl;
	OUT.close();

	return;
}

void createCddIneFile(const dd_MatrixPtr M)
{
	ofstream OUT;
	OUT.open("latte_cdd.ine");
	OUT << "H-representation" << endl;
	OUT << "begin " << endl;
	OUT << M->rowsize << " " << M->colsize << ((M->numbtype == dd_Integer) ? " integer" : " rational") << endl;

	//print out matrix and add equation index to stack array.
	vector<int> stack;
	for(int i = 0; i < M->rowsize; ++i)
	{
		for(int j = 0; j < M->colsize; ++j)
			OUT << M->matrix[i][j] << " ";
		OUT << endl;

		if ( set_member(i+1, M->linset) )
			stack.push_back(i+1);
	}
	OUT << "end" << endl;
	OUT << "adjacency" << endl;
	OUT << "incidence" << endl;


	if ( stack.size())
	{
		OUT << "partial_enumeration " << stack.size();
		for (int i = 0; i < stack.size(); ++i)
		{
			OUT << " " << stack[i];
		}
		OUT << endl;
	}//if equations exist


	OUT.close();

	return;
}//createCddIneFile
/* ----------------------------------------------------------------- */
void createLrsIneFile(listVector* matrix, int numOfVars)
{
	int i;
	listVector* tmp;

	ofstream OUT;
	OUT.open("latte_lrs.ine");
	OUT << "H-representation" << endl;
	OUT << "begin" << endl;
	OUT << lengthListVector(matrix) << " " << numOfVars << " integer" << endl;
	tmp = matrix;
	while (tmp)
	{
		for (i = 0; i < (numOfVars); i++)
			OUT << (tmp->first)[i] << " ";
		OUT << endl;
		tmp = tmp->rest;
	}
	OUT << "end" << endl;
	OUT.close();

	return;
}

/* ----------------------------------------------------------------- */
void createLrsIneFileToPostAnalysys(listVector* matrix, int numOfVars)
{

	ifstream in;
	ofstream OUT;
	OUT.open("latte_cdd.ine");
	in.open("latte_lrs.ine");
	int x, y;
	string tmpString;
	while (tmpString != "begin")
		getline(in, tmpString);
	OUT << "H-representation" << endl;
	OUT << "begin" << endl;
	in >> x >> y;
	getline(in, tmpString);
	OUT << x << " " << y << " integer" << endl;
	while (tmpString != "end")
	{
		getline(in, tmpString);
		OUT << tmpString << endl;
	}

	OUT << "postanalysis " << endl;
	OUT << "adjacency" << endl;

	OUT.close();

	return;
}

/* ----------------------------------------------------------------- */
void createLrsExtFileToPostAnalysys(listVector* matrix, int numOfVars)
{
	ifstream in, in2;
	ofstream OUT;
	OUT.open("latte_cdd.ext");
	in.open("latte_lrs.ext");
	in2.open("latte_lrs.ext");
	int counter = 0;
	string tmpString;
	while (tmpString != "begin")
		getline(in2, tmpString);
	getline(in2, tmpString);
	while (tmpString != "end")
	{
		getline(in2, tmpString);
		counter++;
	}
	counter--;
	while (tmpString != "begin")
		getline(in, tmpString);
	OUT << "V-representation" << endl;
	OUT << "begin" << endl;
	getline(in, tmpString);
	OUT << counter << " " << numOfVars << " rational" << endl;
	while (tmpString != "end")
	{
		getline(in, tmpString);
		OUT << tmpString << endl;
	}

	OUT << "hull " << endl;

	OUT.close();

	return;
}

/* ----------------------------------------------------------------- */

void createCddExtFile(listVector* matrix, int numOfVars)
{
	int i;
	listVector* tmp;

	ofstream OUT;
	OUT.open("latte_cdd.ext1");
	OUT << "V-representation" << endl;
	OUT << "begin" << endl;
	OUT << lengthListVector(matrix) << " " << numOfVars << " integer" << endl;
	tmp = matrix;
	while (tmp)
	{
		for (i = 0; i < (numOfVars); i++)
			OUT << (tmp->first)[i] << " ";
		OUT << endl;
		tmp = tmp->rest;
	}
	OUT << "end" << endl;
	OUT << "hull" << endl;
	OUT.close();

	return;
}

/* ----------------------------------------------------------------- */
void createCddExtFile2(const char* filename)
{
	int i, numOfVec, numOfVars;
	string tmpString;
	// listVector* tmp;
	ifstream in(filename);
	if (!in.good())
	{
		cerr << "Unable to open input file `" << filename << "'" << endl;
		exit(1);
	}
	in >> numOfVec >> numOfVars;
	ofstream OUT;
	getline(in, tmpString);
	OUT.open("latte_cdd.ext");
	OUT << "V-representation" << endl;
	OUT << "begin" << endl;
	OUT << numOfVec << " " << numOfVars << " integer" << endl;
	for (i = 0; i < numOfVec; i++)
	{
		getline(in, tmpString);
		OUT << tmpString << endl;
	}
	OUT << "end" << endl;
	OUT << "hull" << endl;
	OUT.close();

	return;
}
/* ----------------------------------------------------------------- */
/**
 * Creats a cdd ext file from a dd_matrix of rays.
 */
void createCddExtFile2(const dd_MatrixPtr M)
{

	if (M->representation != dd_Generator)
	{
		cerr << "dd_Generator vertex type expected" << endl;
		THROW_LATTE(LattException::pe_UnexpectedRepresentation);
	}

	ofstream OUT;
	OUT.open("latte_cdd.ext");
	OUT << "V-representation" << endl;
	OUT << "begin" << endl;
	OUT << M->rowsize << " " << M->colsize << " rational" << endl;
	for(int i = 0; i < M->rowsize; i++)
	{
		for(int j = 0; j < M->colsize; ++j)
			OUT << M->matrix[i][j] << " ";
		OUT << endl;
	}
	OUT << "end" << endl;
	OUT << "hull" << endl;
	OUT.close();

	return;
}//createCddExtFile2



/* ----------------------------------------------------------------- */
void createCddIneLPFile(listVector* matrix, int numOfVars, vec_ZZ & cost)
{
	int i;
	listVector* tmp;
	// cerr << cost << " " << numOfVars << endl;
	ofstream OUT;

	OUT.open("LP.ine");
	OUT << "H-representation" << endl;
	OUT << "begin " << endl;
	OUT << lengthListVector(matrix) << " " << numOfVars << " integer" << endl;
	tmp = matrix;
	while (tmp)
	{
		for (i = 0; i < (numOfVars); i++)
			OUT << (tmp->first)[i] << " ";
		OUT << endl;
		tmp = tmp->rest;
	}
	OUT << "end" << endl;
	OUT << "maximize" << endl;
	OUT << 0 << " ";
	for (i = 0; i < numOfVars - 1; i++)
		OUT << cost[i] << " ";
	OUT << endl;
	OUT.close();

	return;
}

/* ----------------------------------------------------------------- */

rationalVector* ReadLpsFile(int numOfVars, bool verbose = true)
{
	ifstream in("LP.lps");
	string tmpString;
	ZZ x, y;
	if (verbose)
	{
		cerr << "Reading .lps file...";
		cerr.flush();
	}
	rationalVector* OptVector;
	OptVector = createRationalVector(numOfVars);
	if (!in)
	{
		cerr << "Cannot open input file in ReadLpsFile." << endl;
		exit(1);
	}

	while (tmpString != "begin")
		getline(in, tmpString);
	in >> tmpString;
	for (int i = 0; i < numOfVars; i++)
	{
		in >> tmpString >> tmpString;

		x = 0;
		y = 0;
		ReadCDD(in, x, y);
		OptVector->set_entry(i, x, y);

	}
	if (verbose)
	{
		cerr << "done." << endl;
	}
	return OptVector;
}

/* ----------------------------------------------------------------- */

listVector* createListOfInequalities(listVector* matrix, int numOfVars)
{
	int i, j;
	ZZ g;
	vec_ZZ v;
	listVector *tmp, *inequalities, *endInequalities;

	/* Copying equality constraints as Az<=b, -Az<=-b. */

	tmp = matrix;
	inequalities = createListVector(createVector(numOfVars));
	endInequalities = inequalities;

	while (tmp)
	{
		v = createVector(numOfVars);
		for (i = 0; i < (numOfVars); i++)
			v[i] = -(tmp->first)[i + 1];

		/* Normalize vector by removing gcd of entries. */
		g = v[0];
		for (i = 1; i < numOfVars; i++)
			g = GCD(g, v[i]);
		g = abs(g);
		if (g != 1)
		{
			for (i = 0; i < numOfVars; i++)
				v[i] = (v[i]) / g;
		}

		endInequalities->rest = createListVector(v);
		endInequalities = endInequalities->rest;

		/*    for (i=0; i<(numOfVars); i++) v[i]=(tmp->first)[i+1]);
		 endInequalities->rest=createListVector(v);
		 endInequalities=endInequalities->rest; */
		tmp = tmp->rest;
	}

	/* Writing non-negativity constraints. */
	for (i = 0; i < (numOfVars); i++)
	{
		v = createVector(numOfVars);
		for (j = 0; j < (numOfVars); j++)
		{
			if (i == j)
				v[i] = 1;
			else
				v[j] = 0;
		}
		endInequalities->rest = createListVector(v);
		endInequalities = endInequalities->rest;
	}

	return (inequalities->rest);
}
/* ----------------------------------------------------------------- */
listCone* readCddExtFile(int &numOfVars)
{
	int i, j, numOfVertices;
	ZZ x, y, leadingX, leadingY;
	char cddInFileName[PATH_MAX];
	rationalVector *v;
	listCone *cones, *endCones, *c;
	string tmpString;

	cerr << "Reading .ext file...";
	cerr.flush();

	strcpy(cddInFileName, "latte_cdd.ext");

	ifstream in(cddInFileName);
	if (!in)
	{
		cerr << "Cannot open input file in readCddExtFile." << endl;
		exit(1);
	}

	while (tmpString != "begin")
		getline(in, tmpString);

	in >> numOfVertices >> numOfVars >> tmpString;

	cones = createListCone();
	endCones = cones;
	if (numOfVertices == 0)
	{
		cerr << "readCddExtFile:: Empty Polytope." << endl;
		ofstream OUT("numOfLatticePoints");
		OUT << 0 << endl;
		exit(0);
	}

	if (numOfVertices == 1)
	{
		char read = 'a';
		int flag = 0;
		ofstream OUT("numOfLatticePoints");
		// cout << tmpString << endl;
		//getline(in , tmpString);
		in.get(read);
		//cout << tmpString << endl;

		while (((read == '\n' || read == '\r') || read == ' ') || read == '\t')
		{
			in.get(read);
			if (read == '0')
			{
				cerr << "\n\nreadCddExtFile:: Unbounded polytope!" << endl << endl;
				exit(0);
			}
		}

		while (read != '\n' && read != '\r')
		{
			if (read == '/')
				flag = 1;

			in.get(read);
		}

		if (flag)
		{
			cerr << "Integrally empty Polytope." << endl;
			OUT << 0 << endl;

		}

		else
		{
			cerr << "\n\n*****  Total number of lattice points: " << 1
					<< " ****" << endl << endl;
			OUT << 1 << endl;

		}

		exit(0);
	}

	for (i = 0; i < numOfVertices; i++)
	{
		v = createRationalVector(numOfVars - 1);
		for (j = 0; j < numOfVars; j++)
		{
			x = 0;
			y = 0;
			ReadCDD(in, x, y);

			if (j > 0)
			{
				//the input in in the form (leadingX/leadingY, x/y, ...), we need (1, x/y * leadingY/leadingX, ...)
				v->set_entry(j - 1, x*leadingY, y*leadingX);
			} else
			{
				if (x == 0)
				{
					cerr << "\n\nreadCddExtFile:: Given polyhedron is unbounded!!!(2)\n\n";
					ofstream Empty("numOfLatticePoints");
					Empty << 0 << endl;
					exit(0);
				}
				//Ok, we have just read in the leading "1"...which could be something else if the vertex is not integer.
				//Save the rational number.
				leadingX = x;
				leadingY = y;
			}
		}
		c = createListCone();
		c->vertex = new Vertex(v);
		endCones->rest = c;
		endCones = endCones->rest;
	}

	in.close();

	cerr << "done.\n";

	listCone *result = cones->rest;
	freeCone(cones);
	return result;
}
/* ----------------------------------------------------------------- */
listCone* readCddEadFile(listCone* cones, int numOfVars)
{
	int i, j, k, numOfVertices, numOfRays;
	char cddInFileName[PATH_MAX];
	vec_ZZ v;
	rationalVector **vertices;
	listVector *rays, *endRays;
	listCone *tmp;
	string tmpString;

	cerr << "Reading .ead file...";
	cerr.flush();

	strcpy(cddInFileName, "latte_cdd.ead");

	ifstream in(cddInFileName);
	if (!in)
	{
		cerr << "Cannot open input file in readCddEadFile." << endl;
		exit(1);
	}

	while (tmpString != "begin")
		getline(in, tmpString);

	in >> numOfVertices;
	getline(in, tmpString);

	vertices = createArrayRationalVector(numOfVertices);

	tmp = cones;
	for (i = 0; i < numOfVertices; i++)
	{
		vertices[i] = tmp->vertex->vertex;
		tmp = tmp->rest;
	}
	tmp = cones;
	for (i = 0; i < numOfVertices; i++)
	{
		in >> k;
		if (i != (k - 1))
		{
			cerr
					<< "Vertex numbering in file latte_cdd.ead is not increasing!\n";
			system_with_error_check("rm -f latte_cdd.*");
			exit(1);
		}

		in >> numOfRays;
		in >> tmpString;

		rays = createListVector(createVector(numOfVars));
		endRays = rays;

		if (numOfRays >= 0)
		{
			for (j = 0; j < numOfRays; j++)
			{
				in >> k;
				v = constructRay(vertices[i], vertices[k - 1], numOfVars - 1);
				endRays->rest = createListVector(v);
				endRays = endRays->rest;
			}
		} else
		{
			// numOfRays<0 means that the list of adjacencies is inverted.
			numOfRays = -numOfRays;
			int r = 1;
			for (j = 0; j < numOfVertices - numOfRays; ++j)
			{
				in >> k;
				while (r < k)
				{
					v = constructRay(vertices[i], vertices[r - 1], numOfVars
							- 1);
					endRays->rest = createListVector(v);
					endRays = endRays->rest;
					++r;
				}
				++r;
			}
			while (r <= numOfVertices)
			{
				v = constructRay(vertices[i], vertices[r - 1], numOfVars - 1);
				endRays->rest = createListVector(v);
				endRays = endRays->rest;
				++r;
			}
		}
		tmp->rays = rays->rest;
		delete rays; // only deletes the dummy head
		tmp = tmp->rest;
	}
	delete[] vertices;

	in.close();
	cerr << "done.\n";

	return (cones);
}

/* ----------------------------------------------------------------- */
void CreatExtEadFile()
{
	char cddInFileName[PATH_MAX];
	string tmpString;

	strcpy(cddInFileName, "latte_cdd.out");
	ifstream in(cddInFileName);
	if (!in)
	{
		cerr << "Cannot open input file in readCddEadFile." << endl;
		exit(1);
	}

	while (tmpString != "end")
		getline(in, tmpString);
	getline(in, tmpString);

	ofstream outExt("latte_cdd.ext");
	while (tmpString != "end")
	{
		getline(in, tmpString);
		outExt << tmpString << endl;
	}
	getline(in, tmpString);

	ofstream outEad("latte_cdd.ead");
	while (tmpString != "end")
	{
		getline(in, tmpString);
		outEad << tmpString << endl;
	}

}

/* ----------------------------------------------------------------- */
listCone* readCddEadFileFromVrep(listCone* cones, int numOfVars)
{
	int i, j, k, numOfVertices, numOfRays, counter = 0;
	char cddInFileName[PATH_MAX];
	vec_ZZ v;
	rationalVector **vertices;
	listVector *rays, *endRays;
	listCone *tmp;
	string tmpString;

	cerr << "Reading .ead file...";
	cerr.flush();

	strcpy(cddInFileName, "latte_cdd.ead");
	int tmp_int;
	ifstream in(cddInFileName);
	if (!in)
	{
		cerr << "Cannot open input file in readCddEadFile." << endl;
		exit(1);
	}

	while (tmpString != "begin")
	{
		getline(in, tmpString);
		counter++;
		if (counter > 10)
		{
			cerr << "Redundant vertices!" << endl;
			exit(1);
		}
	}

	in >> numOfVertices >> tmp_int;

	vertices = createArrayRationalVector(numOfVertices);

	tmp = cones;
	for (i = 0; i < numOfVertices; i++)
	{
		vertices[i] = tmp->vertex->vertex;
		tmp = tmp->rest;
	}
	tmp = cones;
	for (i = 0; i < numOfVertices; i++)
	{
		in >> k;
		if (i != (k - 1))
		{
			cerr
					<< "Vertex numbering in file latte_cdd.ead is not increasing!\n";
			system_with_error_check("rm -f latte_cdd.*");
			exit(1);
		}

		in >> numOfRays;
		in >> tmpString;

		rays = createListVector(createVector(numOfVars));
		endRays = rays;

		if (numOfRays >= 0)
		{
			for (j = 0; j < numOfRays; j++)
			{
				in >> k;
				v = constructRay(vertices[i], vertices[k - 1], numOfVars - 1);
				endRays->rest = createListVector(v);
				endRays = endRays->rest;
			}
		} else
		{
			// numOfRays<0 means that the list of adjacencies is inverted.
			numOfRays += numOfVertices;
			int r = 1;
			for (j = 0; j < numOfRays; ++j)
			{
				in >> k;
				while (r < k)
				{
					v = constructRay(vertices[i], vertices[r - 1], numOfVars
							- 1);
					endRays->rest = createListVector(v);
					endRays = endRays->rest;
					++r;
				}
				++r;
			}
			while (r <= numOfVertices)
			{
				v = constructRay(vertices[i], vertices[r - 1], numOfVars - 1);
				endRays->rest = createListVector(v);
				endRays = endRays->rest;
				++r;
			}
		}

		tmp->rays = rays->rest;
		tmp = tmp->rest;
	}

	in.close();
	cerr << "done.\n";

	return (cones);
}

/* ----------------------------------------------------------------- */
listCone* computeVertexCones(const char* fileName, listVector* matrix,
		int numOfVars)
{
	char cddOutFileName[PATH_MAX], command[PATH_MAX];
	listCone *cones;

	/* Compute vertices and edges with cdd. */

	createCddIneFile(matrix, numOfVars + 1);

	cerr << "Computing vertices and edges with cdd...";
	cerr.flush();
	system_with_error_check(relocated_pathname(CDD_PATH)
			+ " latte_cdd.ine > latte_cdd.out");
	cerr << "done." << endl;

	///   strcpy(command,"cp latte_cdd.ext ");
	///   strcat(command,fileName);
	///   strcat(command,".ext");
	///   system_with_error_check(command);

	///   strcpy(command,"cp latte_cdd.ead ");
	///   strcat(command,fileName);
	///   strcat(command,".ead");
	///   system_with_error_check(command);

	{
		int ext_numOfVars;
		cones = readCddExtFile(ext_numOfVars);
		assert(ext_numOfVars == numOfVars+1);
	}
	cones = readCddEadFile(cones, numOfVars + 1);
	system_with_error_check("rm -f latte_cdd.*");

	///   strcpy(cddOutFileName,fileName);
	///   strcat(cddOutFileName,".cdd");
	///   printListConeToFile(cddOutFileName,cones,numOfVars);

	return (cones);
}



/* ----------------------------------------------------------------- */
listCone* computeVertexCones(const char* fileName, const dd_MatrixPtr M)
{
	char cddOutFileName[PATH_MAX], command[PATH_MAX];
	listCone *cones;

	/* Compute vertices and edges with cdd. */

	createCddIneFile(M);

	cerr << "Computing vertices and edges with cdd...";
	cerr.flush();
	system_with_error_check(relocated_pathname(CDD_PATH)
			+ " latte_cdd.ine > latte_cdd.out");
	cerr << "done." << endl;

	{
		int ext_numOfVars;
		cones = readCddExtFile(ext_numOfVars);
		assert(ext_numOfVars == M->colsize);
	}
	cones = readCddEadFile(cones, M->colsize);
	system_with_error_check("rm -f latte_cdd.*");

	return (cones);
}//computeVertexCones
/* ----------------------------------------------------------------- */

listCone* computeVertexConesViaLrs(const char* fileName, listVector* matrix,
		int numOfVars)
{

	char cddOutFileName[PATH_MAX], command[PATH_MAX];
	listCone *cones;

	/* Compute vertices with lrs. */

	createLrsIneFile(matrix, numOfVars + 1);

	cerr << "Computing vertices with lrs...";
	system_with_error_check(LRS_PATH " latte_lrs.ine > latte_lrs.ext");
	cerr << "done.\n\n";

	createLrsIneFileToPostAnalysys(matrix, numOfVars + 1);
	createLrsExtFileToPostAnalysys(matrix, numOfVars + 1);

	cerr << "Computing edges with cdd...";
	system_with_error_check(relocated_pathname(CDD_PATH)
			+ " latte_cdd.ine > latte_cdd.out");
	cerr << "done.\n\n";

	///   strcpy(command,"cp latte_cdd.ext ");
	///   strcat(command,fileName);
	///   strcat(command,".ext");
	///   system_with_error_check(command);

	///   strcpy(command,"cp latte_cdd.ead ");
	///   strcat(command,fileName);
	///   strcat(command,".ead");
	///   system_with_error_check(command);

	{
		int ext_numOfVars;
		cones = readCddExtFile(ext_numOfVars);
		assert(ext_numOfVars == numOfVars+1);
	}
	cones = readCddEadFile(cones, numOfVars + 1);
	system_with_error_check("rm -f latte_cdd.* latte_lrs.*");

	///   strcpy(cddOutFileName,fileName);
	///   strcat(cddOutFileName,".cdd");
	///   printListConeToFile(cddOutFileName,cones,numOfVars);

	return (cones);
}

/* ----------------------------------------------------------------- */

/* from ComputeAdjacency.cpp */

dd_boolean SetInputFile(FILE **f, dd_DataFileType fname)
{
	dd_boolean success = dd_FALSE;
	success = dd_FALSE;

	if ((*f = fopen(fname, "r")) != NULL)
	{
		/*     printf("input file %s is open\n", fname); */
		success = dd_TRUE;
	} else
	{
		printf("The input file %s not found\n", fname);
	}
	return success;
}

dd_boolean SetWriteFile(FILE **f, dd_DataFileType fname)
{
	dd_boolean success = dd_FALSE;

	if ((*f = fopen(fname, "w")) != NULL)
	{
		/*     printf("output file %s is open\n",fname); */
		success = dd_TRUE;
	} else
	{
		printf("The output file %s cannot be opened\n", fname);
	}
	return success;
}

static int compute_adjacency(int argc, char *argv[])
{
	dd_MatrixPtr M = NULL, M2 = NULL, M3 = NULL;
	dd_SetFamilyPtr A = NULL;
	dd_colrange d;
	dd_ErrorType err = dd_NoError;
	dd_rowset redrows, linrows, ignoredrows, basisrows;
	dd_colset ignoredcols, basiscols;
	long rank;
	mytype val;
	FILE* out;
	int flag = 0;
	time_t starttime, endtime;
	dd_DataFileType inputfile;
	FILE *reading = NULL;

	//dd_set_global_constants();  /* First, this must be called. */
	out = fopen("latte_cdd.ead", "w");
	dd_init(val);
	if (argc > 1)
		strcpy(inputfile, argv[1]);
	if (argc <= 1 || !SetInputFile(&reading, argv[1]))
	{
		/*     dd_WriteProgramDescription(stdout); */
		/*     fprintf(stdout,"\ncddlib test program to remove redundancy and compute adjacency of the resulting representation.\n"); */
		dd_SetInputFile(&reading, inputfile, &err);
	}
	if (err == dd_NoError)
	{
		M = dd_PolyFile2Matrix(reading, &err);
	} else
	{
		fprintf(stderr, "Input file not found\n");
		goto _L99;
	}

	if (err != dd_NoError)
		goto _L99;

	if (M->representation == dd_Generator)
		d = M->colsize + 1;
	else
		d = M->colsize;

	/*   fprintf(stdout, "redundant rows:\n");*/
	time(&starttime);
	redrows = dd_RedundantRows(M, &err);
	time(&endtime);
	set_fwrite(out, redrows);
	/*   dd_WriteTimes(stdout,starttime,endtime); */

	M2 = dd_MatrixSubmatrix(M, redrows);
	if (M2->rowsize != M->rowsize)
	{
		fprintf(stderr, "redundant rows.\n");
		goto _L99;
	}

	/*   fprintf(stdout, "Implicit linearity (after removal of redundant rows): "); */
	linrows = dd_ImplicitLinearityRows(M2, &err);

	set_fwrite(stdout, linrows);
	set_uni(M2->linset, M2->linset, linrows);
	/* add the implicit linrows to the explicit linearity rows */

	/* To remove redundancy of the linearity part,
	 we need to compute the rank and a basis of the linearity part. */
	set_initialize(&ignoredrows, M2->rowsize);
	set_initialize(&ignoredcols, M2->colsize);
	set_compl(ignoredrows, M2->linset);
	rank = dd_MatrixRank(M2, ignoredrows, ignoredcols, &basisrows, &basiscols);
	set_diff(ignoredrows, M2->linset, basisrows);

	M3 = dd_MatrixSubmatrix(M2, ignoredrows);
	if (M3->rowsize != M2->rowsize)
	{
		fprintf(stderr, "redundant rows.\n");
		goto _L99;
	}

	A = dd_Matrix2Adjacency(M3, &err);

	dd_WriteSetFamily(out, A);

	dd_clear(val);
	set_free(linrows);
	set_free(basisrows);
	set_free(basiscols);
	set_free(ignoredrows);
	set_free(ignoredcols);
	set_free(redrows);

	if (A != NULL)
		dd_FreeSetFamily(A);
	dd_FreeMatrix(M);
	dd_FreeMatrix(M2);
	dd_FreeMatrix(M3);
	fclose(out);
	_L99: ;
	if (err != dd_NoError)
	{
		dd_WriteErrorMessages(stderr, err);
		return 1;
	} else
		return 0;
}


listCone *computeVertexConesFromExtFile(int &numOfVars)
{
	listCone *cones;

#if 0
	cerr << "Computing vertices and edges with cdd...";
	system_with_error_check(relocated_pathname(COMPUTEADJACENCY_PATH) + " latte_cdd.ext > latte_cdd.jnk 2>&1");
	cerr << "done.\n\n";
#else
	cerr << "Computing vertices and edges with cddlib...";
	// FIXME: This needs to be rewritten properly, avoiding the use of
	// files.
	{
		char *argv[2] =
		{ (char*) "", (char*) "latte_cdd.ext" };
		if (compute_adjacency(2, argv) != 0)
		{
			cerr << "failed." << endl;
			THROW_LATTE(LattException::bug_Unknown);
		};
	}
	cerr << "done.\n\n";
#endif
	//  CreatExtEadFile();

	///   strcpy(command,"cp latte_cdd.ext ");
	///   strcat(command,fileName);
	///   strcat(command,".ext");
	///   system_with_error_check(command);

	///   strcpy(command,"cp latte_cdd.ead ");
	///   strcat(command,fileName);
	///   strcat(command,".ead");
	///   system_with_error_check(command);

	{
		int ext_numOfVars;
		cones = readCddExtFile(ext_numOfVars);
		numOfVars = ext_numOfVars - 1;
	}
	//cout << "*****Start of printing list cone in computeVertexConesFromVrep\n";
	//printListCone(cones, numOfVars);
	//cout << "*****end of printing list cone in computeVertexConesFromVrep\n";

	cones = readCddEadFileFromVrep(cones, numOfVars + 1);

	system_with_error_check("rm -f latte_cdd.*");

	///   strcpy(cddOutFileName,fileName);
	///   strcat(cddOutFileName,".cdd");
	///   printListConeToFile(cddOutFileName,cones,numOfVars);

	return (cones);
}

/**
 * Writes the matrix v-rep to a file and then calls
 * computeVertexConesFromExtFile
 */
listCone* computeVertexConesFromVrep(const dd_MatrixPtr M, int &numOfVars)
{
	createCddExtFile2(M);

	return computeVertexConesFromExtFile(numOfVars);
}//computeVertexConesFromVrep

/**
 * Copes the latte file to a cdd v-rep file and calls
 * computeVertexConesFromExtFile
 */
listCone* computeVertexConesFromVrep(const char* fileName, int &numOfVars)
{
	/* Compute vertices and edges with cdd. */
	createCddExtFile2(fileName);
	return computeVertexConesFromExtFile(numOfVars);
}



/* ----------------------------------------------------------------- */
rationalVector* LP(listVector* matrix, vec_ZZ& cost, int numOfVars,
		bool verbose = true)
{
	rationalVector* Opt_vector;
	createCddIneLPFile(matrix, numOfVars + 1, cost);
	if (verbose)
	{
		cerr << "Computing LP... ";
		cerr.flush();
	}
	system_with_error_check(relocated_pathname(CDD_PATH) + " LP.ine > LP.out");
	if (verbose)
	{
		cerr << "done.";
		cerr.flush();
	}
	Opt_vector = ReadLpsFile(numOfVars, verbose);
	//  cerr << Opt_vector->numerators() << " " << Opt_vector -> denominator << endl;
	system_with_error_check("rm -f LP.*");

	return (Opt_vector);
}
/* ----------------------------------------------------------------- */
