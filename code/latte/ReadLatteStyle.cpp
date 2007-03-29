/* ReadLatteStyle.cpp -- Read input data in the LattE style

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

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cctype>
#include "ReadLatteStyle.h"
#include "latte_gmp.h"
#include "latte_cddlib.h"
#include "vertices/cdd.h"

static void check_stream(const istream &f, const char *fileName, const char *proc)
{
  if (!f.good()) {
    cerr << "Read error on input file " << fileName << " in " << proc << "." << endl;
    exit(1);
  }
};

dd_MatrixPtr ReadLatteStyleMatrix(const char *fileName, bool vrep, bool homogenize)
{
  ifstream f(fileName);
  if (!f) {
    cerr << "Cannot open input file " << fileName << " in ReadLatteStyleMatrix." << endl;
    exit(1);
  }
  return ReadLatteStyleMatrix(f, vrep, homogenize, fileName);
}

dd_MatrixPtr ReadLatteStyleMatrix(istream &f, bool vrep, bool homogenize,
				  const char *fileName)
{
  dd_set_global_constants();
  int numOfVectors, numOfVars;
  f >> numOfVectors >> numOfVars;
  int num_homog = homogenize ? 1 : 0;
  int numOfVars_hom = numOfVars + num_homog;
  check_stream(f, fileName, "ReadLatteStyleMatrix");
  dd_MatrixPtr matrix = dd_CreateMatrix(numOfVectors, numOfVars_hom);
  matrix->numbtype = dd_Rational;
  if (vrep) matrix->representation = dd_Generator;
  else matrix->representation = dd_Inequality;
  int i, j;
  mpq_class x;
  for (i = 0; i<numOfVectors; i++) {
    for (j = 0; j<numOfVars; j++) {
      f >> x;
      check_stream(f, fileName, "ReadLatteStyleMatrix");
      dd_set(matrix->matrix[i][j + num_homog], x.get_mpq_t());
    }
  }
  // Skip whitespace
  while (isspace(f.peek())) {
    char c;
    f.get(c);
  }
  while (!f.eof()) {
    // Interpret keywords
    char buffer[20];
    f.get(buffer, 20, ' ');
    if (strcmp(buffer, "linearity") == 0) {
      int num_linearity;
      f >> num_linearity;
      check_stream(f, fileName, "ReadLatteStyleMatrix");
      int i;
      for (i = 0; i<num_linearity; i++) {
	int index; /* 1-based */
	f >> index;
	check_stream(f, fileName, "ReadLatteStyleMatrix");
	set_addelem(matrix->linset, index /* 1-based */);
      }
    }
    else if (strcmp(buffer, "nonnegative") == 0) {
      if (vrep) {
	cerr << "Keyword `" << buffer << "' not valid in vrep mode." << endl;
	exit(1);
      }
      int num_nonnegative;
      f >> num_nonnegative;
      check_stream(f, fileName, "ReadLatteStyleMatrix");
      // Create a new matrix that includes non-negativity
      // inequalities. 
      dd_MatrixPtr new_matrix
	= dd_CreateMatrix(numOfVectors + num_nonnegative, numOfVars_hom);
      new_matrix->numbtype = dd_Rational;
      new_matrix->representation = dd_Inequality;
      for (i = 0; i<numOfVectors; i++)
	for (j = 0; j<numOfVars_hom; j++)
	  dd_set(new_matrix->matrix[i][j], matrix->matrix[i][j]);
      int k;
      for (k = 0; k<num_nonnegative; k++, i++) {
	int index; /* 1-based */
	f >> index;
	check_stream(f, fileName, "ReadLatteStyleMatrix");
	for (j = 0; j<numOfVars; j++)
	  dd_set_si(new_matrix->matrix[i][j + num_homog], 0);
	dd_set_si(new_matrix->matrix[i][index + num_homog], 1);
      }
      set_copy(new_matrix->linset, matrix->linset);
      dd_FreeMatrix(matrix);
      matrix = new_matrix;
    }
    else {
      cerr << "Unknown keyword `" << buffer << "' in input file " << fileName
	   << " in ReadLatteStyleMatrix." << endl;
      exit(1);
    }
    // Skip whitespace
    while (!f.eof() && isspace(f.peek())) {
      char c;
      f.get(c);
    }
  }
  return matrix;	
}

void WriteLatteStyleMatrix(const char *fileName, dd_MatrixPtr matrix)
{
  ofstream f(fileName);
  WriteLatteStyleMatrix(f, matrix);
}

void WriteLatteStyleMatrix(ostream &f, dd_MatrixPtr matrix)
{
  f << matrix->rowsize << " " << matrix->colsize << endl;
  int i;
  for (i = 0; i<matrix->rowsize; i++) {
    int j;
    for (j = 0; j < matrix->colsize; j++)
      f << matrix->matrix[i][j] << " ";
    f << endl;
  }
  int num_linearity = set_card(matrix->linset);
  if (num_linearity > 0) {
    f << "linearity " << num_linearity << " ";
    for (i = 1; i<=matrix->rowsize; i++)
      if (set_member(i, matrix->linset))
	f << i << " ";
    f << endl;
  }
}

Polyhedron *ReadLatteStyleVrep(const char *filename, bool homogenize)
{
  Polyhedron *P = new Polyhedron;
  if (homogenize) {
    /* Homogenize. */
    dd_MatrixPtr matrix = ReadLatteStyleMatrix(filename,
					       /* vrep: */ true,
					       /* homogenize: */ false);
    //dd_WriteMatrix(stdout, matrix);
    dd_ErrorType error;
    dd_rowset redundant = dd_RedundantRows(matrix, &error);
    check_cddlib_error(error, "ReadLatteStyleVrep");
    /* The non-redundant rows are the rays of the homogenization. */
    int i;
    listCone *cone = createListCone();
    P->numOfVars = matrix->colsize;
    vec_ZZ ray;
    ray.SetLength(matrix->colsize);
    for (i = 1; i<=matrix->rowsize; i++) {
      if (!set_member(i, redundant)) {
	int j;
	/* CDD has homogenization in the 0-th,
	   LattE expects it in the last coordinate. */
	for (j = 0; j < matrix->colsize - 1; j++)
	  ray[j] = convert_mpq_to_ZZ(matrix->matrix[i - 1][j + 1]);
	ray[matrix->colsize-1] = convert_mpq_to_ZZ(matrix->matrix[i - 1][0]);
	cone->rays = appendVectorToListVector(ray, cone->rays);
	cone->vertex = new Vertex(createRationalVector(P->numOfVars));
      }
    }
    dd_FreeMatrix(matrix);
    P->cones = cone;
    P->dualized = false;
    P->homogenized = true;
  }
  else {
    /* Don't homogenize. */
    P->cones = computeVertexConesFromVrep(filename, P->numOfVars);
    P->dualized = false;
    P->homogenized = false;
  }
  return P;
}
