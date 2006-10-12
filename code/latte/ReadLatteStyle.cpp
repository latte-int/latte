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

static void check_stream(const ifstream &f, const char *fileName, const char *proc)
{
  if (f.bad()) {
    cerr << "Read error on input file " << fileName << " in " << proc << "." << endl;
    exit(1);
  }
};

dd_MatrixPtr ReadLatteStyleMatrix(const char *fileName, bool vrep)
{
  dd_set_global_constants();
  ifstream f(fileName);
  if (!f) {
    cerr << "Cannot open input file " << fileName << " in ReadLatteStyleMatrix." << endl;
    exit(1);
  }
  int numOfVectors, numOfVars_hom;
  f >> numOfVectors >> numOfVars_hom;
  dd_MatrixPtr matrix = dd_CreateMatrix(numOfVectors, numOfVars_hom);
  matrix->numbtype = dd_Rational;
  if (vrep) matrix->representation = dd_Generator;
  else matrix->representation = dd_Inequality;
  int i, j;
  mpq_class x;
  for (i = 0; i<numOfVectors; i++) {
    for (j = 0; j<numOfVars_hom; j++) {
      f >> x;
      check_stream(f, fileName, "ReadLatteStyleMatrix");
      dd_set(matrix->matrix[i][j], x.get_mpq_t());
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
    else {
      cerr << "Unknown keyword " << buffer << " in input file " << fileName
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

static void check_cddlib_error(dd_ErrorType error, const char *proc)
{
  if (error != dd_NoError) {
    cerr << "CDDLIB error in " << proc << ": " << endl;
    dd_WriteErrorMessages(stderr, error);
    exit(1);
  }    
}

static ZZ
convert_mpq_to_ZZ(mpq_t mpq)
{
  mpq_class elt(mpq);
  assert(elt.get_den() == 1);
  return convert_mpz_to_ZZ(elt.get_num());
}

Polyhedron *ReadLatteStyleVrep(const char *filename, bool homogenize)
{
  Polyhedron *P = new Polyhedron;
  if (homogenize) {
    /* Homogenize. */
    dd_MatrixPtr matrix = ReadLatteStyleMatrix(filename,
					       true /* vrep */);
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
    P->cones = cone;
    P->dualized = false;
    P->homogenized == true;
  }
  else {
    /* Don't homogenize. */
    P->cones = computeVertexConesFromVrep(filename, P->numOfVars);
    P->dualized = false;
    P->homogenized = false;
  }
  return P;
}
