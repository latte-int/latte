/* latte_cddlib.cpp -- Interface to cddlib

   Copyright 2006, 2007 Matthias Koeppe

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

#include <cassert>

#include "latte_cddlib.h"
#include "latte_gmp.h"

init_cddlib_class::init_cddlib_class()
{
  dd_set_global_constants();
}

init_cddlib_class init_cddlib;

void check_cddlib_error(dd_ErrorType error, const char *proc)
{
  if (error != dd_NoError) {
    cerr << "CDDLIB error in " << proc << ": " << endl;
    dd_WriteErrorMessages(stderr, error);
    exit(1);
  }    
}

dd_MatrixPtr
rays_to_cddlib_matrix(listVector *rays, int numOfVars,
		      int num_homogenization_vars,
		      int num_extra_rows)
{
  dd_set_global_constants();
  int num_rays = lengthListVector(rays);
  int numOfVars_hom = numOfVars + num_homogenization_vars;
  dd_MatrixPtr matrix = dd_CreateMatrix(num_rays + num_extra_rows, numOfVars_hom);
  matrix->numbtype = dd_Rational;
  matrix->representation = dd_Generator;
  int i, j;
  listVector *ray;
  mpq_class x;
  for (i = 0, ray = rays; i<num_rays; i++, ray = ray->rest) {
    for (j = 0; j<numOfVars; j++) {
      x = convert_ZZ_to_mpq(ray->first[j]);
      dd_set(matrix->matrix[i][j + num_homogenization_vars], x.get_mpq_t());
    }
  }
  return matrix;
}

dd_PolyhedraPtr
cone_to_cddlib_polyhedron(listCone *cone, int numOfVars)
{
  dd_MatrixPtr matrix = rays_to_cddlib_matrix(cone->rays, numOfVars);
  dd_ErrorType error;
  dd_PolyhedraPtr poly = dd_DDMatrix2Poly(matrix, &error);
  check_cddlib_error(error, "cone_to_cddlib_polyhedron");
  return poly;
}

listCone *
cddlib_matrix_to_cone(dd_MatrixPtr matrix)
{
  int numOfVars = matrix->colsize - 1;
  assert(matrix->representation == dd_Generator);
  listCone *result = createListCone();
  result->vertex = new Vertex(new rationalVector(numOfVars));
  int i;
  for (i = matrix->rowsize - 1; i>=0; i--) {
    vec_ZZ ray;
    ray.SetLength(numOfVars);
    int j;
    {
      /* Check generator is homogeneous */
      mpq_class x(matrix->matrix[i][0]);
      assert(x == 0);
    }
    for (j = 0; j<numOfVars; j++) {
      ray[j] = convert_mpq_to_ZZ(matrix->matrix[i][j + 1]);
    }
    result->rays = new listVector(ray, result->rays);
  }
  return result;
}



mat_ZZ
cddlib_matrix_to_mat_ZZ(dd_MatrixPtr matrix)
{
  int col = matrix->colsize - 1;
  int row = matrix->rowsize - 1;
  mat_ZZ new_matrix;
  new_matrix.SetDims(row,col);

  for (int i = 1; i <= row; i++) {
    
    for (int j = 0; j<col; j++) {
      new_matrix(i,j) = convert_mpq_to_ZZ(matrix->matrix[i][j + 1]);
    }
  }
  return new_matrix;
}






void
cddlib_matrix_to_equations_and_inequalities(dd_MatrixPtr matrix,
					    listVector **equations,
					    listVector **inequalities)
{
  assert(matrix->representation == dd_Inequality);
  int numOfVars = matrix->colsize - 1;
  *equations = NULL;
  *inequalities = NULL;
  int i;
  for (i = matrix->rowsize - 1; i>=0; i--) {
    vec_ZZ ineq;
    ineq.SetLength(numOfVars + 1);
    int j;
    for (j = 0; j<=numOfVars; j++) {
      ineq[j] = convert_mpq_to_ZZ(matrix->matrix[i][j]);
    }
    if (set_member(i + 1, matrix->linset))
      (*equations) = new listVector(ineq, *equations);
    else
      (*inequalities) = new listVector(ineq, *inequalities);
  }
}

