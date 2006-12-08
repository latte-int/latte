/* latte_cddlib.cpp -- Interface to cddlib

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
rays_to_cddlib_matrix(listVector *rays, int numOfVars)
{
  dd_set_global_constants();
  int num_rays = lengthListVector(rays);
  int numOfVars_hom = numOfVars + 1;
  dd_MatrixPtr matrix = dd_CreateMatrix(num_rays, numOfVars_hom);
  matrix->numbtype = dd_Rational;
  matrix->representation = dd_Generator;
  int i, j;
  listVector *ray;
  mpq_class x;
  for (i = 0, ray = rays; i<num_rays; i++, ray = ray->rest) {
    for (j = 0; j<numOfVars; j++) {
      x = convert_ZZ_to_mpq(ray->first[j]);
      dd_set(matrix->matrix[i][j + 1], x.get_mpq_t());
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

