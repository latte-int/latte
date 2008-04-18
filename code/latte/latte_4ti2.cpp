/* latte_4ti2.cpp -- Interface with 4ti2
	       
   Copyright 2007 Matthias Koeppe

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

#include "latte_4ti2.h"
#include "latte_gmp.h"

using namespace _4ti2_;

VectorArray *
rays_to_4ti2_VectorArray(listVector *rays, int numOfVars,
			 int num_homogenization_vars,
			 int num_extra_rows)
{
  int num_rays = lengthListVector(rays);
  int numOfVars_hom = numOfVars + num_homogenization_vars;
  VectorArray *va = new VectorArray(num_rays + num_extra_rows, numOfVars_hom);
  int i, j;
  listVector *ray;
  for (i = 0, ray = rays; i<num_rays; i++, ray = ray->rest) {
    for (j = 0; j<numOfVars; j++) {
      convert_ZZ_to_mpz(ray->first[j], (*va)[i][j + num_homogenization_vars]);
    }
  }
  return va;
}

VectorArray *
rays_to_transposed_4ti2_VectorArray(listVector *rays, int numOfVars,
				    int num_extra_rows)
{
  int num_rays = lengthListVector(rays);
  VectorArray *va = new VectorArray(numOfVars + num_extra_rows, num_rays);
  int i, j;
  listVector *ray;
  for (i = 0, ray = rays; i<num_rays; i++, ray = ray->rest) {
    for (j = 0; j<numOfVars; j++) {
      convert_ZZ_to_mpz(ray->first[j], (*va)[j][i]);
    }
  }
  return va;
}

