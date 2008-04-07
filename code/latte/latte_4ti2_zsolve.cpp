/* latte_4ti2_zsolve.cpp -- Interface to 4ti2's zsolve component
	       
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

#include <iostream>
#include "latte_gmp.h"
#include "latte_4ti2_zsolve.h"

using namespace std;

LinearSystem
facets_to_4ti2_zsolve_LinearSystem(listVector *facets,
				   listVector *equalities,
				   int numOfVars)
{
  int num_facets = lengthListVector(facets);
  int num_equalities = lengthListVector(equalities);
  Matrix matrix = createMatrix(numOfVars, num_facets + num_equalities);
  listVector *f;
  int row;
  for (f = facets, row = 0; f!=NULL; f=f->rest, row++) {
    int col;
    for (col = 0; col<numOfVars; col++) {
      matrix->Data[row * matrix->Width + col] = convert_ZZ_to_int(f->first[col]);
    }
  }
  for (f = equalities; f!=NULL; f=f->rest, row++) {
    int col;
    for (col = 0; col<numOfVars; col++) {
      matrix->Data[row * matrix->Width + col] = convert_ZZ_to_int(f->first[col]);
    }
  }
  LinearSystem ls = createLinearSystem();
  setLinearSystemMatrix(ls, matrix);
  deleteMatrix(matrix);
  Vector rhs = createZeroVector(num_facets + num_equalities);
  setLinearSystemRHS(ls, rhs);
  deleteVector(rhs);
  setLinearSystemLimit(ls, -1, -MAXINT, MAXINT, true);
  setLinearSystemEquationType(ls, -1, EQUATION_LESSEREQUAL, 0);
  for (row = num_facets; row<num_facets + num_equalities; row++)
    setLinearSystemEquationType(ls, row, EQUATION_EQUAL, 0);
  return ls;
}
