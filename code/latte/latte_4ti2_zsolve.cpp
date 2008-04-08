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
using namespace _4ti2_zsolve_;

LinearSystem<int> *
facets_to_4ti2_zsolve_LinearSystem(listVector *facets,
				   listVector *equalities,
				   int numOfVars)
{
  int num_facets = lengthListVector(facets);
  int num_equalities = lengthListVector(equalities);
  VectorArray<int> matrix(/*height:*/ num_facets + num_equalities,
			  /*width:*/ numOfVars);
  listVector *f;
  int row;
  for (f = facets, row = 0; f!=NULL; f=f->rest, row++) {
    int col;
    for (col = 0; col<numOfVars; col++) {
      matrix[row][col] = convert_ZZ_to_int(f->first[col]);
    }
  }
  for (f = equalities; f!=NULL; f=f->rest, row++) {
    int col;
    for (col = 0; col<numOfVars; col++) {
      matrix[row][col] = convert_ZZ_to_int(f->first[col]);
    }
  }
  int *rhs = new int[num_facets + num_equalities];
  for (row = 0; row<num_facets + num_equalities; row++)
    rhs[row] = 0;
  LinearSystem<int> *ls
    = new LinearSystem<int> (/*matrix:*/ matrix, /*rhs:*/ rhs,
			     /*free:*/ true,
			     /*lower:*/ 1,
			     /*upper:*/ -1);
  delete[] rhs;
  for (row = 0; row<num_facets; row++)
    ls->get_relation(row) = Relation<int>::LesserEqual;
  for (; row<num_facets + num_equalities; row++)
    ls->get_relation(row) = Relation<int>::Equal;
  return ls;
}
