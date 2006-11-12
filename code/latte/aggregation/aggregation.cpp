/* aggregation.cpp -- Aggregate equations a la Kannan (1983)

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

#include <iostream>
#include <cstdlib>

#include "ReadLatteStyle.h"
#include "latte_gmp.h"
#include "todd/gmp_pow.h"

using namespace std;

static void usage()
{
  cerr << "usage: aggregation ROW-INDICES...< INPUTFILE > OUTPUTFILE" << endl;
  exit(1);
}

dd_MatrixPtr KannanAggregation(dd_MatrixPtr matrix,
			       dd_rowset aggrows)
{
  int num_rows = matrix->rowsize;
  int num_cols = matrix->colsize;
  int num_aggrows = set_card(aggrows);
  int new_num_rows = num_rows - num_aggrows + 1;
  dd_MatrixPtr new_matrix
    = dd_CreateMatrix(new_num_rows, num_cols);
  matrix->numbtype = dd_Rational;
  matrix->representation = dd_Inequality;
  int aggrow_index = new_num_rows - 1; // 0-based
  set_addelem(new_matrix->linset, aggrow_index + 1);
  mpz_class b_max, A_max;
  int i, k;
  for (i = 0, k = 0; i<num_rows; i++) {
    if (set_member(i + 1, aggrows)) {
      // Check that rows are equations with non-negative entries,
      // and compute the maximum entry.
      if (!set_member(i + 1, matrix->linset)) {
	cerr << "Row " << i + 1 << " must be an equation." << endl;
	exit(1);
      }
      mpq_class b_coeff(matrix->matrix[i][0]);
      if (b_coeff < 0) {
	cerr << "RHS of row " << i + 1 << " must be non-negative." << endl;
	exit(1);
      }
      if (b_coeff.get_den() != 1) {
	cerr << "RHS of the equation given in row " << i + 1 << " must be an integer." << endl;
	exit(1);
      }
      mpz_class b_coeff_z = b_coeff.get_num();
      if (b_coeff_z > b_max)
	b_max = b_coeff_z;
      int j;
      for (j = 1; j<num_cols; j++) {
	mpq_class A_coeff(matrix->matrix[i][j]);
	if (A_coeff > 0) {
	  cerr << "All coefficients of the equation given in row " << i + 1 << " must be non-negative (i.e., non-positive in the input format b | -A)." << endl;
	  exit(1);
	}
	if (A_coeff.get_den() != 1) {
	  cerr << "All coefficients of the equation given in row " << i + 1 << " must be integers." << endl;
	  exit(1);
	}
	mpz_class A_coeff_z = A_coeff.get_num();
	if (-A_coeff_z > A_max)
	  A_max = -A_coeff_z;
      }
    }
    else {
      // Copy non-aggregated rows.
      int j;
      for (j = 0; j<num_cols; j++)
	dd_set(new_matrix->matrix[k][j], matrix->matrix[i][j]);
      if (set_member(i + 1, matrix->linset))
	set_addelem(new_matrix->linset, k + 1);
      k++;
    }
  }
  mpz_class lambda = 3 * num_aggrows * (num_cols - 1) * A_max * b_max;
  cerr << "Kannan aggregation: Using lambda = " << lambda << endl;
  mpz_class lmp1 = pow(lambda, num_aggrows + 1);
  mpz_class li = lambda;
  for (i = 0; i<num_rows; i++) {
    if (set_member(i + 1, aggrows)) {
      mpz_class multiplier = lmp1 + li;
      int j;
      for (j = 0; j<num_cols; j++) {
	mpq_class A_coeff(matrix->matrix[i][j]);
	mpq_class coeff(new_matrix->matrix[aggrow_index][j]);
	coeff += multiplier * A_coeff;
	dd_set(new_matrix->matrix[aggrow_index][j],
	       coeff.get_mpq_t());
      }
      li *= lambda;
    }
  }
  return new_matrix;
}

int main(int argc, const char **argv)
{
  if(argc < 2) usage();

  dd_MatrixPtr matrix = ReadLatteStyleMatrix(cin, false, "standard input");

  dd_rowset aggset;
  set_initialize(&aggset, matrix->rowsize);
  int i;
  for (i = 1; i<argc; i++) {
    int index;
    int result = sscanf(argv[i], "%d", &index);
    if (result != 1) usage();
    if (index >= 1 && index <= matrix->rowsize)
      set_addelem(aggset, index);
    else {
      cerr << "Bad row index: " << index << endl;
      exit(1);
    }
  }

  dd_MatrixPtr agg_matrix = KannanAggregation(matrix, aggset);
  dd_FreeMatrix(matrix);
  set_free(aggset);
  
  WriteLatteStyleMatrix(cout, agg_matrix);
  dd_FreeMatrix(agg_matrix);
  
  return 0;
}
