/* matrix_ops.cpp -- Smith normal form

   Copyright 2006 Susan Margulies, Matthias Koeppe
   Copyright 2009 Stanislav Moreinis, Matthias Koeppe

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

#include "matrix_ops.h"
#include "config.h"


BarvinokParameters::SmithFormType
smith_form_type_from_name(const char *name)
{
  if (strcmp(name, "ilio") == 0) return BarvinokParameters::IlioSmithForm;
  else if (strcmp(name, "lidia") == 0) return BarvinokParameters::LidiaSmithForm;
  else {
    cerr << "Unknown Smith form type name: " << name << endl;
    exit(1);
  }
}

void
show_standard_smith_option(ostream &stream)
{
  stream << "  --smith-form={ilio,lidia}" << endl;
}

bool
parse_standard_smith_option(const char *arg,
			    BarvinokParameters *params)
{
  if (strncmp(arg, "--smith-form=", 13) == 0) {
    params->smithform = smith_form_type_from_name(arg + 13);
  }
  else return false;
  return true;
}

mat_ZZ
SmithNormalForm(const mat_ZZ& A, mat_ZZ& L, mat_ZZ& R,
		BarvinokParameters *params)
{
  switch (params->smithform) {
  case BarvinokParameters::LidiaSmithForm:
    return SmithNormalFormLidia(A, L, R);
  case BarvinokParameters::IlioSmithForm:
    return SmithNormalFormIlio(A, L, R);
  default:
    cerr << "Unknown Smith form type" << endl;
    exit(1);
  }
}

#ifdef HAVE_LIDIA

#include "NTL_to_LiDIA.h"
#include "lidia-include.h"

using namespace LiDIA;

/*
 * Converts rays to a LiDIA bigint_matrix, and then finds the SmithNormalForm
 * and accompanying matrices such that Smith(A) = BAC. Converts all matrices
 * back to NTL matrices. This will eventually be overwritten.
 */
mat_ZZ
SmithNormalFormLidia(const mat_ZZ & U, mat_ZZ & B, mat_ZZ & C) {
   bigint_matrix lidia_U;
   bigint_matrix lidia_snf_U;
   bigint_matrix lidia_B;
   bigint_matrix lidia_C;
   mat_ZZ snf_U;

   lidia_U = convert_mat_ZZ_to_bigint_matrix(U);
   lidia_snf_U = snf(lidia_U, lidia_B, lidia_C);
   snf_U = convert_bigint_matrix_to_mat_ZZ(lidia_snf_U);
   B = convert_bigint_matrix_to_mat_ZZ(lidia_B);
   C = convert_bigint_matrix_to_mat_ZZ(lidia_C);

   //print_debug_matrix(lidia_U);
   //print_debug_matrix(snf_U);

   return (snf_U);
}

/*
 * Converts rays to a LiDIA bigint_matrix, and then finds the SmithNormalForm
 * and accompanying matrices such that Smith(A) = BAC. Converts all matrices
 * back to NTL matrices. This will eventually be overwritten.
 */
// mat_ZZ
// SmithNormalForm(listVector *list, mat_ZZ & B, mat_ZZ & C) {
//    bigint_matrix lidia_U;
//    bigint_matrix lidia_snf_U;
//    bigint_matrix lidia_B;
//    bigint_matrix lidia_C;
//    mat_ZZ snf_U;

//    lidia_U = convert_listVector_to_bigint_matrix(list);
//    lidia_snf_U = snf(lidia_U, lidia_B, lidia_C);
//    snf_U = convert_bigint_matrix_to_mat_ZZ(lidia_snf_U);
//    B = convert_bigint_matrix_to_mat_ZZ(lidia_B);
//    C = convert_bigint_matrix_to_mat_ZZ(lidia_C);

//    //print_debug_matrix(lidia_U);
//    //print_debug_matrix(snf_U);

//    return (snf_U);
// }

#else
mat_ZZ SmithNormalForm(const mat_ZZ & U, mat_ZZ & B, mat_ZZ & C) {
  cout << "SmithNormalForm: This build is configured without LiDIA. " << endl;
  exit(1);
}
/// mat_ZZ
/// SmithNormalFormLidia(listVector *list, mat_ZZ & B, mat_ZZ & C) {
///   cout << "SmithNormalFormLidia: This build is configured without LiDIA. " << endl;
///   abort();
/// }
#endif


/* The remainder of this file contains an implementation
   of the Iliopoulos algorithm, derived from linbox.
   Author: Stanislav Moreinis
*/
using namespace NTL;

static bool check(const mat_ZZ& A, int offset)
{
    if (IsZero(A[offset][offset])) return true;
    
    ZZ tmp;
    for (int row_p = offset + 1; row_p < A.NumCols(); ++ row_p ) {
        rem(tmp, A[offset][row_p], A[offset][offset]);
        if (!IsZero(tmp)) { return false; }
    }
    return true;
}

static void eliminationRow(mat_ZZ& A, mat_ZZ& R, int offset)
{
  ZZ r, x1, x2, y1, y2, g, s;
  if (!IsZero(A[offset][offset]))
  {
    if (IsZero(A[offset][offset + 1]))
    {
      for (int i = 0; i < A.NumRows(); i++)
      {
	r = A[i][offset]; A[i][offset] = A[i][offset + 1]; A[i][offset + 1] = r;
	r = R[i][offset]; R[i][offset] = R[i][offset + 1]; R[i][offset + 1] = r;
      }
    }
    else
    {
      XGCD(r, x1, x2, A[offset][offset], A[offset][offset + 1]);
      
      y1 = A[offset][offset + 1] / r;
      y2 = to_ZZ(-1) * A[offset][offset] / r;

      for (int i = 0; i < A.NumRows(); i ++)
      {
	
	g = y1 * R[i][offset] + y2 * R[i][offset + 1];
	s = x1 * R[i][offset] + x2 * R[i][offset + 1];
	R[i][offset] = g;
	R[i][offset + 1] = s;
	
	g = y1 * A[i][offset] + y2 * A[i][offset + 1]; 
	s = x1 * A[i][offset] + x2 * A[i][offset + 1];
	A[i][offset] = g;
	A[i][offset + 1] = s;
      }
    }
  }
  vec_ZZ t_Coeffs; t_Coeffs.SetLength(A.NumCols() - offset);

  g = A[offset][A.NumCols() - 1]; //last element in row
  t_Coeffs[A.NumCols() - offset - 1] = to_ZZ(1);
  
  for (int i = 2; i < A.NumCols() - offset; i++) //from second last element in row to first element off the diagonal
  {
    XGCD(g, r, s, A[offset][A.NumCols() - i], g);
    t_Coeffs[A.NumCols() - offset - i] = r;
    if (!IsOne(s)) { for (int j = A.NumCols() - offset - i + 1; j < A.NumCols() - offset; j++) {t_Coeffs[j] *= s;} } //multiply everything after the current element's coefficients
  }
    
  for (int i = 0; i < A.NumRows(); i++)
  {
    for (int j = 1; j < A.NumCols() - offset; j++)
    {
      R[i][offset] += t_Coeffs[j] * R[i][offset + j];
      A[i][offset] += t_Coeffs[j] * A[i][offset + j];
    }
  }

  for (int j = offset + 1; j < A.NumCols(); j++)
  {
    if (!IsZero(A[offset][j]))
    {
      ZZ tmp = (A[offset][j] / g) * sign(A[offset][offset]);
      for (int i = 0; i < A.NumRows(); i++)
      {
	R[i][j] -= tmp * R[i][offset];
	A[i][j] -= tmp * A[i][offset];
      }
    }
  }
  t_Coeffs.kill();
}

static void eliminationCol(mat_ZZ& A, mat_ZZ& L, int offset)
{
  ZZ r, x1, x2, y1, y2;
  vec_ZZ tmp1, tmp2;
  if (!IsZero(A[offset][offset]))
  {
    if (IsZero(A[offset + 1][offset]))
    {
      tmp1 = L[offset]; L[offset] = L[offset + 1]; L[offset + 1] = tmp1;
      tmp1 = A[offset]; A[offset] = A[offset + 1]; A[offset + 1] = tmp1;
    }
    else
    {
      XGCD(r, x1, x2, A[offset][offset], A[offset + 1][offset]);
      y1 = A[offset + 1][offset] / r;
      y2 = to_ZZ(-1) * A[offset][offset] / r;
  
      tmp1 = y1 * L[offset] + y2 * L[offset + 1];
      tmp2 = x1 * L[offset] + x2 * L[offset + 1];
      L[offset] = tmp1;
      L[offset + 1] = tmp2;
      
      tmp1 = y1 * A[offset] + y2 * A[offset + 1];
      tmp2 = x1 * A[offset] + x2 * A[offset + 1];
      A[offset] = tmp1;
      A[offset + 1] = tmp2;
    }
  }
  
  ZZ g, s;
  vec_ZZ t_Coeffs; t_Coeffs.SetLength(A.NumRows() - offset);
  
  g = A[A.NumRows() - 1][offset]; //last element in row
  t_Coeffs[A.NumRows() - offset - 1] = to_ZZ(1);
  
  for (int i = 2; i < A.NumRows() - offset; i++) //from second last element in row to first element off the diagonal
  {
    XGCD(g, r, s, A[A.NumRows() - i][offset], g);
    t_Coeffs[A.NumRows() - offset - i] = r;
    if (!IsOne(s)) { for (int j = A.NumRows() - offset - i + 1; j < A.NumRows() - offset; j++) {t_Coeffs[j] *= s;} } //multiply everything after the current element's coefficients
  }
    
  for (int j = 1; j < A.NumCols() - offset; j++)
  {
    L[offset] += t_Coeffs[j] * L[offset + j];
    A[offset] += t_Coeffs[j] * A[offset + j];
  }
  
  for (int i = offset + 1; i < A.NumCols(); i++)
  {
    if (!IsZero(A[i][offset]))
    {
      s = (A[i][offset] / g) * sign(A[offset][offset]);
      L[i] -= s * L[offset];
      A[i] -= s * A[offset];
    }
  }
  t_Coeffs.kill();
}

static void diagonalizationIn(mat_ZZ& A, mat_ZZ& L, mat_ZZ& R)
{
  for (int offset = 0; offset < A.NumRows() - 1; offset++)
  {
    do
    {
        eliminationRow(A, R, offset);
        eliminationCol(A, L, offset);
    }
    while (!check(A, offset));
    
    if(A[offset][offset] < 0)
    {
      A[offset] *= to_ZZ(-1);
      L[offset] *= to_ZZ(-1);
    }
    
    if (!IsZero(A[offset][offset]))
    {
      for (int i = offset + 1; i < A.NumCols(); i++)
      {
	if (A[offset][i] != 0)
	{
	  ZZ tmp = (A[offset][i] / A[offset][offset]);
	  for (int j = 0; j < A.NumRows(); j++)
	  {
	    R[j][i] -= tmp * R[j][offset];
	    A[j][i] -= tmp * A[j][offset];
	  }
	}
      }
    }
  }
  if(A[A.NumRows() - 1][A.NumRows() - 1] < 0)
  {
    A[A.NumRows() - 1] *= to_ZZ(-1);
    L[A.NumRows() - 1] *= to_ZZ(-1);
  }
}

mat_ZZ SmithNormalFormIlio(const mat_ZZ& copy, mat_ZZ& L, mat_ZZ& R)
{
  mat_ZZ A; A = copy;
  ident(L, A.NumRows());
  ident(R, A.NumRows());

  diagonalizationIn(A, L, R);
  int min = A.NumRows() <= A.NumCols() ? A.NumRows() : A.NumCols();
  int i, j;
  vec_ZZ temp;  
  
  for (i = 0; i < min; ++ i) {
    
    for ( j = i + 1; j < min; ++ j) {
				    
      if (IsOne(A[i][i]))  break;
	      
      else if (IsZero(A[j][j])) continue;
      
      else if (IsZero(A[i][i])) {
	std::swap (A[i][i], A[j][j]);
	//switching columns and rows of A effectively, since non diagonal entries are zero
	temp = L[j];
	L[j] = L[i];
	L[i] = temp;
	for (int k = 0; k < R.NumRows(); k++)
	{
	  temp[k] = R[k][j];
	  R[k][j] = R[k][i];
	}
	for (int k = 0; k < R.NumRows(); k++)
	{
	  R[k][i] = temp[k];
	}
      }
      
      else {

	ZZ x, y, g, myNeg;
							
	XGCD (g, y, x, A[j][j], A[i][i]);
	
	for (int k = 0; k < A.NumRows(); k++)
	{
	  R[k][j] += x * R[k][i];
	  myNeg = R[k][i];
	  R[k][i] = R[k][j];
	  R[k][j] = myNeg;
	  R[k][j] -= (A[i][i] / g) * R[k][i];
	}
	
	L[i] += y * L[j];
	L[j] *= to_ZZ(-1);
	L[j] += (A[j][j] / g) * L[i];
	
	A[j][j] = A[j][j] / g;
	
	A[j][j] = A[j][j] * A[i][i];
  
	A[i][i] = g;
	
      }
    }
  }
  
  return A;
}

