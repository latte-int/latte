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

static void eliminationRow(mat_ZZ& A, mat_ZZ& R, long offset); 
static void eliminationCol(mat_ZZ& A, mat_ZZ& L, long offset); 
static void diagonalizationIn(mat_ZZ& A, mat_ZZ& L, mat_ZZ& R, long offset); 
static bool check(const mat_ZZ& A, long offset); 


mat_ZZ SmithNormalFormIlio(const mat_ZZ& A, mat_ZZ& L, mat_ZZ& R)
{
   mat_ZZ copy1;
   copy1 = A;

   diagonalizationIn(copy1, L, R, 0);

   int min = copy1.NumRows();
   int i, j;
   
   ZZ g;
   
   for (i = 0; i < min; ++ i) {
              
      for ( j = i + 1; j < min; ++ j) {
              
         if (IsOne(copy1[i][i]))  break;
                 
         else if (IsZero(copy1[j][j])) continue;
         
         else if (IsZero(copy1[i][i])) {
            std::swap (copy1[i][i], copy1[j][j]);
            //switching columns and rows of A effectively, since non diagonal entries are zero
            vec_ZZ temp = L[j];
            L[j] = L[i];
            L[i] = temp;
            for (int k = 0; k < R.NumRows(); k++)
            {
               temp[k] = R[k][j];
               R[k][j] = R[k][i];
               R[k][i] = temp[k];
            }
         }
         
         else {

            ZZ x, y, myNeg;
                                                            
            XGCD (g, y, x, copy1[j][j], copy1[i][i]);
            x /= g;
            y /= g;

            for (int k = 0; k < R.NumRows(); k++)
            {
               R[k][j] += x * R[k][i];
               R[k][i] -= (copy1[i][i] / g) * R[k][j];
            }

            L[i] += y * L[j];
            L[j] *= to_ZZ(-1);
            L[j] -= (copy1[j][j] / g) * L[i];
            L[j] *= to_ZZ(-1);
            
            copy1[j][j] = copy1[j][j] / g;
            copy1[j][j] = copy1[j][j] * copy1[i][i];
            copy1[i][i] = g;
         }
      }
   }
   
   return copy1;
}

static void diagonalizationIn(mat_ZZ& A, mat_ZZ& L, mat_ZZ& R, long stuff)
{
   for (int offset = 0; offset < A.NumRows(); offset++)
   {
      do {
         eliminationRow(A, R, offset);
         eliminationCol(A, L, offset);
         
         if (!IsZero(A[offset][offset]))
         {
            for (long i = offset + 1; i < A.NumRows(); i++)
            {
               if (!IsZero(A[i][offset]))
               {
                  L[i] -= (A[i][offset] / A[offset][offset]) * L[offset];
                  A[i] -= (A[i][offset] / A[offset][offset]) * A[offset];
               }
                         
               if (!IsZero(A[offset][i]))
               {
                  ZZ temp = A[offset][i] / A[offset][offset];
                  for (int j = 0; j < R.NumRows(); j++)
                  {
                     R[j][i] -= temp * R[j][offset];
                     A[j][i] -= temp * A[j][offset];
                  }
               }
            }
         }
      }
      while (!check(A, offset));
   }
}

static void eliminationCol(mat_ZZ& A, mat_ZZ& L, long offset)
{
   if((A.NumRows() >= offset)) { return; }

   if (!IsZero(A[offset][offset])) 
   {
      ZZ g, s, t, y1, y2, h;		
                
      XGCD(g, s, t, A[offset][offset], A[offset + 1][offset]);
      XGCD(h, y2, y1, s, t);

      y1 /= h*to_ZZ(-1); // gamma
      y2 /= h; // alpha
      if (IsZero(y1)) { y2 *= to_ZZ(-1); } //make sure one of them is negated ..

      vec_ZZ tmp1, tmp2, tmp1B, tmp2B;
      tmp1.SetLength(A.NumCols()); tmp2.SetLength(A.NumCols());
      tmp1B.SetLength(L.NumCols()); tmp2B.SetLength(L.NumCols());	
      tmp1 = y1 * A[offset] + y2 * A[offset + 1];
      tmp2 = s * A[offset] + t * A[offset + 1];
      A[offset] = tmp1;
      A[offset + 1] = tmp2;

      tmp1B = y1 * L[offset] + y2 * L[offset + 1];
      tmp2B = s * L[offset] + t * L[offset + 1];
      L[offset] = tmp1B;
      L[offset + 1] = tmp2B;


      if (!IsZero (A[offset][offset])) 
      {
         ZZ q = to_ZZ(-1) * A[offset][offset] / g;
         A[offset] += q * A[offset + 1];
         L[offset] += q * L[offset + 1]; 
      }
   
      vec_ZZ tmp_v, tempCol;
      tmp_v.SetLength(A.NumRows()); tempCol.SetLength(A.NumRows());			
           
      g = A[offset + 1][offset];
      for (int i = 0; i < offset; i++) { tmp_v[i] = 0; }
      tmp_v[offset] = to_ZZ(1);           
      tmp_v[offset + 1] = to_ZZ(1);
   
      long p1 = offset + 2;
      for(long col_p2 = offset + 2; col_p2 < A.NumRows(); ++ col_p2, ++ p1)
      {
         XGCD (g, s, tmp_v[p1], g, A[col_p2][offset]);
         //at the end, g = A[0][m]*tmp_v[m] + A[0][m-1]*tmp_v[m-1] + ... + A[0][2]*tmp_v[2] + A[0][1] + 0					
         if (!IsOne(s))
         {
            for (long p2 = offset + 1; p2 != p1; ++ p2) 
            {
                tmp_v[p2] *= (s / g);
            }
         }
      }
     
      if (IsZero(g)) { return; }
   
      for (long tmp_c = offset; tmp_c < A.NumCols(); ++ tmp_c) 
      {
         tempCol[tmp_c] = 0;
         for (int tmp_r = offset; tmp_r < A.NumRows(); ++ tmp_r)
         {
             tempCol[tmp_c] += A[tmp_r][tmp_c] * tmp_v[tmp_r];
         }
      }
      for (long tmp_c = offset; tmp_c < A.NumCols(); ++ tmp_c)
      {
         A[offset][tmp_c] = tempCol[tmp_c];
      }

      L[offset] = tmp_v[offset] * L[offset];		
      for (long tmp_rB = offset + 1; tmp_rB < L.NumRows(); ++ tmp_rB)
      {
         L[offset] += tmp_v[tmp_rB] * L[tmp_rB];
      }
   }
   // A pivot is found
   
   ZZ g, tmp;
   g = A[offset][offset];
   if(IsZero(g)) { return; }
   
   long tmp_rB = offset + 1;
   for (long tmp_r = offset + 1; tmp_r < A.NumCols(); ++ tmp_r, ++ tmp_rB)
   {      
      if (!IsZero(A[tmp_r][offset]))
      {
         tmp = A[tmp_r][offset] / g * to_ZZ(-1);
         A[tmp_r] += tmp * A[offset];
         L[tmp_rB] += tmp * L[offset];
      }
   }
}

static void eliminationRow(mat_ZZ& A, mat_ZZ& R, long offset)
{
   if (A.NumCols() >= offset) { return; }		

   if (!IsZero(A[offset][offset])) {
      ZZ y1, y2, g, h, s, t;
      
      XGCD(g, s, t, A[offset][offset], A[offset][offset + 1]);
      XGCD(h, y2, y1, s, t);

      y1 /= h*to_ZZ(-1); // gamma
      y2 /= h; // alpha
      if (IsZero(y1)) { y2 *= to_ZZ(-1); } //make sure one of them is negated ..
      
      vec_ZZ tmp1, tmp2, tmp1B, tmp2B;
      tmp1.SetLength(A.NumRows()); tmp2.SetLength(A.NumRows());
      tmp1B.SetLength(R.NumRows()); tmp2B.SetLength(R.NumRows());
      for (long i = offset; i < A.NumRows(); i++)
      {
          tmp1[i] = A[i][offset] * y1 + A[i][offset + 1] * y2;
          tmp2[i] = A[i][offset] * s + A[i][offset + 1] * t;
      }
      for (long i = 0; i < R.NumRows(); i++)
      {
          tmp1B[i] = R[i][offset] * y1 + R[i][offset + 1] * y2;
          tmp2B[i] = R[i][offset] * s + R[i][offset + 1] * t;
      }

      for (long i = offset; i < A.NumRows(); i++)
      {
          A[i][offset] = tmp1[i];
          A[i][offset + 1] = tmp2[i];
      }
      
      for (long i = 0; i < R.NumRows(); i++)
      {
          R[i][offset] = tmp1B[i];
          R[i][offset + 1] = tmp2B[i];
      }

      if (!IsZero (A[offset][offset])) {

          ZZ q = (A[offset][offset] / g) * to_ZZ(-1);

          for (long i = offset; i < A.NumRows(); i++)
          {
              A[i][offset] += q * A[i][offset + 1];
          }
          for (long i = 0; i < R.NumRows(); i++)
          {
              R[i][offset] += q * R[i][offset + 1];
          }
      }

      vec_ZZ tmp_v;
      tmp_v.SetLength(A.NumCols());
      
      g = A[offset][offset + 1];
      for (long i = 0; i < R.NumCols(); ++ i) { if (i < offset) {tmp_v[i] = 0;} }
      tmp_v[offset] = to_ZZ(1);
      tmp_v[offset + 1] = to_ZZ(1);
      long p1 = offset + 2;
      for (long row_p2 = offset + 2; row_p2 < A.NumCols(); ++ row_p2, ++ p1) {
              
         XGCD(g, s, tmp_v[p1], g, A[offset][row_p2]); //g = s * A[0][1] + tmp_v[p1]*A[0][row_p2]
         if (!IsOne(s))
         {        
            for (long p2 = offset + 1; p2 != p1; ++ p2) 
            {
               tmp_v[p2] *= (s / g);
            }
         }   
      }
      
      if (IsZero(g)) { return; }
      
      for (long tmp_r = offset; tmp_r < A.NumRows(); ++ tmp_r) 
      {
         InnerProduct(A[tmp_r][offset], A[tmp_r], tmp_v);
      }
      
      for (int i = 0; i < R.NumRows(); i++)
      {
         R[i][offset] = R[i][offset]*tmp_v[offset];
         for (int j = offset + 1; j < R.NumCols(); j++)
         {
            R[i][offset] += R[i][j] * tmp_v[j];
         }
      }
   }             
   // after finding the pivot

   ZZ g, tmp;
   
   g = A[offset][offset];
   if(IsZero(g)) { return; }
   long tmp_cB = offset + 1;
   for (long tmp_c = offset + 1; tmp_c < A.NumCols(); ++ tmp_c, ++ tmp_cB)
   {
      // test if needing to update 
      if (!IsZero (A[offset][tmp_c]))
      {
         tmp = A[offset][tmp_c] / g * to_ZZ(-1);
         for (long i = offset; i < A.NumRows(); i++)
         {
            A[i][tmp_c] += tmp * A[i][offset];
         }
         for (long i = 0; i < R.NumRows(); i++)
         {
            R[i][tmp_cB] += tmp * R[i][offset];
         }
      }
   }
}

static bool check(const mat_ZZ& A, long offset)
{
   ZZ tmp;
   if (IsZero(A[offset][offset])) return true;

   for (int i = offset + 1; i < A.NumRows(); i++)
   {
     if (A[i][offset] != 0 || A[offset][i] != 0) { return false; }
   }
   return true;
}

