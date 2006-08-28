#include "NTL_to_LiDIA.h"
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include "lidia-include.h"
#include "latte_gmp.h"

using namespace std;

/* NTL to LiDIA conversions */

/* converts an NTL mat_ZZ matrix to a bigint_matrix */
bigint_matrix
convert_mat_ZZ_to_bigint_matrix(const mat_ZZ & ntl_m) {
   int rows = ntl_m.NumRows();
   int cols = ntl_m.NumCols();
   bigint *lidia_v;
   bigint_matrix lidia_m;
   lidia_m.set_no_of_rows(rows);
   lidia_m.set_no_of_columns(cols);

   for(int i = 0; i < rows; i++) {
      /* add vector to lattice basis */
      lidia_v = convert_vec_ZZ_to_bigint_array(ntl_m[i]);
      //print_debug_vector(lidia_v, ntl_m[i].length());
      /* LiDIA calls don't seem to work here...simply copy in by hand */        
      for (int j = 0; j < cols; j++) {
         lidia_m.sto(i,j, lidia_v[j]);
      }
      delete [] lidia_v;
   }
   return (lidia_m);
}


/* converts a latte listVector (an array of vec_ZZ) to a bigint_matrix */
bigint_matrix
convert_listVector_to_bigint_matrix(listVector *rays) {
   int rows = rays->first.length();
   int cols = lengthListVector(rays);
   int cur_col = 0;
   listVector *iter_list = rays;
   bigint *lidia_v;
   bigint_matrix lidia_m;
   lidia_m.set_no_of_rows(rows);
   lidia_m.set_no_of_columns(cols);

   while (iter_list) {
      /* add vector to lattice basis */
      lidia_v = convert_vec_ZZ_to_bigint_array(iter_list->first);
      lidia_m.sto_column(lidia_v, cols, cur_col);
      cur_col++;
      iter_list = iter_list->rest;
      delete [] lidia_v;
   }
   return (lidia_m);
}

/* converts a NTL vec_ZZ data type to a bigint * array */
bigint*
convert_vec_ZZ_to_bigint_array(const vec_ZZ& ntl_v)
{
  bigint element;
  /* allocate memory for the bigint array */
  bigint *lidia_v = new bigint[ntl_v.length()];
  copy_vec_ZZ_to_bigint_array(lidia_v, ntl_v);
  return lidia_v;
}

/* copies a NTL vec_ZZ data type to a bigint * array */
void
copy_vec_ZZ_to_bigint_array(bigint *lidia_v, const vec_ZZ& ntl_v)
{
  for (int i = 0; i < ntl_v.length(); i++) {
    /* copy from vector to bigint */
    mpz_class mpz = convert_ZZ_to_mpz(ntl_v[i]);
    lidia_v[i].assign(*mpz.get_mpz_t());
  }
}

/* LiDIA to NTL conversions */

/* converts a LiDIA bigint_matrix to a mat_ZZ */ 
mat_ZZ
convert_bigint_matrix_to_mat_ZZ(const bigint_matrix & lidia_m) {
   int rows = lidia_m.get_no_of_rows();
   int cols = lidia_m.get_no_of_columns();
   mat_ZZ ntl_m;

   ntl_m.SetDims(rows, cols);
   for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++) {
	//(lidia_m.member(i,j)).longify(elm);
	mpz_class mpz(&lidia_m.member(i,j).bigint_rep());
	ntl_m[i][j] = convert_mpz_to_ZZ(mpz);
      }
   }
   return (ntl_m);
}

void
print_debug_vector(const bigint * v, int len) {
   cout << "Begin vector: ["; 
   for (int i = 0; i < len; i++) {
      cout << v[i] << ","; 
   }
   cout << "]: End vector\n"; 
}

void
print_debug_matrix(const bigint_matrix & m) {
   int rows = m.get_no_of_rows(); 
   int cols = m.get_no_of_columns(); 

   cout << "Begin matrix:\n"; 
   for (int i = 0; i < rows; i++) {
   cout << "["; 
      for (int j = 0; j < cols; j++) {
         cout << m.member(i,j) << ","; 
      }
      cout << "]\n"; 
   }
   cout << ":End matrix\n"; 
}

