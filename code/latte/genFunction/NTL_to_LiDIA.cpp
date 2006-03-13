#include "NTL_to_LiDIA.h"
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include <LiDIA/bigint.h>
#include <LiDIA/bigint_matrix.h>

using namespace std;

/* NTL to LiDIA conversions */

/* converts an NTL mat_ZZ matrix to a bigint_matrix */
bigint_matrix
convert_mat_ZZ_to_bigint_matrix(const mat_ZZ & ntl_m) {
cout << "convert_mat_ZZ_to_bigint_matrix...\n";
   int rows = ntl_m.NumRows();
   int cols = ntl_m.NumCols();
cout << "convert_mat_ZZ_to_bigint_matrix: dimensions [" << rows << ", " << cols << "]\n";
   bigint *lidia_v;
   bigint_matrix lidia_m;
   lidia_m.set_no_of_rows(rows);
   lidia_m.set_no_of_columns(cols);

   for(int i = 0; i < rows; i++) {
      /* add vector to lattice basis */
      lidia_v = convert_vec_ZZ_to_bigint_array(ntl_m[i]);
      lidia_m.sto_row(lidia_v, cols, i);
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
convert_vec_ZZ_to_bigint_array(const vec_ZZ& ntl_v) {
   bigint element;
   /* allocate memory for the bigint array */
   bigint *lidia_v = new bigint[ntl_v.length()];
   for (int i = 0; i < ntl_v.length(); i++) {
      /* copy from ZZ to bigint */
      element.assign(to_long(ntl_v[i]));
      lidia_v[i] = element;
   }
   return lidia_v;
}

/* copies a NTL vec_ZZ data type to a bigint * array */
void
copy_vec_ZZ_to_bigint_array(bigint *lidia_v, const vec_ZZ& ntl_v) {
   bigint vec_comp;

   for (int i = 0; i < ntl_v.length(); i++) {
      /* copy from vector to bigint */
      vec_comp.assign(to_long(ntl_v[i]));
      lidia_v[i] = vec_comp;
   }
}

/* LiDIA to NTL conversions */

/* converts a LiDIA bigint_matrix to a mat_ZZ */ 
mat_ZZ
convert_bigint_matrix_to_mat_ZZ(const bigint_matrix & lidia_m) {
   int rows = lidia_m.get_no_of_rows();
   int cols = lidia_m.get_no_of_columns();
   mat_ZZ ntl_m;
   long elm;

   ntl_m.SetDims(rows, cols);
   for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++) {
         (lidia_m.member(i,j)).longify(elm);
         ntl_m[i][j] = elm;
      }
   }
   return (ntl_m);
}

/* latte to NTL conversions */

/* converts a latte listVector to a mat_ZZ */ 
mat_ZZ
convert_listVector_to_mat_ZZ(listVector *list) {
   listVector *tmp_list = list;
   int rows = tmp_list->first.length();
   int cols = lengthListVector(tmp_list);
   int cur_col = 0;
   mat_ZZ m;
                                                                                
   m.SetDims(cols, rows);
                                                                                
   /* add columns as rows and then take the transpose */
   while (tmp_list) {
      m[cur_col] = copyVector(tmp_list->first, rows);
      cur_col++;
      tmp_list = tmp_list->rest;
   }
   return (transpose(m));
}
