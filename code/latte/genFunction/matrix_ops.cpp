#include "matrix_ops.h"
#include "config.h"

#ifdef HAVE_LIDIA

#include "NTL_to_LiDIA.h"
#include <LiDIA/bigint_matrix.h>
#include <LiDIA/bigint.h>

using namespace LiDIA;

/*
 * Converts rays to a LiDIA bigint_matrix, and then finds the SmithNormalForm
 * and accompanying matrices such that Smith(A) = BAC. Converts all matrices
 * back to NTL matrices. This will eventually be overwritten.
 */
mat_ZZ SmithNormalForm(const mat_ZZ & U, mat_ZZ & B, mat_ZZ & C) {
   bigint_matrix lidia_U;
   bigint_matrix lidia_snf_U;
   bigint_matrix lidia_B;
   bigint_matrix lidia_C;
   mat_ZZ snf_U;

cout << "SmithNormalForm(const mat_ZZ& ...)\n"; 
   lidia_U = convert_mat_ZZ_to_bigint_matrix(U);
cout << "SmithNormalForm:: convert_mat_ZZ_to_bigint_matrix\n";
   lidia_snf_U = snf(lidia_U, lidia_B, lidia_C);
cout << "SmithNormalForm:: snf\n"; 
   snf_U = convert_bigint_matrix_to_mat_ZZ(lidia_snf_U);
cout << "SmithNormalForm:: convert_bigint_matrix_to_mat_ZZ\n"; 
   B = convert_bigint_matrix_to_mat_ZZ(lidia_B);
   C = convert_bigint_matrix_to_mat_ZZ(lidia_C);

   return (snf_U);
}

/*
 * Converts rays to a LiDIA bigint_matrix, and then finds the SmithNormalForm
 * and accompanying matrices such that Smith(A) = BAC. Converts all matrices
 * back to NTL matrices. This will eventually be overwritten.
 */
mat_ZZ
SmithNormalForm(listVector *list, mat_ZZ & B, mat_ZZ & C) {
   bigint_matrix lidia_U;
   bigint_matrix lidia_snf_U;
   bigint_matrix lidia_B;
   bigint_matrix lidia_C;
   mat_ZZ snf_U;

   lidia_U = convert_listVector_to_bigint_matrix(list);
   lidia_snf_U = snf(lidia_U, lidia_B, lidia_C);
   snf_U = convert_bigint_matrix_to_mat_ZZ(lidia_snf_U);
   B = convert_bigint_matrix_to_mat_ZZ(lidia_B);
   C = convert_bigint_matrix_to_mat_ZZ(lidia_C);

   return (snf_U);
}

#else

mat_ZZ SmithNormalForm(const mat_ZZ & U, mat_ZZ & B, mat_ZZ & C) {
  cerr << "SmithNormalForm: This build is configured without LiDIA, " << endl
       << "which we require for computing the Smith normal form." << endl;
  abort();
}
mat_ZZ
SmithNormalForm(listVector *list, mat_ZZ & B, mat_ZZ & C) {
  cerr << "SmithNormalForm: This build is configured without LiDIA, " << endl
       << "which we require for computing the Smith normal form." << endl;
  abort();
}

#endif
