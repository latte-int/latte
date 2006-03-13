#ifndef NTL_TO_LIDIA_H 
#define NTL_TO_LIDIA_H 

#include "latte_ntl.h"
#include <LiDIA/bigint.h>
#include <LiDIA/bigint_matrix.h>
#include "../ramon.h"

using namespace LiDIA;

/* NTL to LiDIA conversions */
bigint_matrix
convert_mat_ZZ_to_bigint_matrix(const mat_ZZ &);
bigint_matrix
convert_listVector_to_bigint_matrix(listVector *);
bigint*
convert_vec_ZZ_to_bigint_array(const vec_ZZ &);
void
copy_vec_ZZ_to_bigint_array(bigint *, const vec_ZZ &);

/* LiDIA to NTL conversions */
mat_ZZ
convert_bigint_matrix_to_mat_ZZ(const bigint_matrix &);

/* latte to NTL conversions */
mat_ZZ
convert_listVector_to_mat_ZZ(listVector *);
 
#endif

