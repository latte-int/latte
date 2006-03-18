#ifndef MATRIX_OPS_H 
#define MATRIX_OPS_H 

#include "latte_ntl.h"
#include "cone.h"

mat_ZZ
SmithNormalForm(const mat_ZZ &, mat_ZZ &, mat_ZZ &);

mat_ZZ
SmithNormalForm(listVector *, mat_ZZ &, mat_ZZ &);

#endif

