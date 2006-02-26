// This is a -*- C++ -*- header file.

#ifndef LATTE_GMP_H
#define LATTE_GMP_H

#include <gmpxx.h>
#include <vector>

typedef std::vector<mpq_class> mpq_vector;
typedef std::vector<mpz_class> mpz_vector;

#include <NTL/ZZ.h>
using namespace NTL_NAMESPACE;

/* Even though NTL is based on GMP, there are no functions that allow
   to convert to GMP.  Here are our own functions.
*/

mpz_class
convert_ZZ_to_mpz(const ZZ &zz);

ZZ
convert_mpz_to_ZZ(const mpz_class &mpz);

#endif

