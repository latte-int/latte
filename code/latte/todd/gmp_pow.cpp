#include "todd/gmp_pow.h"

mpq_class 
pow(mpq_class x, unsigned int n)
{
  mpz_class result_num;
  mpz_class result_den;

  mpz_pow_ui(result_num.get_mpz_t(), x.get_num_mpz_t(), n);
  mpz_pow_ui(result_den.get_mpz_t(), x.get_den_mpz_t(), n);

  return mpq_class(result_num,
		   result_den);
}

mpz_class 
pow(mpz_class x, unsigned int n)
{
  mpz_class result;
  mpz_pow_ui(result.get_mpz_t(), x.get_mpz_t(), n);
  return result;
}
