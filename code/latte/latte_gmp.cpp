#include "latte_gmp.h"

mpz_class
convert_ZZ_to_mpz(const ZZ &zz)
{
  long size = NumBytes(zz);
  unsigned char *data = new unsigned char[size];
  int sig = sign(zz);
  BytesFromZZ(data, zz, size);

  mpz_class mpz;
  mpz_import(mpz.get_mpz_t(), size, -1, 1, 1, 0, data);
  if (sig == -1)
    mpz = -mpz;
  delete[] data;
  return mpz;
}

ZZ
convert_mpz_to_ZZ(const mpz_class &mpz)
{
  size_t size;
  int sig = sgn(mpz);
  void *data = mpz_export(NULL, &size,
			  -1, 1, 1, 0,
			  mpz.get_mpz_t());
  ZZ result = ZZFromBytes((unsigned char *)data, size);
  if (sig == -1)
    result = -result;
  free(data);
  return result;
}

