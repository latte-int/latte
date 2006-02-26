#include "gmp-interface.h"
#include <iostream>
using namespace std;

int main()
{
  ZZ zz;
  cin >> zz;
  cout << convert_ZZ_to_mpz(zz);
  mpz_class mpz;
  cin >> mpz;
  cout << convert_mpz_to_ZZ(mpz);
  return 0;
}
