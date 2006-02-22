#include <iostream>
#include <cstdlib>

#include "todd/todd.h"

using namespace std;

void usage()
{
  cerr << "usage: test-todd DIMENSION K X_1 ... X_DIMENSION" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  if (argc < 3) usage();
  int dim = atoi(argv[1]);
  int k = atoi(argv[2]);
  if (argc != dim + 3) usage();
  mpz_class *x = new mpz_class[dim];
  int i;
  for (i = 0; i<dim; i++) {
    x[i] = mpz_class(atoi(argv[3 + i]));
  }
  cout << todd(dim, k, x) << endl;
}
