#include <iostream>
#include <cstdlib>

#include "todd/todd.h"
#include "todd/todd-expansion.h"

using namespace std;

void usage()
{
  cerr << "usage: test-todd X_1 ... X_d" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  if (argc < 2) usage();
  int dim = argc - 1;
  vector<mpz_class> x(dim);
  int i;
  for (i = 0; i<dim; i++) {
    x[i] = mpz_class(atoi(argv[1 + i]));
  }
  int k;
  for (k = 0; k<=dim; k++)
    cout << todd(dim, k, x) << " ";
  cout << endl;
  vector<mpq_class> todds = evaluate_todd(x);
  vector<mpq_class>::iterator it;
  for (it = todds.begin(); it!=todds.end(); it++)
    cout << *it << " ";
  cout << endl;
}
