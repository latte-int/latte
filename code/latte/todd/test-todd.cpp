/* test-todd.cpp -- Evaluation of the Todd function

   Copyright 2006 Matthias Koeppe

   This file is part of LattE.
   
   LattE is free software; you can redistribute it and/or modify it
   under the terms of the version 2 of the GNU General Public License
   as published by the Free Software Foundation.

   LattE is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with LattE; if not, write to the Free Software Foundation,
   Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
*/

#include <iostream>
#include <cstdlib>

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
  vector<mpq_class> todds = evaluate_todd(x);
  vector<mpq_class>::iterator it;
  for (it = todds.begin(); it!=todds.end(); it++)
    cout << *it << " ";
  cout << endl;
}
