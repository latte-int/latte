/* ReductionTestWithCPLEX.cpp -- Check reducibility of vectors using CPLEX
	       
   Copyright 2007 Matthias Koeppe

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

#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cassert>

#include "ReductionTestWithCPLEX.h"

using namespace std;

static void
open_matrix_file(const string &filename, ifstream &file,
		 int &num, int &dim)
{
  file.open(filename.c_str());
  if (!file.good()) {
    cerr << "Failed to open " << filename << endl;
    exit(1);
  }
  file >> num >> dim;
  if (!file.good()) {
    cerr << "Parse error reading file " << filename << endl;
    exit(1);
  }
}

ReductionTestWithCPLEX::ReductionTestWithCPLEX(const ReductionTestFactory &data) :
  verbose(data.verbose)
{
  if (data.rays_filename.size() == 0) {
    cerr << "--reduction=cplex requires the use of --reduction-rays-file=FILE" << endl;
    exit(1);
  }
  string generators_filename = data.rays_filename;
  int dim_generators;
  ifstream generators_file;
  open_matrix_file(generators_filename, generators_file, num_generators, dim_generators);
  dim = dim_generators;
  
  int status;
  env = CPXopenCPLEX(&status);
  if (status != 0) {
    cerr << "Failed to obtain CPLEX environent." << endl;
    abort();
  }

  lp = CPXcreateprob(env, &status, "repr");
  if (status != 0) abort();
  
  status = CPXnewrows(env, lp, dim, /*rhs:*/ NULL, /*sense:*/ NULL,
		      /*rngval:*/ NULL, /*rownames:*/ NULL);
  if (status != 0) abort();

  status = CPXnewcols(env, lp, num_generators, /*obj:*/ NULL,
		      /*lb:*/ NULL, /*ub:*/ NULL,
		      /*ctype:*/ NULL, /*colname:*/NULL);
  if (status != 0) abort();
  
  int i;
  for (i = 0; i<num_generators; i++) {
    status = CPXchgcoef(env, lp, -1, i, 1.0);
    if (status != 0) abort();
    int j;
    for (j = 0; j<dim; j++) {
      double x;
      generators_file >> x;
      status = CPXchgcoef(env, lp, j, i, x);
      if (status != 0) abort();
      
    }
    if (!generators_file.good()) {
      cerr << "Parse error reading generator file" << endl;
      exit(1);
    }
  }

  multipliers = new double[num_generators];
  indices = new int[num_generators];
  ctype = new char[num_generators];
  {
    int i;
    for (i = 0; i<num_generators; i++) {
      indices[i] = i;
      ctype[i] = CPX_INTEGER;
    }
  }
}

ReductionTestWithCPLEX::~ReductionTestWithCPLEX()
{
  int status;
  status = CPXfreeprob(env, &lp);
  if (status != 0) abort();
  status = CPXcloseCPLEX(&env);
  if (status != 0) abort();
}

bool
ReductionTestWithCPLEX::IsReducible(const vector<int> &v)
{
  bool is_reducible;
  int status;
  int j;
  for (j = 0; j<dim; j++) {
    double x = v[j];
    status = CPXchgcoef(env, lp, j, -1, x);
    if (status != 0) abort();
  }
  status = CPXdualopt(env, lp);
  if (status != 0) abort();
  
  int lpstat = CPXgetstat(env, lp);
  if (lpstat != CPX_STAT_OPTIMAL) {
    cerr << "LP solution status code: " << lpstat << endl;
    cout << "LP written out as repr.lp" << endl;
    status = CPXwriteprob(env, lp, "repr.lp", "LP");
    if (status != 0) abort();
    exit(1);
  }
  else {
    if (verbose) {
      cout << "CPLEX says: Vector lies in cone." << endl
	   << "Multipliers: ";
      status = CPXgetx(env, lp, multipliers, 0, num_generators - 1);
      if (status != 0) abort();
	
      for (j = 0; j<num_generators; j++)
	cout << multipliers[j] << " ";
      cout << endl;
    }
  }

  status = CPXchgprobtype(env, lp, CPXPROB_MILP);
  if (status != 0) abort();
  status = CPXchgctype(env, lp, num_generators, indices, ctype);
  if (status != 0) abort();
  status = CPXmipopt(env, lp);
  if (status != 0) abort();
  int mipstat = CPXgetstat(env, lp);
  if (mipstat != CPXMIP_OPTIMAL) {
    cerr << "MIP solution status code: " << mipstat << endl;
    status = CPXwriteprob(env, lp, "repr-ip.lp", "LP");
    if (status != 0) {
      cerr << "Cannot write repr-ip.lp" << endl;
      abort();
    }
    cout << "MIP written out as repr-ip.lp" << endl;
    is_reducible = false;
  }
  else {
    if (verbose) {
      cout << "CPLEX says: Vector is generated over the nonnegative integers." << endl
	   << "Multipliers: ";
      status = CPXgetmipx(env, lp, multipliers, 0, num_generators - 1);
      if (status != 0) abort();
	
      for (j = 0; j<num_generators; j++)
	cout << multipliers[j] << " ";
      cout << endl;
    }
    is_reducible = true;
  }
  status = CPXchgprobtype(env, lp, CPXPROB_LP);
  if (status != 0) abort();
  return is_reducible;
}

int ReductionTestWithCPLEX::GetDimension()
{
  return dim;
}

