/* SpecialSimplex.cpp -- Check for a special simplex using CPLEX
	       
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

#include "latte_gmp.h"

#include <cplex.h>

using namespace std;

static double
convert_ZZ_to_double(const ZZ &zz)
{
  mpz_class mpz = convert_ZZ_to_mpz(zz);
  return mpz.get_d();
}

void
FindSpecialSimplex(listCone *cone, int numOfVars)
{
  CPXENVptr env;
  int status;
  env = CPXopenCPLEX(&status);
  if (status != 0) {
    cerr << "Failed to obtain CPLEX environent." << endl;
    abort();
  }

  int num_rays = lengthListCone(cone->rays);

  lp = CPXcreateprob(env, &status, "repr");
  if (status != 0) abort();
  
  status = CPXchgprobtype(env, CPXPROB_MIP);
  if (status != 0) abort();


  status = CPXnewrows(env, lp, numOfVars + 2 * num_rays + 1, /*rhs:*/ NULL, /*sense:*/ NULL,
		      /*rngval:*/ NULL, /*rownames:*/ NULL);
  if (status != 0) abort();

  status = CPXnewcols(env, lp, 2 * numOfVars + 1, /*obj:*/ NULL,
		      /*lb:*/ NULL, /*ub:*/ NULL,
		      /*ctype:*/ NULL, /*colname:*/NULL);
  if (status != 0) abort();
  
  listVector *ray;
  int j;
  for (ray = cone->rays, j = 0; ray!=NULL; ray = ray->rest, j++) {
    int i;
    /* Fill equations that express that e_n is a linear combination of
       the rays. */
    for (i = 0; i<numOfVars; i++) {
      status = CPXchgcoef(env, lp, i, j, convert_ZZ_to_double(ray->first[i]));
      if (status != 0) abort();
    }
    /* Add constraints x_i <= y_i <= M x_i */
    int beg[2];
    int ind[4];
    double val[4];
    char sense[2];
    beg[0] = 0;  ind[0] = j; val[0] = +1; ind[1] = num_rays + j; val[1] = -M; sense[0] = 'L';
    beg[1] = 2;  ind[2] = j; val[2] = -1; ind[3] = num_rays + j; val[3] = +1; sense[1] = 'L';
    status = CPXaddrows(env, lp, 
    
  }
}
