/* SubspaceAvoidingDecomposition.cpp -- Barvinok decomposition that avoids a prescribed subspace
	       
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

#include <cassert>
#include <iostream>
#include "barvinok/SubspaceAvoidingDecomposition.h"
#include "barvinok/Cone.h"

using namespace std;

/* The maximum norm. */
static ZZ
max_norm(const vec_ZZ& x)
{
  ZZ norm;
  int m = x.length();
  int i;
  for(i = 1; i <= m; i++) {
    ZZ absx = abs(x(i));
    if (norm < absx) norm = absx;
  }
  return norm;
}

vec_ZZ
ComputeShortVectorAvoidingSubspace(const mat_ZZ & B, const mat_ZZ &Dual)
{
  int dim = B.NumRows();
  assert(dim == B.NumCols());
   
  mat_ZZ L = -transpose(Dual);
  mat_ZZ U; 
  U.SetDims(dim, dim);

  ZZ det2;
  LLL(det2, L, U, 1, 1);

  // Find shortest basis vector that avoids the subspace.

  int i;
  ZZ least_max_norm;
  int least_max_norm_index;
  bool anyone = false;
  for(i = 1; i <= dim; i++) {
    if (U(i, dim) != 0) {
      ZZ n = max_norm(L(i));
      if (!anyone || n < least_max_norm) {
	least_max_norm = n;
	least_max_norm_index = i;
	anyone = true;
      }
    }
  }
  assert(anyone); // Full-dimensional, so we must have one vector that
		  // avoids the subspace.

  vec_ZZ result_U = U(least_max_norm_index);
  vec_ZZ result_L = L(least_max_norm_index);

#if 0
  // Try to greedily reduce the norm of result_L by adding or subtracting basis
  // vectors that lie in the subspace.
  bool change;
  cerr << "L = " << result_L << endl;
  do {
    change = false;
    for(i = 1; i <= dim; i++) {
      if (U(i, dim) == 0) {
	vec_ZZ tL = result_L - L(i);
	ZZ tn = max_norm(tL);
	if (tn < least_max_norm) {
	  least_max_norm = tn;
	  result_L = tL;
	  cerr << "L = " << result_L << endl;
	  result_U -= U(i);
	  change = true;
	}
	else {
	  tL = result_L + L(i);
	  tn = max_norm(tL);
	  if (tn < least_max_norm) {
	    least_max_norm = tn;
	    result_L = tL;
	    cerr << "L = " << result_L << endl;
	    result_U += U(i);
	    change = true;
	  }
	}
      }
    }
  } while (change);
#endif
  return result_U;
}

  
