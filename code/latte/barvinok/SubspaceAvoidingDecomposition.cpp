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

vec_ZZ
ComputeShortVectorAvoidingSubspace(const mat_ZZ & B, const mat_ZZ &Dual)
{
  int dim = B.NumRows();
  assert(dim == B.NumCols());
  vec_ZZ result = ComputeOmega(B, Dual, dim, 0, 0);
  if (result[dim-1] == 0) {
    std::cout << "*";
    result[dim-1] = 1;
  }
  return result;
}

  
