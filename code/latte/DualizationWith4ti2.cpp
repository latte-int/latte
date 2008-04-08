/* DualizationWith4ti2.cpp -- Compute dual cones with 4ti2
	       
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

// From 4ti2:
#include "groebner/BitSet.h"
#include "groebner/VectorArrayStream.h"
#include "groebner/LatticeBasis.h"
#include "groebner/RayAlgorithm.h"

// From LattE:
#include "latte_gmp.h"
#include "latte_4ti2.h"
#include "dual.h"

using namespace _4ti2_;

static listVector *
listVectors_from_VectorArray(VectorArray *facets,
			     int numOfVars,
			     int index_offset)
{
  int num_facets = facets->get_number();
  vec_ZZ facet;
  facet.SetLength(numOfVars);
  listVector *result = NULL;
  int i;
  for (i = 0; i<num_facets; i++) {
    int j;
    for (j = 0; j<numOfVars; j++) {
      mpz_class z = (*facets)[i][j + index_offset];
      facet[j] = convert_mpz_to_ZZ(z);
    }
    result = new listVector(facet, result);
  }
  return result;
}

void
dualizeCone_with_4ti2(listCone *cone, int numOfVars)
{
  assert(cone->facets == NULL);
  assert(cone->subspace_generators == NULL);
  int num_rays = lengthListVector(cone->rays);
  /* Create a matrix from the rays, with extra coordinates 
     at the front for slack variables.
     (4ti2 does not use a homogenization coordinate.) */
  int lifted_dim = numOfVars + num_rays;
  BitSet *rs = new BitSet(lifted_dim);
  VectorArray *matrix
    = rays_to_4ti2_VectorArray(cone->rays, numOfVars,
			       /* num_homogenization_vars: */ num_rays,
			       /* num_extra_rows: */ 0);
  /* Add identity matrix for the slack variables. */
  {
    int i;
    for (i = 0; i<num_rays; i++) {
      (*matrix)[i][i] = 1;
      rs->set(i);
    }
  }
  VectorArray *facets = new VectorArray(0, matrix->get_size());
  lattice_basis(*matrix, *facets);
  VectorArray* subspace = new VectorArray(0, matrix->get_size());
  RayAlgorithm algorithm;
  algorithm.compute(*matrix, *facets, *subspace, *rs);

  cone->facets = listVectors_from_VectorArray(facets, numOfVars,
					      /*index_offset:*/ num_rays);
  cone->equalities = listVectors_from_VectorArray(subspace, numOfVars,
						  /*index_offset:*/ num_rays);
  
  delete facets;
  delete subspace;
  delete matrix;
  delete rs;
  /* Final swap data. */
  swap(cone->determinant, cone->dual_determinant);
  swap(cone->rays, cone->facets);
  swap(cone->subspace_generators, cone->equalities);
}
