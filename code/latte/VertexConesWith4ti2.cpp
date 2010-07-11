/* VertexConesWith4ti2.cpp -- Compute vertex cones with 4ti2
	       
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

#include "latte_gmp.h"
#include "latte_4ti2.h"
#include "VertexConesWith4ti2.h"
#include "dual.h"
#include "DualizationWith4ti2.h"

#include "print.h"

using namespace _4ti2_;

void
computeVertexConesWith4ti2(listVector* ineqs, int numOfVars,
			   bool &unbounded,
			   ConeConsumer &consumer)
{
  unbounded = false;
  listCone *cones = NULL;
  int num_ineqs = lengthListVector(ineqs);
  /* Create a matrix from the facets, with extra coordinates
     at the front for slack variables.  The 1 is related to the fact
     that the vectors of MATRIX already have 1 + numOfVars entries.
  */
  int lifted_dim = 1 + numOfVars + num_ineqs;
  BitSet *rs = new BitSet(lifted_dim);
  VectorArray *matrix
    = rays_to_4ti2_VectorArray(ineqs, 1 + numOfVars,
			       /* num_homogenization_vars: */ num_ineqs,
			       /* num_extra_rows: */ 0);
  /* Add negative identity matrix for the slack variables. */
  {
    int i;
    for (i = 0; i<num_ineqs; i++) {
      (*matrix)[i][i] = -1;
      rs->set(i);
    }
  }
  /* Make the homogenization coordinate (corresponding to the right-hand side) non-negative */
  rs->set(num_ineqs);

#if 0
  {
    std::ofstream file("cone_for_4ti2_vertex_cones_computation");
    file << matrix->get_number() << " " << lifted_dim << endl;
    print(file, *matrix, 0, lifted_dim);
    cerr << "Created file `cone_for_4ti2_vertex_cones_computation'" << endl;
  }
#endif
  
  VectorArray *rays = new VectorArray(0, matrix->get_size());
  lattice_basis(*matrix, *rays);
  VectorArray* subspace = new VectorArray(0, matrix->get_size());
  RayAlgorithm algorithm;
  algorithm.compute(*matrix, *rays, *subspace, *rs);
  delete rs;

#if 0
  {
    std::ofstream file("4ti2_vertex_cones_computation_output");
    file << rays->get_number() << " " << lifted_dim << "\n";
    print(file, *rays, 0, lifted_dim);
    cerr << "Created file `4ti2_vertex_cones_computation_output'" << endl;
  }
#endif
  
  assert(subspace->get_number() == 0); /* We assume polytopes,
					  thus a pointed
					  homogenization */
  delete matrix;
  delete subspace;
  
  int num_rays = rays->get_number();
  
  /* Walk through all rays of the homogenization; each gives a vertex
     of the polytope. */
  int i;
  for (i = 0; i<num_rays; i++) {
    ZZ denominator = convert_mpz_to_ZZ((*rays)[i][num_ineqs]);
    if (denominator == 0) {
      /* This ray corresponds to a ray of the (unbounded!) polyhedron;
	 ignore it. */
      unbounded = true;
    }
    else {
      listCone *cone = createListCone();
      vec_ZZ numerator;
      numerator.SetLength(numOfVars);
      int j;
      for (j = 0; j<numOfVars; j++)
	numerator[j] = convert_mpz_to_ZZ((*rays)[i][num_ineqs + 1 + j]);
      rationalVector *vertex_vector = new rationalVector(numerator, denominator);
      cone->vertex = new Vertex(vertex_vector);

      /* Compute the facets: */
      /* Find incident facets.
	 They are the facets whose corresponding slack variables are
	 zero. */
      listVector *ineq;
      for (j = 0, ineq = ineqs; j<num_ineqs; j++, ineq = ineq->rest) {
	if ((*rays)[i][j] == 0) {
	  /* Incident! */
	  vec_ZZ facet_vector;
	  facet_vector.SetLength(numOfVars);
	  int k;
	  for (k = 0; k<numOfVars; k++)
	    facet_vector[k] = - ineq->first[k + 1];
	  /* Cancel GCD: */
	  ZZ gcd;
	  for (k = 0; k<numOfVars; k++)
	    gcd = GCD(gcd, facet_vector[k]);
	  if (gcd != 0 && gcd != 1) {
	    for (k = 0; k<numOfVars; k++)
	      facet_vector[k] /= gcd;
	  }
	  cone->facets = new listVector(facet_vector, cone->facets);
	}
      }

      /* FIXME: To cheaply compute the rays of the vertex cone, we need
	 to get hold of the adjacent rays of the current ray of the
	 homogenization.

	 For the moment, don't compute the rays.  Later code will construct
	 them when needed.
      */

      printCone(cone, numOfVars);
      
      consumer.ConsumeCone(cone);
    }
  }
  delete rays;
}

listCone *
computeVertexConesWith4ti2(listVector* matrix, int numOfVars,
			   bool &unbounded)
{
  CollectingConeConsumer ccc;
  computeVertexConesWith4ti2(matrix, numOfVars, unbounded, ccc);
  return ccc.Collected_Cones;
}
