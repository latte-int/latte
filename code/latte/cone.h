/* This is a -*- C++ -*- header file. */

/* cone.h -- Linked list of cones

   Copyright 2002-2004 Jesus A. De Loera, David Haws, Raymond
      Hemmecke, Peter Huggins, Jeremy Tauzer, Ruriko Yoshida
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

#ifndef CONE_H
#define CONE_H

#include "rational.h"

struct listVector {
  vec_ZZ first;
  struct listVector *rest;
  listVector() : first(), rest(0) {}
  listVector(const vec_ZZ &a_first, struct listVector *a_rest = 0) :
    first(a_first), rest(a_rest) {}
};

/* Return the length of the linked list of vectors. */
int lengthListVector(const listVector*);

listVector *appendVectorToListVector(const vec_ZZ &, listVector*);

/* Copy the linked list of vectors. */
listVector *copyListVector(listVector *);

/* Free the whole linked list of vectors. */
void freeListVector(listVector *p);

struct Vertex {
  rationalVector *vertex;
  vec_ZZ ehrhart_vertex; // for exponential Ehrhart computation
  Vertex(rationalVector *v) : vertex(v) {}
  Vertex(const Vertex &v) : vertex(new rationalVector(*v.vertex)),
			    ehrhart_vertex(v.ehrhart_vertex) {}
  ~Vertex() { delete vertex; }
};

struct listCone {
  int coefficient;
  Vertex *vertex;
  ZZ determinant;		// determinant of the matrix formed by
				// the RAYS, with sign
  listVector *rays;
  listVector *subspace_generators; // for non-pointed cones
  ZZ dual_determinant;		// determinant of the matrix formed by
				// the FACETS, with sign
  listVector *facets;
  listVector *equalities;	// for non-fulldimensional cones
  // For simplicial cones where RAYS and FACETS are both computed, we
  // guarantee that < RAY_i, FACET_j > = -FACET_DIVISOR_i * DELTA_{i,j}.
  vec_ZZ facet_divisors;	
  listVector *latticePoints;
  vec_ZZ lattice_points_scalar_products;
  struct listCone *rest;
};

/* Allocate a single listCone element and initialize all members. */
listCone* createListCone();

/* Return the length of the linked list. */
int lengthListCone(listCone*);

/* Free the first cone, not the whole list. */
void freeCone(listCone *cone);

/* Free the whole list of cones. */
void freeListCone(listCone *list);

/* Destructively concatenate two linked lists of cones.  Return the
   resulting list. */
listCone *appendListCones(listCone *A, listCone *B);

/* Copy a single cone. */
listCone *copyCone(listCone *cone);

/* Deduce the ambient dimension of CONE from its data. */
int ambient_cone_dimension(listCone *cone);

#endif
