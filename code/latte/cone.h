/* This is a -*- C++ -*- header file. */

#ifndef CONE_H
#define CONE_H

#include "rational.h"

struct listVector {
  vec_ZZ first;
  struct listVector *rest;
  listVector() : first(), rest(0) {}
  listVector(const vec_ZZ a_first) : first(a_first), rest(0) {}
};

typedef struct listCone {
  int coefficient;
  ZZ determinant;		// determinant of the matrix formed by
				// the RAYS, with sign
  rationalVector* vertex;
  listVector *rays;
  ZZ dual_determinant;		// determinant of the matrix formed by
				// the FACETS, with sign
  listVector *facets;
  // For simplicial cones where RAYS and FACETS are both computed, we
  // guarantee that < RAY_i, FACET_j > = -FACET_DIVISOR_i * DELTA_{i,j}.
  vec_ZZ facet_divisors;	
  listVector *latticePoints;
  vec_ZZ lattice_points_scalar_products; 
  struct listCone *rest;
} listCone;

listCone* createListCone();
int lengthListCone(listCone*);

/* Free the first cone, not the whole list. */
void freeCone(listCone *cone);

/* Free the whole list of cones. */
void freeListCone(listCone *list);

#endif
