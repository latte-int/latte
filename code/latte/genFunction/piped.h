// This is a -*- C++ -*- header file.

#ifndef GENFUNCTION_PIPED_H
#define GENFUNCTION_PIPED_H

#include "cone.h"
#include <vector>
using namespace std;

class PointsInParallelepipedGenerator {
  const listCone *cone;
  vector<int> max_multipliers;
  mat_ZZ B_inv;
  mat_ZZ U;
  vec_ZZ beta;
  ZZ facet_divisor_common_multiple;
  vec_ZZ facet_scale_factors;
public:
  PointsInParallelepipedGenerator(const listCone *a_cone, int numOfVars);
  const vector<int> &GetMaxMultipliers();
  /* Let n be the vector obtained from GetMaxMultipliers().
     Then all points in the fundamental parallelepiped can be obtained
     by calling GeneratePoint for all integer multiplier vectors m
     such that 0 <= m_i < n_i. */
  vec_ZZ GeneratePoint(int *multipliers);
  vec_ZZ GeneratePoint(const vec_ZZ &multipliers);
private:
  vec_ZZ translate_lattice_point(const vec_ZZ& m);
};

listVector* pointsInParallelepiped(listCone *cone, int numOfVars);
listVector* pointsInParallelepipedOfUnimodularCone(rationalVector*, 
						   listVector*, int);

/* Compute the latticePoints slot of CONE. */
void computePointsInParallelepiped(listCone *cone, int numOfVars);

/* For all cones in the linked list CONES, compute their latticePoints
   slot. */
void computePointsInParallelepipeds(listCone *cones, int numOfVars);

#endif
