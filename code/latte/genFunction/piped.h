// This is a -*- C++ -*- header file.

/* piped.h -- Computing all lattice points in the fundamental parallelepiped

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

#ifndef GENFUNCTION_PIPED_H
#define GENFUNCTION_PIPED_H

#include "cone.h"
#include "cone_consumer.h"
#include <vector>
using namespace std;

class PointsInParallelepipedGenerator {
protected:
  const listCone *cone;
  vec_ZZ max_multipliers;
  mat_ZZ B_inv;
  mat_ZZ U;
  vec_ZZ beta;
  ZZ facet_divisor_common_multiple;
  vec_ZZ facet_scale_factors;
  mat_ZZ Tau;
  vec_ZZ beta_mod_facet_divisors;
public:
  PointsInParallelepipedGenerator(const listCone *a_cone, int numOfVars);
  const vec_ZZ &GetMaxMultipliers();
  int *GetMaxMultipliers_int();
  /* Let n be the vector obtained from GetMaxMultipliers().
     Then all points in the fundamental parallelepiped can be obtained
     by calling GeneratePoint for all integer multiplier vectors m
     such that 0 <= m_i < n_i. */
  vec_ZZ GeneratePoint(const int *multipliers);
  vec_ZZ GeneratePoint(const vec_ZZ &multipliers);
protected:
  void compute_scaled_fundamental_multiplier_from_multipliers
  (ZZ &multiplier, const vec_ZZ &multipliers,
   const vec_ZZ &facet, int facet_index);
  void compute_scaled_fundamental_multiplier_from_multipliers
  (ZZ &multiplier, const int *multipliers,
   const vec_ZZ &facet, int facet_index);
};

listVector* pointsInParallelepiped(listCone *cone, int numOfVars);

/* Compute the latticePoints slot of CONE. */
void computePointsInParallelepiped(listCone *cone, int numOfVars);

/* For all cones in the linked list CONES, compute their latticePoints
   slot. */
void computePointsInParallelepipeds(listCone *cones, int numOfVars);

class PointsScalarProductsGenerator : public PointsInParallelepipedGenerator {
  vec_ZZ generic_vector;
  vec_ZZ scaled_ray_scalar_products;
public:
  PointsScalarProductsGenerator(const listCone *a_cone, int numOfVars,
				const vec_ZZ &a_generic_vector);
  ZZ GeneratePointScalarProduct(int *multipliers);
};

/* Compute the lattice_points_scalar_products slot of CONE. */
void computeLatticePointsScalarProducts(listCone *cone, int numOfVars,
					const vec_ZZ &generic_vector);

class PointsInParallelepipedComputingConeTransducer : public ConeTransducer
{
public:
  int ConsumeCone(listCone *cone);
};

#endif
