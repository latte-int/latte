/* piped.cpp -- Computing all lattice points in the fundamental parallelepiped

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

#include "config.h"
#include "../cone.h"
#include "../print.h"
#include "../ramon.h"
#include "../rational.h"
#include "../vertices/cdd.h"
#include <cassert>
#include <string>

#include "IntCombEnum.h"
#include "matrix_ops.h" 
#include <NTL/mat_ZZ.h>
#include "convert.h"
#include "piped.h"

using namespace std;

/*
 * Description:
 * pointsInParallelepiped: 1) Find Smith Normal form of cone basis. In other
 * words, find a basis v and w such that v_i = n_i*w_i. 2) Enumerate all
 * lattice points in fund||piped by taking all (bounded: 0 <= k <= n_{i-1})
 * integer combinations of w_i, and translating them accordingly
 *
 * Parameters:
 * listCone *cone : data pertaining to the cone
 * int numOfVars: dimension of cone
 *
 * Return: listVector* of lattice points
 */

void
PointsInParallelepipedGenerator::compute_scaled_fundamental_multiplier_from_multipliers
(ZZ &scaled_fundamental_multiplier, const vec_ZZ &multipliers,
 const vec_ZZ &facet, int facet_index)
{
  int dim = facet.length();
  scaled_fundamental_multiplier = beta_mod_facet_divisors[facet_index];
  ZZ t;
  int j;
  for (j = 0; j<dim; j++) {
    mul(t, Tau[facet_index][j], multipliers[j]);
    scaled_fundamental_multiplier -= t;
  }
  scaled_fundamental_multiplier %= cone->facet_divisors[facet_index];
  scaled_fundamental_multiplier -= beta[facet_index];
}

void
PointsInParallelepipedGenerator::compute_scaled_fundamental_multiplier_from_multipliers
(ZZ &scaled_fundamental_multiplier, const int *multipliers,
 const vec_ZZ &facet, int facet_index)
{
  int dim = facet.length();
  scaled_fundamental_multiplier = beta_mod_facet_divisors[facet_index];
  ZZ t;
  int j;
  for (j = 0; j<dim; j++) {
    mul(t, Tau[facet_index][j], multipliers[j]);
    scaled_fundamental_multiplier -= t;
  }
  scaled_fundamental_multiplier %= cone->facet_divisors[facet_index];
  scaled_fundamental_multiplier -= beta[facet_index];
}

/*
 * If S is the Smith Normal Form of A, then the multipliers n_i such that
 * v_i = n_i w_i are on the diagonal
 */
static vec_ZZ
get_multipliers_from_snf(const mat_ZZ & snf)
{
   int j;
   int row = snf.NumRows();
   vec_ZZ n;
   n.SetLength(row);
   for (j = 0; j < row; j++) {
     //cout << "get_multipliers_from_snf:: snf[" << j << "," << j << "] = " << snf[j][j] << "\n"; 
     n[j] = snf[j][j];
   }
   return n;
}

PointsInParallelepipedGenerator::PointsInParallelepipedGenerator(const listCone *a_cone, int numOfVars) :
  cone(a_cone)
{
  U = convert_listVector_to_mat_ZZ(cone->rays);

  if (abs(cone->determinant) == 1) {
    /* Unimodular case: Id = Smith(U) = B U C.
       Thus, for instance, B = det U * U^{-1}, C = Id.
     */
    max_multipliers.SetLength(numOfVars);
    int i;
    for (i = 0; i<numOfVars; i++)
      max_multipliers[i] = 1;
    B_inv = cone->determinant * U;
  }
  else
    {
    /* General case: */
    mat_ZZ snf_U, B, C;

    //cout << "Computing Smith-Normal form...\n";
    /* get Smith Normal form of matrix, Smith(U) = BUC */
    snf_U = SmithNormalForm(U, B, C);

#ifdef CHECK_SNF
    {
      mat_ZZ SU = B * U * C;
      assert(SU == snf_U);
      int i, j;
      for (i = 0; i<numOfVars; i++) {
	for (j = 0; j<numOfVars; j++)
	  if (i != j) assert(IsZero(snf_U[i][j]));
      }
      assert(abs(determinant(B)) == 1);
      assert(abs(determinant(C)) == 1);
    }
#endif
    /* extract n_i such that v_i = n_i*w_i from Smith Normal form */
    max_multipliers = get_multipliers_from_snf(snf_U);

    /*
     * Smith(U) = BUC, and the diagonal entries of Smith(U) are the multipliers
     * n_i such that w_i*n_i= v_i. Therefore, UC = V, and B^-1 = W. Since we
     * must take the integer combinations of W, we must calculate B^-1.
     */
#if 0
    /* However, this is expensive. --mkoeppe */
    B_inv = inv(B);
#else
    /* We also have B^-1 = UC Smith(U)^{-1}; this is faster to
       compute because Smith(U) is diagonal. --mkoeppe */
    B_inv = U * C;
    int i, j;
    ZZ q, r;
    for (i = 0; i<B_inv.NumRows(); i++) {
      for (j = 0; j<B_inv.NumCols(); j++) {
	DivRem(q, r, B_inv[i][j], snf_U[j][j]);
	assert(IsZero(r));
	B_inv[i][j] = q;
      }
    }
#endif
  }
  /* We compute beta_i = floor(<v, facet_i>) and beta_i modulo facet_divisor_i. */
  {
    ZZ v_scale_factor;
    vec_ZZ v_scaled = scaleRationalVectorToInteger(cone->vertex,
						   numOfVars, v_scale_factor);
    beta.SetLength(numOfVars);
    beta_mod_facet_divisors.SetLength(numOfVars);
    int i;
    listVector *facet;
    ZZ sp;
    for (i = 0, facet = cone->facets; i<numOfVars; i++, facet=facet->rest) {
      InnerProduct(sp, v_scaled, facet->first);
      div(beta[i], sp, v_scale_factor);
      /* assert(beta[i] * v_scale_factor <= sp); */
      beta_mod_facet_divisors[i] = beta[i] % cone->facet_divisors[i];
    }
  }
  /* Precompute facet_scale_factors. */
  facet_scale_factors.SetLength(numOfVars);
  facet_divisor_common_multiple = abs(cone->determinant);
  {
    int i;
    for (i = 0; i<numOfVars; i++)
      facet_scale_factors[i] = facet_divisor_common_multiple / cone->facet_divisors[i];
  }
  /* Precompute Tau_{i,j} = < j-th column of B_inv, facet_i > modulo facet_divisor_i. */
  Tau.SetDims(numOfVars, numOfVars);
  {
    int i;
    listVector *facet;
    for (i = 0, facet = cone->facets; i<numOfVars; i++, facet=facet->rest) {
      int j;
      for (j = 0; j<numOfVars; j++) {
	ZZ tau;
	/* Scalar multiplication modulo facet_divisors[i]. */
	int k;
	for (k = 0; k<numOfVars; k++) {
	  tau += B_inv[k][j] * facet->first[k];
	}
	tau %= cone->facet_divisors[i];
	Tau[i][j] = tau;
      }
    }
  }
}

const vec_ZZ &
PointsInParallelepipedGenerator::GetMaxMultipliers()
{
  return max_multipliers;
}

int *
PointsInParallelepipedGenerator::GetMaxMultipliers_int()
{
  int length = max_multipliers.length();
  int *n = new int[length];
  int i;
  for (i = 0; i<length; i++) {
     if (max_multipliers[i] > INT_MAX) {
       cerr << "Implementation restriction hit:  Smith normal form has entries larger than INT_MAX\n";
       abort();
     }
     n[i] = to_int(max_multipliers[i]);
  }
  return n;
}

vec_ZZ
PointsInParallelepipedGenerator::GeneratePoint(const int *multipliers)
{
  int dimension = max_multipliers.length();
  int i;
  listVector *facet, *ray;
  vec_ZZ result;
  result.SetLength(dimension);
  for (i = 0, facet = cone->facets, ray = cone->rays;
       i<dimension;
       i++, facet=facet->rest, ray=ray->rest) {
    ZZ scaled_fundamental_multiplier;
    compute_scaled_fundamental_multiplier_from_multipliers
      (scaled_fundamental_multiplier, multipliers, facet->first, i);
    result += scaled_fundamental_multiplier * facet_scale_factors[i] * ray->first;
  }
  for (i = 0; i<dimension; i++)
    result[i] /= facet_divisor_common_multiple;
  return result;
}

vec_ZZ
PointsInParallelepipedGenerator::GeneratePoint(const vec_ZZ &multipliers)
{
  int dimension = max_multipliers.length();
  int i;
  listVector *facet, *ray;
  vec_ZZ result;
  result.SetLength(dimension);
  for (i = 0, facet = cone->facets, ray = cone->rays;
       i<dimension;
       i++, facet=facet->rest, ray=ray->rest) {
    ZZ scaled_fundamental_multiplier;
    compute_scaled_fundamental_multiplier_from_multipliers
      (scaled_fundamental_multiplier, multipliers, facet->first, i);
    result += scaled_fundamental_multiplier * facet_scale_factors[i] * ray->first;
  }
  for (i = 0; i<dimension; i++)
    result[i] /= facet_divisor_common_multiple;
  return result;
}

static void
check_point_in_parallelepiped(listCone *cone, const vec_ZZ &point)
{
  int numOfVars = point.length();
  ZZ v_scale_factor;
  vec_ZZ v_scaled = scaleRationalVectorToInteger(cone->vertex,
						 numOfVars, v_scale_factor);
  listVector *facet, *ray;
  for (facet = cone->facets, ray = cone->rays;
       facet;
       facet = facet->rest, ray = ray->rest) {
    vec_ZZ &f = facet->first;
    vec_ZZ &r = ray->first;
    ZZ psp;
    InnerProduct(psp, f, point);
    psp *= v_scale_factor;
    ZZ vsp;
    InnerProduct(vsp, f, v_scaled);
    ZZ rsp;
    InnerProduct(rsp, f, r);
    rsp *= v_scale_factor;
    assert(psp <= vsp);
    assert(psp > vsp + rsp);
  }
}

static bool
are_points_unique(listVector *points)
{
  listVector *a, *b;
  for (a = points; a; a = a->rest) {
    for (b = a->rest; b; b = b->rest) {
      if (a->first == b->first) return false;
    }
  }
  return true;
}

listVector*
pointsInParallelepiped(listCone *cone, int numOfVars)
{
  PointsInParallelepipedGenerator generator(cone, numOfVars);
  int *n = generator.GetMaxMultipliers_int();
  /*
   * enumerate lattice points by taking all integer combinations
   * 0 <= k <= (n_i - 1) of each vector.
   */
  IntCombEnum iter_comb(n, numOfVars);
  //cout << "Enumerating lattice points...\n";
  iter_comb.decrementUpperBound();
  listVector *lat_points = NULL;
  int *next;
  while((next = iter_comb.getNext())) {
    vec_ZZ trans_lat_pt = generator.GeneratePoint(next);
    lat_points = appendVectorToListVector(trans_lat_pt, lat_points);
  }
  delete[] n;
  return (lat_points);
}

/* ----------------------------------------------------------------- */

void computePointsInParallelepiped(listCone *cone, int numOfVars)
{
#if 0
  if (abs(cone->determinant) != 1)
    cout << "Processing cone with determinant " << cone->determinant << endl;
#endif
  cone->latticePoints = pointsInParallelepiped(cone, numOfVars);
}

void computePointsInParallelepipeds(listCone *cones, int numOfVars)
{
  listCone *tmp = cones;
  int Cones_Processed_Count = 0;
  while (tmp) {
    computePointsInParallelepiped(tmp, numOfVars);
    tmp=tmp->rest;
    Cones_Processed_Count++;
    if ((Cones_Processed_Count % 1000) == 0 )
      cout << Cones_Processed_Count << " cones processed." << endl;
  }
}

/* ----------------------------------------------------------------- */

PointsScalarProductsGenerator::PointsScalarProductsGenerator
(const listCone *a_cone, int numOfVars, const vec_ZZ &a_generic_vector) :
  PointsInParallelepipedGenerator(a_cone, numOfVars),
  generic_vector(a_generic_vector)
{
  scaled_ray_scalar_products.SetLength(numOfVars);
  int i;
  ZZ inner;
  listVector *ray;
  for (i = 0, ray = cone->rays; i<numOfVars; i++, ray=ray->rest) {
    InnerProduct(inner, generic_vector, ray->first);
    scaled_ray_scalar_products[i] = facet_scale_factors[i] * inner;
  }
}

ZZ PointsScalarProductsGenerator::GeneratePointScalarProduct(int *multipliers)
{
  ZZ result;
  result = 0;
  int dim = beta.length();
  int i;
  listVector *facet;
  listVector *ray;
  ZZ multiplier;
  for (i = 0, facet = cone->facets, ray = cone->rays;
       i<dim;
       i++, facet=facet->rest, ray=ray->rest) {
    compute_scaled_fundamental_multiplier_from_multipliers(multiplier, multipliers, facet->first, i);
    result += multiplier * scaled_ray_scalar_products[i];
  }
  ZZ q, r;
  DivRem(q, r, result, facet_divisor_common_multiple);
  assert(IsZero(r));
  return q;
}

void computeLatticePointsScalarProducts(listCone *cone, int numOfVars,
					const vec_ZZ &generic_vector)
{
  ZZ index = abs(cone->determinant);
  if (index > INT_MAX) {
    cerr << "Implementation restriction hit:  Attempt to enumerate a fundamental parallelepiped of index greater than INT_MAX.  (Probably not a good idea anyway.)" << endl;
    abort();
  }
  int num_lattice_points = to_long(index);
  cone->lattice_points_scalar_products.SetLength(num_lattice_points);
#if 0
  /* Slow method... */
  computePointsInParallelepiped(cone, numOfVars);
#endif
  if (cone->latticePoints) {
    // Lattice points already computed.
    listVector *point;
    int i;
    for (point = cone->latticePoints, i = 0; point != NULL; point = point->rest, i++) {
      InnerProduct(cone->lattice_points_scalar_products[i],
		   generic_vector, point->first);
    }
  }
  else {
    PointsScalarProductsGenerator generator(cone, numOfVars, generic_vector);
    int *n = generator.GetMaxMultipliers_int();
    /*
     * enumerate lattice points by taking all integer combinations
     * 0 <= k <= (n_i - 1) of each vector.
     */
    IntCombEnum iter_comb(n, numOfVars);
    //cout << "Enumerating lattice points...\n";
    iter_comb.decrementUpperBound();
    int *next;
    int num_scalar = 0;
    while((next = iter_comb.getNext())) {
      cone->lattice_points_scalar_products[num_scalar]
	= generator.GeneratePointScalarProduct(next);
      num_scalar++;
    }
    delete[] n;
  }
}
