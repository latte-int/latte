/* ExponentialEhrhart.cpp -- Computing Ehrhart polynomials
                             using the exponential substitution

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
#include "ExponentialSubst.h"
#include "ExponentialEhrhart.h"

void
Exponential_Ehrhart_Parameters::InitializeComputation()
{
  Generic_Vector_Single_Cone_Parameters::InitializeComputation();
  int i;
  for (i = 0; i<=Number_of_Variables; i++)
    ehrhart_coefficients[i] = 0;
}

int
Exponential_Ehrhart_Parameters::ConsumeCone(listCone *cone)
{
  assert(cone->rest == NULL);
  int status = 1;
  try {
    int numOfVars = Number_of_Variables;
    mpq_vector weights
      = computeExponentialResidueWeights(generic_vector, cone, numOfVars);
    mpz_vector sums_of_scalar_powers
      = compute_sums_of_scalar_powers_mpz(cone, numOfVars, generic_vector);
    ZZ ehrhart_vertex_scalar_zz;
    InnerProduct(ehrhart_vertex_scalar_zz, generic_vector,
		 cone->vertex->ehrhart_vertex);
    mpz_class ehrhart_vertex_scalar
      = convert_ZZ_to_mpz(ehrhart_vertex_scalar_zz);
    mpz_class ehrhart_vertex_scalar_power = 1;
    int n;
    for (n = 0; n<=numOfVars; n++) {
      int l;
      mpq_class x;
      for (l = n; l<=numOfVars; l++) {
	mpz_class binomial;
	mpz_bin_uiui(binomial.get_mpz_t(), l, n);
	x += binomial * weights[l] * sums_of_scalar_powers[l-n];
      }
      ehrhart_coefficients[n] += cone->coefficient * ehrhart_vertex_scalar_power * x;
      ehrhart_vertex_scalar_power *= ehrhart_vertex_scalar;
    }
  } catch (NotGenericException) {
    status = -1;
  }
  freeListCone(cone);
  return status;
}

mpq_vector 
decomposeAndComputeEhrhartPolynomial(listCone *cones,
				     Exponential_Ehrhart_Parameters &param)
{
  barvinokDecomposition_List(cones, param);
  return param.ehrhart_coefficients;
}
