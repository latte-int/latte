// This is a -*- C++ -*- header file.

/* ExponentialSubst.h -- Computing Ehrhart polynomials
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

#ifndef EXPONENTIALEHRHART_H
#define EXPONENTIALEHRHART_H

#include "latte_gmp.h"
#include "barvinok/dec.h"

class Exponential_Ehrhart_Parameters
  : public Generic_Vector_Single_Cone_Parameters {
public:
  mpq_vector ehrhart_coefficients;
  Exponential_Ehrhart_Parameters(const BarvinokParameters &params) :
    Generic_Vector_Single_Cone_Parameters(params),
    ehrhart_coefficients(params.Number_of_Variables + 1) {};
  virtual void InitializeComputation();
  virtual int ConsumeCone(listCone *cone);
};

mpq_vector 
decomposeAndComputeEhrhartPolynomial(listCone *cones,
				     Exponential_Ehrhart_Parameters &param);

#endif
