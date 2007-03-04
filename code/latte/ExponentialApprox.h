// This is a -*- C++ -*- header file.

/* ExponentialApprox.h -- Sample lattice points in fundamental
   parallelepipeds, using the exponential substitution

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

#ifndef EXPONENTIALAPPROX_H
#define EXPONENTIALAPPROX_H

#include <fstream>
#include "latte_gmp.h"
#include "barvinok/dec.h"
#include "heap.h"

using namespace std;

class Write_Exponential_Sample_Formula_Single_Cone_Parameters;

class ConeApproximationData {
public:
  listCone *cone;
  mpq_vector weights;
  vec_ZZ ray_scalar_products;
  mpq_class this_total_lower_bound;
  mpq_class this_total_upper_bound;
  
  ConeApproximationData(listCone *a_cone, const Write_Exponential_Sample_Formula_Single_Cone_Parameters &params);
  double GetWeight();
};
  
class Write_Exponential_Sample_Formula_Single_Cone_Parameters :
  public Exponential_Single_Cone_Parameters
{
public:
  ofstream stream;
  double sampling_factor;
  long int num_samples;
  mpq_class total_lower_bound;
  mpq_class total_upper_bound;
  mpq_class approximate_result;
  heap *cone_heap;
  Write_Exponential_Sample_Formula_Single_Cone_Parameters
  (const BarvinokParameters &params,
   const char *filename, double a_sampling_factor, long int a_num_samples) :
    Exponential_Single_Cone_Parameters(params),
    stream(filename), sampling_factor(a_sampling_factor), num_samples(a_num_samples) {};
  virtual void InitializeComputation();
  virtual int ConsumeCone(listCone *cone);
public:
  void ShowStats();
  double EvaluateCone(listCone *cone);
};

void
decomposeAndWriteExponentialSampleFormula(listCone *cones,
					  Write_Exponential_Sample_Formula_Single_Cone_Parameters &param);

#endif
