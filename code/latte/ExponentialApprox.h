// This is a -*- C++ -*- header file.

#ifndef EXPONENTIALAPPROX_H
#define EXPONENTIALAPPROX_H

#include <fstream>
#include "myheader.h"
#include "latte_gmp.h"
#include "barvinok/dec.h"

using namespace std;

class Write_Exponential_Sample_Formula_Single_Cone_Parameters :
  public Exponential_Single_Cone_Parameters
{
public:
  ofstream stream;
  double sampling_factor;
  double total_lower_bound;
  double total_upper_bound;
  mpq_class approximate_result;
  Write_Exponential_Sample_Formula_Single_Cone_Parameters
  (const BarvinokParameters &params,
   const char *filename, double a_sampling_factor) :
    Exponential_Single_Cone_Parameters(params),
    stream(filename), sampling_factor(a_sampling_factor) {};
  virtual void InitializeComputation();
  virtual int ConsumeCone(listCone *cone);
};

void
decomposeAndWriteExponentialSampleFormula(listCone *cones,
					  Write_Exponential_Sample_Formula_Single_Cone_Parameters &param);

#endif
