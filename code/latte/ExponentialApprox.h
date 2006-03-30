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
  double total_lower_bound;
  double total_upper_bound;
  Write_Exponential_Sample_Formula_Single_Cone_Parameters
  (const char *filename, int a_max_determinant) :
    Exponential_Single_Cone_Parameters(),
    stream(filename), total_lower_bound(0.0), total_upper_bound(0.0)
  {
    max_determinant = a_max_determinant;
  }
  virtual void InitializeComputation();
  virtual int ConsumeCone(listCone *cone);
};

void
decomposeAndWriteExponentialSampleFormula(listCone *cones,
					  Write_Exponential_Sample_Formula_Single_Cone_Parameters &param);

#endif
