// This is a -*- C++ -*- header file.

#ifndef EXPONENTIALAPPROX_H
#define EXPONENTIALAPPROX_H

#include <fstream>
#include "myheader.h"
#include "latte_gmp.h"
#include "barvinok/dec.h"

using namespace std;

class Write_Exponential_Sample_Formula_Single_Cone_Parameters :
  public Generic_Vector_Single_Cone_Parameters
{
public:
  ofstream stream;
  Write_Exponential_Sample_Formula_Single_Cone_Parameters
  (const char *filename, int a_max_determinant) :
    stream(filename)
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
