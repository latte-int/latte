/// Experimental code.

#include <cassert>
#include "dual.h"
#include "print.h"
#include "ExponentialSubst.h"
#include "ExponentialApprox.h"

void Write_Exponential_Sample_Formula_Single_Cone_Parameters::InitializeComputation()
{
  Generic_Vector_Single_Cone_Parameters::InitializeComputation();
  stream << "*** Computation with new generic vector " << endl;
}

int
Write_Exponential_Sample_Formula_Single_Cone_Parameters::ConsumeCone(listCone *cone)
{
  assert(cone->rest == NULL);

  printConeToFile(stream, cone, Number_of_Variables);
  try {
    mpz_class prod;
    mpq_vector weights
      = computeExponentialResidueWeights(generic_vector,
					 prod,
					 cone,
					 Number_of_Variables);
    int dimension = weights.size();
    stream << "Approximate Weights: ";
    double d_root_prod = pow(abs(prod.get_d()), 1.0/dimension);
    mpq_vector::const_iterator i;
    int k;
    for (k = 0, i = weights.begin(); i!=weights.end(); ++i, k++) {
      double weight = (*i).get_d();
      // scale
      double scaled_weight = weight * pow(d_root_prod, k);
      stream << scaled_weight << " ";
    }
    stream << endl << endl;
    return 1;
  }
  catch (NotGenericException) {
    return -1;
  }
}

void
decomposeAndWriteExponentialSampleFormula(listCone *cones,
					  Write_Exponential_Sample_Formula_Single_Cone_Parameters &param)
{
  barvinokDecomposition_List(cones, param);
}
