#include "ConeDeterminant.h"

int computeDeterminant(int**Mat, int n)
{
  mat_ZZ tmp;
  int D = 0;
  ZZ DD;
  tmp.SetDims(n, n);

  for(int i = 0; i < n; i++)
     for(int j = 0; j < n; j++)
        conv(tmp[i][j], Mat[i][j]);

  DD = determinant(tmp);
  if(DD > 1000000000){
    cerr << "Big Integer in ComputeDet!" << endl;
    exit(1);
   }

   conv(D, DD);
   return D;
}
