
#include <string.h>
#include <stdio.h>
#include <NTL/mat_ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_RR.h>
#include <NTL/vec_RR.h>
#include <NTL/config.h>
#include <NTL/LLL.h>
#include <NTL/HNF.h>
#include <NTL/ZZ.h>

#include "myheader.h"
#include "barvinok/dec.h"
#include "barvinok/barvinok.h"
#include "barvinok/Cone.h"
#include "barvinok/ConeDecom.h"
#include "barvinok/Triangulation.h"
#include "vertices/cdd.h"
#include "genFunction/maple.h"
#include "genFunction/piped.h"
#include "cone.h"
#include "ConeDeterminant.h"
#include "dual.h"
#include "RudyResNTL.h"
#include "Grobner.h"
//  #include "jesus.h"
#include "preprocess.h"
#include "print.h"
#include "ramon.h"
#include "rational.h"
#include "timing.h"
#include "flags.h"
//#include "testing.h"
#include "IntegralHull.h"
#include "ReadingFile.h"



/* ----------------------------------------------------------------- */
int main(int argc, char *argv[]) {
  mat_ZZ A, W, AA;
  int m, n;
  ifstream in(argv[1]);
  in >> m >> n;
  A.SetDims(m, n);
  AA.SetDims(m, m);
  for(int i = 0; i < m; i++)
    for(int j = 0; j < n; j++) in >> A[i][j];
  /* for(int i = 0; i < m; i++)
     for(int j = 0; j < m; j++) AA[i][j] = A[i][j]; */
  cout << A << endl;
  vec_ZZ b;
  b.SetLength(m);

  ZZ D;
  //  A=transpose(A);
 D = 3; 
 //   HNF(W, AA, D);
   //LLL(D, A, W);

 A=inv(A);
 cout << b*A << endl;
 //  cout << A << "\n" << W << endl;
 return(0);
}
/* ----------------------------------------------------------------- */




