/**************************************************************************
 Estimating the number of lattice points in the knapsack simplex { x : A x = B, x >= 0 }
 using an easy variant of what Persi presented at the Tables meeting in Davis:
   [Chen-Diaconis-Holmes-Liu: Sequential Monte Carlo Methods for
   Statistical Analysis of tables, Stanford Statistics Report, August 2003]

Author: Ruriko Yoshida
Date: October 15th, 2004.

Note: I modified the maple code from prof Sturmfels.

****************************************************************************/

#include <string.h>
#include <stdio.h>
#include <iostream>
#include <fstream>

#include <NTL/mat_ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_RR.h>
#include <NTL/vec_RR.h>
#include <NTL/config.h>
#include <NTL/LLL.h>
#include <NTL/HNF.h>
#include <NTL/ZZ.h>

/* #include "myheader.h"
#include "barvinok/dec.h"
#include "barvinok/barvinok.h"
#include "barvinok/cone.h"
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
#include "binarySearchIP.h" */

using namespace std;

/* ----------------------------------------------------------------- */
int main(int argc, char *argv[]) {
  int Dim, tmp, i, run, die, seed = 213, tmpbound;
  char command[127], filename[127];
  //if you are very picky...
  if(argc > 2){
    strcpy(filename, argv[2]);}
  else{
    cerr << "Enter the input file name: " << endl;
    cin >> filename;
    }

  ifstream in(filename);
  if(!in){
    cerr << "No file exits!" << endl;
    exit(1);
  }
  

  in >> tmp >> Dim;
  Dim--;
  seed = atoi(argv[1]);
  vec_ZZ A;
  A.SetLength(Dim);
  ZZ B, totalP, P, b, bound, answer;
  srand(seed);

  /* Reading the input. */
  in >> B;
  for(i = 0; i < Dim; i++)
    {in >> A[i]; A[i] = -A[i];}

  /* Estimating via MCMC. */
  for(run = 1; run < 20001; run++){
    P = 1; b = B;
    for(i = 1; i < Dim; i++){
      bound = b / A[Dim - i];
      conv(tmpbound, bound); 
      die = rand() % (tmpbound + 1);
      P = P * (bound + 1);
      b = b - die * A[Dim - i];
    }
    totalP = totalP + P;
  }

  /* Writing the output. */
  strcpy(command, filename);
  strcat(command, ".out");
  ofstream out(command);
  answer = totalP / run;
  out << "The estimated number of lattice points: " << answer << endl;
  out << GetTime() << " sec." << endl;
  cout << "The estimated number of lattice points: " << answer << endl;
  cout << GetTime() << " sec." << endl;
  return 0;
}
