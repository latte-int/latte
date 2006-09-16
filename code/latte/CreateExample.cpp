/* CreateExample.cpp -- This is for creating examples of IP.

   Copyright 2003 Ruriko Yoshida

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

#include <list>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <math.h>
#include <algorithm>
#include <time.h>

#include <NTL/ZZ.h>
#include <NTL/mat_ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/RR.h>
#include <NTL/mat_RR.h>
#include <NTL/vec_RR.h>
#include <NTL/LLL.h>

using namespace std;

int main(int argc, char *argv[]){
  int dim, seed;
  dim = atoi(argv[1]);
  seed = atoi(argv[2]);
  ZZ s;
  conv(s, seed);
  srand(seed);
  vec_RR OriginalVec;
  mat_RR W;
  mat_ZZ randomVec, Kernel;
  randomVec.SetDims(dim, 1);
  W.SetDims(dim - 1, dim);
  OriginalVec.SetLength(dim);
  mat_RR M;
  M.SetDims(dim - 1, dim - 1);
  //SetSeed(s);

  RR detM, detWW;

  for(int i = 0; i < dim; i++) {randomVec[i][0]=rand()%1000000;
  conv(OriginalVec[i], randomVec[i][0]);}
  
  for(int i = 0; i < dim - 1; i++)
    for(int j = 0; j < dim - 1; j++) M[i][j] = 1;

  for(int i = 0; i < dim - 1; i++) conv(M[i][i], (OriginalVec[dim - 1]/OriginalVec[i]) * (OriginalVec[dim - 1]/OriginalVec[i]) + 1);

  ZZ det;
  cout << randomVec << endl; 
  LLL(det, randomVec, Kernel);
    
  //   cout << Kernel << endl;
  for(int i = 0; i < dim; i++) 
    for(int j = 0; j < dim - 1; j++) conv(W[j][i],  Kernel[j][i]);
  
  detM = determinant(M);
  detWW = determinant(W * transpose(W));

  RR c1, c2, V;
  V = power(OriginalVec[dim - 1], dim - 1);
  // cout << detM <<" " << M << endl;
  c1 = sqrt(detWW)/2;
  c2 = sqrt(detM);
  RR fact;
  fact = 1;
  for(int i = 2; i < dim; i++) fact = fact * i;
  c2 = c2 / fact;
  c2 = c2 / V;
  //  cout << c1 << " " << c2 << endl;
  RR num;
  RR P;
  P = dim - 1;
  P = inv(P);
  num = pow(c1/c2, P);
  cout << num << endl; 

  ZZ L, U;
  L = CeilToZZ(3 * num / 4);
  U = CeilToZZ(5 * num / 4);
  cout << L << " " << U << endl;

  return 0;
}
