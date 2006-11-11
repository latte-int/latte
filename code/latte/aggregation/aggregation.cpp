/*********************************************************************
 Author:  Ruriko Yoshida
 Date:    February 21th, 2005
 Update:  February 22nd, 2005

 Program: This is for aggregation transformation of a system of linear 
 equations. 
 This is from Theorem 2.1 on "Polynomial-Time Aggregation
 of Integer Programming Problems" by RAVINDRAN KANNAN, Journal of the ACM
  Volume 30 ,  Issue 1,  Pages: 133 - 145  (1983).

***********************************************************************/

#include <fstream.h>
#include <iostream.h>
#include <cstdlib>
#include <cstring>
#include <math.h>
#include <algorithm>
#include <time.h>
#include <string>

#include <NTL/ZZ.h>
#include <NTL/mat_ZZ.h>
#include <NTL/vec_ZZ.h>

#include "matrix.h"
#include "read.h"

using namespace std;

// ZZ power(int x, int exp){
//   ZZ answer, tmp;
//   conv(answer, 1);
//   conv(tmp, x);
//   for(int i = 1; i <= exp; i++){ answer *= tmp;}
//   return answer;
// }

int main(int argc, char **argv)
{
  int n, m, *RHS, bigB = 0, bigA = 0;
  Matrix Mat;
  ZZ newRHS, gcd, lambda;
  vec_ZZ newMat;
  char file[127];

  if(argc < 2) {
    cerr << "Please enter the input file." << endl;
    exit(0);
  }

  strcpy(file, argv[1]);
  strcat(file, ".lat"); 
  ReadInput(argv[1], &Mat, &RHS, n, m);

  mat_ZZ A;
  A.SetDims(m, n);
  // cout << Mat.numOfRows << " " << Mat.numOfCols << endl;
//   for(int i = 0; i < m; i++){cout << RHS[i] << " ";
//     for(int j = 0; j < n; j++)
//       cout << -Mat.elements[i][j] << " ";
//     cout << endl;
//   }

  newMat.SetLength(n);
 
   cout << m << " " << n << endl;
  for(int i = 0; i < m; i++)
    for(int j = 0; j < n; j++) {
      conv(A[i][j], Mat.elements[i][j]);
      if(bigA < Mat.elements[i][j]) bigA = Mat.elements[i][j];
    }

  for(int i = 0; i < m; i++){cout << RHS[i] << " ";
    for(int j = 0; j < n; j++)
      cout << -A[i][j] << " ";
    cout << endl;
  }

  for(int i = 0; i < m; i++) 
    if(bigB < RHS[i]) bigB = RHS[i];

  if(bigB == 0){
    cerr << "Error!  The RHS must contains nonnegative elements!" << endl;
    exit(0);
  }

  if(bigA == 0){
    cerr << "Error!  The LHS must contains nonnegative elements!" << endl;
    exit(0);
  }
  //cout << m << " " << n << " " << bigA << " " << bigB << endl;
  conv(lambda, 3 * ((n * m) * (bigA * bigB)));
  //  lambda = 3 * ((n * m) * (bigA * bigB));
  
  for(int j = 0; j < n; j++){
    for(int i = 0; i < m; i++){ 
      newMat[j] += (power(lambda, m + 1) + power(lambda, i + 1)) * A[i][j];
      }
  }
  
  for(int i = 0; i < m; i++) 
    newRHS += (power(lambda, m + 1) + power(lambda, i + 1)) * RHS[i];
  
  gcd = GCD(newRHS,  newMat[0]);

  for(int j = 1; j < n; j++) 
    if(gcd > GCD(newRHS,  newMat[j]))
      gcd = GCD(newRHS,  newMat[j]);
  
  //cout << lambda << " " <<gcd <<endl;

  for(int j = 0; j < n; j++) newMat[j] /= gcd;
  newRHS /= gcd;

  ofstream out(file);
  if(!out) {
    cerr << "Error on writing the latte output." << endl;
    exit(0);
  }
  out << 1 << " " << n + 1 << endl;
  out << newRHS << " ";
  for(int j = 0; j < n; j++) out << -newMat[j] << " ";
  out << endl;
  out << "linearity 1 1" << endl;
  out << "nonnegative " << n << " ";
  for(int i = 0; i < n; i++) out << i + 1 << " ";
  out << endl;

  return 0;
}
