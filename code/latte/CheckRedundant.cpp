
/*                                                                    *
 * Author: Ruriko Yoshida                                             *
 * Date: October 25th, 2003                                           *
 * Update: October 25th, 2003                                         *
 * This is for checking redundant inequalities and hidden equations.  *
 *
 */

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
#include "binarySearchIP.h"

listVector* CheckRedIneq(mat_ZZ S){
  mat_ZZ SS;
  SS = S;

  listVector *answer, *endanswer;
  int numOfRows, numOfColms, NewNumOfRows = 0;
  numOfRows = S.NumRows();
  numOfColms = S.NumCols();
  ZZ rhs;
  RR optVal, Rat_rhs;
  vec_ZZ cost;
  vec_RR Rat_solution, tmp_num, tmp_den;

  cost.SetLength(numOfColms - 1);
  tmp_num.SetLength(numOfColms - 1);
  tmp_den.SetLength(numOfColms - 1);
  Rat_solution.SetLength(numOfColms - 1);

  for(int i = 0; i < numOfRows; i++)
    SS[i][0] = S[i][0] + 1;

  answer = createListVector(S[0]);
  endanswer = answer;
  rationalVector* LP_vertex;
  vec_RR Rat_cost;  Rat_cost.SetLength(numOfColms - 1);

  for(int i = 0; i < numOfRows; i++){
    listVector *basis, *endBasis;
    basis = createListVector(S[0]);
    endBasis = basis;
    if(i == 0) endBasis -> first[0]++;
    for(int j = 1; j < numOfRows; j++){
      endBasis->rest = createListVector(S[j]);
      endBasis=endBasis->rest;
      if(j == i) endBasis -> first[0]++;
    }
    rhs = S[i][0];
    for(int j = 0; j < numOfColms - 1; j++) cost[j] = -S[i][j + 1]; 
    //printListVector(basis, numOfColms);// cout << endl << cost << endl;
    LP_vertex = LP(basis, cost, numOfColms - 1);
    // cout << LP_vertex->enumerator << endl << LP_vertex->denominator << endl;
    for (int s = 0; s < numOfColms - 1; s++){
      conv(tmp_num[s], LP_vertex->enumerator[s]);
      conv(tmp_den[s], LP_vertex->denominator[s]);
      Rat_solution[s] = tmp_num[s]/tmp_den[s];
      conv(Rat_cost[s], cost[s]);
    }
    optVal = Rat_solution * Rat_cost;
    conv(Rat_rhs, rhs); //cout << optVal << " " << rhs << " " << Rat_rhs << endl; 
    if(optVal > Rat_rhs)  {
      NewNumOfRows++;
      endanswer -> rest = createListVector(S[i]);
      endanswer = endanswer  -> rest;}
    else ;

  }

  return (answer ->rest);
}


listVector* CheckHidEqs(mat_ZZ S, listVector* equ){
  mat_ZZ SS;
  SS = S;

  listVector *answer, *endanswer, *endEqu;
  int numOfRows, numOfColms;
  numOfRows = S.NumRows();
  numOfColms = S.NumCols();
  ZZ rhs;
  RR optVal, optVal2, Rat_rhs;
  vec_ZZ cost;
  vec_RR Rat_solution, tmp_num, tmp_den;
  ZZ OptVal, OptVal2;

  cost.SetLength(numOfColms - 1);
  tmp_num.SetLength(numOfColms - 1);
  tmp_den.SetLength(numOfColms - 1);
  Rat_solution.SetLength(numOfColms - 1);

  for(int i = 0; i < numOfRows; i++)
    SS[i][0] = S[i][0] + 1;

  answer = createListVector(S[0]);
  endanswer = answer;
  equ = createListVector(S[0]);
  endEqu = equ;

  rationalVector* LP_vertex;
  vec_RR Rat_cost;  Rat_cost.SetLength(numOfColms - 1);

  for(int i = 0; i < numOfRows; i++){
    listVector *basis, *endBasis;
    basis = createListVector(S[0]);
    endBasis = basis;
    for(int j = 1; j < numOfRows; j++){
      endBasis->rest = createListVector(S[j]);
      endBasis=endBasis->rest;
    }
    rhs = S[i][0];
    for(int j = 0; j < numOfColms - 1; j++) cost[j] = -S[i][j + 1]; 
    //   printListVector(basis, numOfColms); //cout << endl << cost << endl;
    LP_vertex = LP(basis, cost, numOfColms - 1);
    for (int s = 0; s < numOfColms - 1; s++){
      conv(tmp_num[s], LP_vertex->enumerator[s]);
      conv(tmp_den[s], LP_vertex->denominator[s]);
      Rat_solution[s] = tmp_num[s]/tmp_den[s];
      conv(Rat_cost[s], cost[s]);
    }
    optVal = Rat_solution * Rat_cost;
    conv(OptVal, optVal);

    for(int j = 0; j < numOfColms - 1; j++) cost[j] = S[i][j + 1];
    // printListVector(basis, numOfColms); //cout << endl << cost << endl;
    LP_vertex = LP(basis, cost, numOfColms - 1);
    //cout << LP_vertex->enumerator << endl << LP_vertex->denominator << endl;
    for (int s = 0; s < numOfColms - 1; s++){
      conv(tmp_num[s], LP_vertex->enumerator[s]);
      conv(tmp_den[s], LP_vertex->denominator[s]);
      Rat_solution[s] = tmp_num[s]/tmp_den[s];
      conv(Rat_cost[s], cost[s]);
    }
    optVal2 = Rat_solution * Rat_cost;
    conv(OptVal2, optVal2);

    if(optVal != -optVal2)  {
      endanswer -> rest = createListVector(S[i]);
      endanswer = endanswer  -> rest;}
    else {
      endEqu -> rest = createListVector(S[i]);
      endEqu = endEqu -> rest;
     }

  }
  //  equ = equ -> rest;
  return (answer ->rest);
}
