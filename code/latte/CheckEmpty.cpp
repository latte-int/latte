/*                                                                    *
 * Author: Ruriko Yoshida                                             *
 * Date: Febrary 10th, 2004                                           *
 * Update: Febrary 10th, 2004                                         *
 * This is for checking whether the input polytope is empty or not.   *
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

#include "preprocess.h"
#include "print.h"
#include "ramon.h"
#include "rational.h"
#include "timing.h"
#include "flags.h"

#include "IntegralHull.h"
#include "ReadingFile.h"
#include "binarySearchIP.h"

void CheckEmpty(char * Filename){
  int numOfConsts, numOfDims, numOfEqu = 0, flag = 0;
  string tmpString;
  char equ[127], nonneg[127];
  equ[0] = 'n';
  nonneg[0] = 'n';
  int numOfNonneg = 0, hold = 0;
  vector Index, Index2;
  cout << "Checking whether the input polytope is empty or not..." << endl;
  ifstream IN(Filename);
  if(!IN){
    cerr << "Input file is missing!!  Please check input file." << endl;
    exit(0);
  }
  while(IN >> tmpString){
    if(tmpString == "linearity"){
      strcpy(equ, "yes");
      IN >> numOfEqu; 
      Index.SetLength(numOfEqu);
      for(int i = 0; i < numOfEqu; i++) IN >> Index[i];
      flag = 1;
    }
    else if((tmpString == "NONNEGATIVE") || (tmpString == "nonnegative")){
      strcpy(nonneg, "yes");
      flag = 2;
      IN >> numOfNonneg;
      Index2.SetLength(numOfNonneg);
      for(int i = 0; i < numOfNonneg; i++) IN >> Index2[i];
    }
  }
  ifstream in(Filename);
  // cout << "here" << endl;
  in >> numOfConsts >> numOfDims;
  //  cout << numOfConsts << " " << numOfDims << endl; 
  int equs[numOfConsts];
  //ZZ cost[numOfDims - 1];
  mat_ZZ entries;
  entries.SetDims(numOfConsts, numOfDims);
  // cout << "here2" << endl;
  if(flag == 2) hold = numOfDims - 1;
  for(int i = 0; i < numOfEqu; i++) conv(equs[i], Index[i]);
  //  cout << "ther" << endl; 
 for(int i = 0; i < numOfConsts; i++)
    for(int j = 0; j < numOfDims; j++)
      { in >> entries[i][j]; }
 //  cout << "here3" << endl;
 // printListVector(CheckRedIneq(entries), numOfDims);
 if((equ[0] == 'y') && (flag == 0)){ in >> numOfEqu;
 
 for(int i = 0; i < numOfEqu; i++) in >> equs[i];}
 // for(int i = 0; i < numOfEqu; i++) cout << equs[i] << endl;
 /*  if(max[0] == 'y') 
      for(int i = 0; i < numOfDims - 1; i++) in >> cost[i]; */
 mat_ZZ NONNEG;
 NONNEG.SetDims(numOfDims, numOfDims);
 
  int tmpInt;
  for(int i = 0; i < numOfNonneg; i++){
    conv(tmpInt, Index2[i]); //cout << Index2[i] << " ";
    NONNEG[tmpInt-1][tmpInt-1] = 1;
  }
  
  ofstream out("Check_emp.lp");
  if(!out){
    cerr << "Cannot write Check_red file..." << endl;
    exit(0);
  }

  out << "H-representation" << endl;
  out << "begin" << endl;
  out << numOfConsts + hold + numOfEqu << " " << numOfDims << " rational" << endl;
  for(int i = 0; i < numOfConsts; i++){
    for(int j = 0; j < numOfDims; j++){
      out << entries[i][j] << " ";}//  cout << entries[i][j] << " ";}
    out << endl;
  }//cout << "here" << endl;
  if(equ[0] == 'y'){
    for(int i = 0; i < numOfEqu; i++){out << -entries[equs[i]-1][0] << " ";
    for(int j = 1; j < numOfDims; j++){
      out << -entries[equs[i]-1][j] << " "; }//cout << entries[i][j] << " ";}
    out << endl;
    }
  }
  if(nonneg[0] == 'y'){
    for(int i = 0; i < numOfDims - 1; i++){ out << 0 << " ";
      for(int j = 0; j < numOfDims - 1; j++) out << NONNEG[i][j] << " ";
      out << endl;
    }
  }
  out << "end" << endl;
  /* if(equ[0] == 'y'){
    out << "linearity " << numOfEqu << " ";
    for(int i = 0; i < numOfEqu; i++) out << equs[i] << " ";
    out << endl;
    }*/
  out<< "maximize" << endl; out << 0 << " " ;
  for(int j = 1; j < numOfDims; j++){
    out << entries[0][j] << " "; }
  out << endl;

  system("./cdd Check_emp.lp > Check_emp.out");

  int FLAG = 0;

  ifstream IN3("Check_emp.lps");
  if(!IN3){
    cerr << "Check_emp.lps is missing!!  Please check input file." << endl;
    exit(0);
  }
  while(IN3 >> tmpString){
    //  cout << tmpString << endl;
    if((tmpString == "primal_solution") || ( tmpString == "primal_direction")){
      FLAG = 1;
    }
  }

  //  system("rm Check_emp.*");

  if(FLAG == 0) {
    cerr << "Empty polytope or unbounded polytope!"<< endl;
    ofstream NOL("numOfLatticePoints");
    NOL << 0 << endl;
    exit (0);
  }
  cout << "done." << endl;
}
