/* binarySearch.cpp -- Binary Search Integer Programming method

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

#include <iostream>
#include <string.h>
#include <stdio.h>

#include "cone.h"
#include "print.h"
#include "ramon.h"
#include "rational.h"
#include "timing.h"
#include "flags.h"
#include "binarySearchIP.h"
#include "latte_system.h"
#include "latte_relocatable.h"

#include "../config.h"

using namespace std;

ZZ  NOT_FOUND;

/* ----------------------------------------------------------------- */
void createCddIneLPFile2(listVector* matrix, listVector* ineq, int numOfVars, vec_ZZ & cost) {
  int i;
  listVector* tmp, *tmp2;
  // cerr << cost << " " << numOfVars << endl;
  ofstream OUT;

  OUT.open("LP.ine");
  OUT << "H-representation" << endl;
  /*if(lengthListVector(matrix) != 0) {
     OUT << "linearity " << lengthListVector(matrix) << " ";
     for(i = 0; i < lengthListVector(ineq); i++) OUT << i+1 << " ";
     OUT << endl;
     }*/
  OUT << "begin " << endl;
  OUT << 2*lengthListVector(matrix) + lengthListVector(ineq) << " " << numOfVars << " integer" << endl;
  tmp=matrix;  tmp2 = ineq;
  while (tmp) {
    for (i=0; i<(numOfVars); i++) OUT << (tmp->first)[i] << " ";
    OUT << endl;
    for (i=0; i<(numOfVars); i++) OUT << -(tmp->first)[i] << " ";
    OUT << endl;
    tmp=tmp->rest;
  }
    while (tmp2) {
    for (i=0; i<(numOfVars); i++) OUT << (tmp2->first)[i] << " ";
    OUT << endl;
    tmp2=tmp2->rest;
  }
  OUT << "end" << endl;
  OUT << "maximize" << endl;
  OUT << 0 << " ";
  for(i = 0; i < numOfVars - 1; i++)
    OUT << cost[i] << " ";
  OUT << endl;
  OUT.close();

  return ;
}

/* ----------------------------------------------------------------- */

rationalVector* ReadLpsFile2(int numOfVars)
{
  ifstream in("LP.lps");
  string tmpString;
  ZZ x, y;
  cerr << "Reading .lps file...";
  rationalVector* OptVector;
  OptVector = createRationalVector(numOfVars);
  if(!in){
    cerr << "Cannot open input file in ReadLpsFile." << endl;
    exit(1);
  }

  while (tmpString!="begin") getline(in,tmpString);
  in >> tmpString;
  for(int i = 0; i < numOfVars; i++)
    {
      in >> tmpString >> tmpString;

      x=0;
      y=0;
      ReadCDD(in,x,y);
      OptVector->set_entry(i, x, y);

    }
  cerr <<"done." << endl;
  return OptVector;
}

/* ----------------------------------------------------------------- */
rationalVector* LP2(listVector* matrix, listVector* ineq, vec_ZZ& cost, int numOfVars) {

  rationalVector* Opt_vector;
   createCddIneLPFile2(matrix, ineq, numOfVars+1, cost);
   cerr << "Computing LP...";
   system_with_error_check(shell_quote(relocated_pathname(CDD_PATH)) + " LP.ine > LP.out");
  cerr << "done.\n\n";
  Opt_vector = ReadLpsFile2(numOfVars);
  //  cerr << Opt_vector->numerators() << " " << Opt_vector -> denominator << endl;
  system_with_error_check("rm -f LP.*");

  return(Opt_vector);
}
/* ----------------------------------------------------------------- */

void createLatteFileEqu(listVector* matrix, listVector* ineq, int numOfVars, ZZ rhs, vec_ZZ lhs) {
  int i;
  listVector* tmp, *tmp2;
  //cerr << rhs << endl;
  ofstream OUT;
  OUT.open("latte_BS");
  //  OUT << "H-representation" << endl;
  // OUT << "begin " << endl;
  OUT << lengthListVector(matrix)+1 + lengthListVector(ineq)<< " " << numOfVars << endl; //" integer" << endl;
  tmp=matrix, tmp2 = ineq;
  while (tmp) {
    for (i=0; i<(numOfVars); i++) OUT << (tmp->first)[i] << " ";
    OUT << endl;
    tmp=tmp->rest; //cerr << "here" << endl;
  }
  OUT << -rhs << " ";
  for (i=0; i<(numOfVars) - 1; i++) OUT << lhs[i] << " ";
  OUT << endl;
    while (tmp2) {
    for (i=0; i<(numOfVars); i++) OUT << (tmp2->first)[i] << " ";
    OUT << endl;
    tmp2=tmp2->rest; //cerr << "here" << endl;
  }
  // OUT << "end" << endl;
  // OUT << "linearity" << " "<< 1 << " " << 1 << endl;
    OUT << "linearity " << 1 + lengthListVector(matrix) << " " << 1 << " ";
    for(i = 0; i < lengthListVector(matrix); i++) OUT << 2 + i << " ";
    OUT << endl;
  OUT.close();

  return ;
}
/* ----------------------------------------------------------------- */

ZZ OptimalCheckEqu(listVector* matrix, listVector* ineq, int numOfVars, ZZ rhs, vec_ZZ lhs)
{
  ZZ NumOfLatticePoints;

  //system_with_error_check("rm numOfLatticePoints");
  createLatteFileEqu(matrix, ineq, numOfVars + 1, rhs, lhs);
  system_with_error_check(shell_quote(relocated_pathname(string(INSTALLPREFIX) + "/bin/count")) + " latte_BS > latte_BS.out");

  ifstream in("numOfLatticePoints");
  in >> NumOfLatticePoints;
  system_with_error_check("rm -f latte_BS*");

  return NumOfLatticePoints;

}

/* ----------------------------------------------------------------- */

void createLatteFile(listVector* matrix, listVector* ineq, int numOfVars, ZZ rhs, vec_ZZ lhs) {
  int i;
  listVector* tmp, *tmp2;
  //cerr << rhs << endl;
  ofstream OUT;
  OUT.open("latte_BS");
  //  OUT << "H-representation" << endl;
  // OUT << "begin " << endl;
  OUT << lengthListVector(matrix)+1 + lengthListVector(ineq)<< " " << numOfVars << endl; //" integer" << endl;
  tmp=matrix, tmp2 = ineq;
  while (tmp) {
    for (i=0; i<(numOfVars); i++) OUT << (tmp->first)[i] << " ";
    OUT << endl;
    tmp=tmp->rest; //cerr << "here" << endl;
  }
  OUT << -rhs << " ";
  for (i=0; i<(numOfVars) - 1; i++) OUT << lhs[i] << " ";
  OUT << endl;
    while (tmp2) {
    for (i=0; i<(numOfVars); i++) OUT << (tmp2->first)[i] << " ";
    OUT << endl;
    tmp2=tmp2->rest; //cerr << "here" << endl;
  }

  if(lengthListVector(matrix) != 0){
    OUT << "linearity " << lengthListVector(matrix) << " ";
    for(i = 0; i < lengthListVector(matrix); i++) OUT << 1 + i << " ";
    OUT << endl;     }
  OUT.close();

  return ;
}
/* ----------------------------------------------------------------- */

ZZ OptimalCheck(listVector* matrix, listVector* ineq, int numOfVars, ZZ rhs, vec_ZZ lhs, ZZ & TotalNumOfUniCones)
{
  ZZ NumOfLatticePoints;

  //system_with_error_check("rm -f numOfLatticePoints");
  createLatteFile(matrix, ineq, numOfVars + 1, rhs, lhs);
  if(lengthListVector(matrix) != 0)
    system_with_error_check("time " + shell_quote(relocated_pathname(string(INSTALLPREFIX) + "/bin/count")) + " latte_BS > latte_BS.out");
  else
    system_with_error_check("time " + shell_quote(relocated_pathname(string(INSTALLPREFIX) + "/bin/count")) + " latte_BS > latte_BS.out");

  ifstream in("numOfLatticePoints");
  in >> NumOfLatticePoints;

  ifstream in2("numOfUnimodularCones");
  ZZ numOfUniCones;
  in2 >> numOfUniCones;
  TotalNumOfUniCones += numOfUniCones;
  cerr << "Number of Unimodular cones: " << numOfUniCones << endl;

  system_with_error_check("rm -f latte_BS*");

  return NumOfLatticePoints;

}

/* ----------------------------------------------------------------- */

ZZ binarySearch(listVector* matrix, listVector* ineq, vec_ZZ cost, int numOfVars, char * min)
{
      int i = 0, counter = 0;
      ZZ low, high, mid, Opt, TotalNumOfUniCones;
      rationalVector *Low_opt, *High_opt;
      vec_RR Low_solution, High_solution, tmp_den, tmp_num, Rat_cost, Rat_cost2;
      RR Rat_low, Rat_high;
      vec_ZZ low_cost;
      // cerr << cost << endl;
      low_cost = -cost;
      //      cerr << cost << endl;
      //  printListVector(matrix, numOfVars); exit(0);
      ofstream tmpFile("numOfLatticePoints");
      tmpFile << "Junk." << endl;
      system_with_error_check("rm -f numOfLatticePoints");
      Low_opt =  LP2(matrix, ineq, low_cost, numOfVars);
      High_opt =  LP2(matrix, ineq, cost, numOfVars);

      //cerr << Low_opt ->numerators() << " " << High_opt->numerators() << endl;
      //cerr << Low_opt ->denominators() << " " << High_opt->denominators()<< endl;
      cerr << "A optimal solution for LP relax.: ";
      printRationalVector(High_opt, numOfVars);
      Low_solution.SetLength(numOfVars);
      High_solution.SetLength(numOfVars);
      tmp_den.SetLength(numOfVars);
      tmp_num.SetLength(numOfVars);
      Rat_cost.SetLength(numOfVars);
      Rat_cost2.SetLength(numOfVars);

      for (i = 0; i < numOfVars; i++){
      	 conv(tmp_num[i], High_opt->numerators()[i]);
	 conv(tmp_den[i], High_opt->denominators()[i]);
       	 High_solution[i] = tmp_num[i]/tmp_den[i];
	 conv(Rat_cost[i], cost[i]);
     		}

      for (i = 0; i < numOfVars; i++){
      	 conv(tmp_num[i], Low_opt->numerators()[i]);
          conv(tmp_den[i], Low_opt->denominators()[i]);
       	 Low_solution[i] = tmp_num[i]/tmp_den[i];
     		}

      Rat_high = Rat_cost * High_solution;
      Rat_low = (Rat_cost) * Low_solution;
      cerr << "The optimal value for LP relax.: " << Rat_high << endl << endl;
      high = CeilToZZ(Rat_high);
      Opt = OptimalCheckEqu(matrix, ineq, numOfVars, high, cost);

      high += 1;
      low = CeilToZZ(Rat_low);

      if(high < low) {
        ZZ tmp_ZZ;
        tmp_ZZ = high;
        high = low;
        low = tmp_ZZ;
      }

      if(IsZero(Opt) == 0){
        cerr << "The optimal value: " << high - 1 << endl << endl;
	return Opt;
	}
       else{
	while(IsOne(high - low) == 0)
	  {
	    mid = (low + high) / 2;
	    //cerr << mid << endl;
	    Opt = OptimalCheck(matrix, ineq, numOfVars, mid, cost, TotalNumOfUniCones);
	    if(IsZero(Opt) == 0)
	      low = mid;

	    else if(IsZero(Opt) == 1)
	      high = mid;
	    // else return mid; //found
	    counter++;
	    if((counter % 10) == 0){
	      cerr << "Iterations: " << counter << endl;
	    }
	  }
      }
      Opt = OptimalCheck(matrix, ineq, numOfVars, low, cost, TotalNumOfUniCones);
      cerr << endl << "Total of Iterations: " << counter << endl;
      cerr << "The total number of unimodular cones: " << TotalNumOfUniCones << endl;
      if(min[0] == 'y')
	cerr << "The optimal value: " << -low << endl << endl;
      else
	cerr << "The optimal value: " << low << endl << endl;
      return Opt;

}







