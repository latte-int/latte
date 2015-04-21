/* ReadingFile.cpp -- Reading file and check the input is correct.

   Copyright 2002, 2003 Raymond Hemmecke, Ruriko Yoshida
   Copyright 2006, 2007 Matthias Koeppe

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

#include <string.h>
#include <stdio.h>

#include "config.h"
#include "cone.h"
#include "print.h"
#include "ramon.h"
#include "rational.h"
#include "timing.h"
#include "latte_system.h"
#include "latte_relocatable.h"
#include "ReadingFile.h"
#include "LattException.h"

/* ----------------------------------------------------------------- */

void
CheckRed(char *Filename, char *equ, char* max, char* nonneg, char* interior, char* dil, int dilation)
{
  string fn = Filename;
  CheckRed(fn, equ, max, nonneg, interior, dil, dilation);
  strcpy(Filename, fn.c_str());
}

void
CheckRed(string &Filename, char *equ, char * max, char* nonneg, char* interior, char* dil, int dilation)
{
  int numOfConsts, numOfDims, numOfEqu = 0, flag = 0;
  string tmpString;
  int numOfNonneg = 0, hold = 0;
  vec_ZZ Index, Index2;
  cerr << "Removing redundant inequalities and finding hidden equalities....";
  cerr.flush();
  ifstream IN(Filename.c_str());
  if(!IN){
    cerr << "Input file is missing!!  Please check input file." << endl;
    THROW_LATTE(LattException::ue_FileNameMissing, 0);
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
  ifstream in(Filename.c_str());
  // cerr << "here" << endl;
  in >> numOfConsts >> numOfDims;
  //  cerr << numOfConsts << " " << numOfDims << endl; 
  int equs[numOfConsts];
  //ZZ cost[numOfDims - 1];
  mat_ZZ entries;
  entries.SetDims(numOfConsts, numOfDims);
  // cerr << "here2" << endl;
  if(flag == 2) hold = numOfDims - 1;
  for(int i = 0; i < numOfEqu; i++) conv(equs[i], Index[i]);
  //  cerr << "ther" << endl; 
 for(int i = 0; i < numOfConsts; i++)
    for(int j = 0; j < numOfDims; j++)
      { in >> entries[i][j]; }
 //  cerr << "here3" << endl;
 // printListVector(CheckRedIneq(entries), numOfDims);
 if((equ[0] == 'y') && (flag == 0)){ in >> numOfEqu;
 
 for(int i = 0; i < numOfEqu; i++) in >> equs[i];}
 
 if(dil[0] == 'y')
  for(int i = 0; i < numOfConsts; i++)
    entries[i][0] *= dilation;

 if(interior[0] == 'y'){
   for(int i = 0; i < numOfConsts; i++)
     for(int j = 0; j < numOfDims; j++)
       entries[i][j] *= 10;
  for(int i = 0; i < numOfConsts; i++)
    entries[i][0]--;
 }

 /*  if(max[0] == 'y') 
      for(int i = 0; i < numOfDims - 1; i++) in >> cost[i]; */
 mat_ZZ NONNEG;
 NONNEG.SetDims(numOfDims, numOfDims);
 
  int tmpInt;
  for(int i = 0; i < numOfNonneg; i++){
    conv(tmpInt, Index2[i]); //cerr << Index2[i] << " ";
    NONNEG[tmpInt-1][tmpInt-1] = 1;
  }
  
  ofstream out("Check_red");
  if(!out){
    cerr << "Cannot write Check_red file..." << endl;
    exit(1);
  }
  out << "H-representation" << endl;
  out << "begin" << endl;
  out << numOfConsts + hold << " " << numOfDims << " rational" << endl;
  for(int i = 0; i < numOfConsts; i++){
    for(int j = 0; j < numOfDims; j++){
      out << entries[i][j] << " "; }//cerr << entries[i][j] << " ";}
    out << endl;
  }
  if(nonneg[0] == 'y'){
    for(int i = 0; i < numOfDims - 1; i++){ 
      if(interior[0] == 'y') out << -1 << " ";
      else
	out << 0 << " ";
      for(int j = 0; j < numOfDims - 1; j++) out << NONNEG[i][j] << " ";
      out << endl;
    }
  }
  out << "end" << endl;
  if(equ[0] == 'y'){
    out << "linearity " << numOfEqu << " ";
    for(int i = 0; i < numOfEqu; i++) out << equs[i] << " ";
    out << endl;
  }
 
  int * NewIndex = NULL;
  system_with_error_check(shell_quote(relocated_pathname(REDCHECK_PATH)) + " Check_red > Check_red.out 2>&1");
  int numOfEqu2 = 0, numOfConsts2;
  ifstream in2("Check_red.out");
  if(!in2){
    cerr << "Missing Check_red.out file..." << endl;
    exit(1);
  }
  while(tmpString != "H-representation") {
    in2 >> tmpString;
    if (in2.eof()) {
      cerr << relocated_pathname(REDCHECK_PATH) << " failed to create a non-redundant H-representation; "
	   << "see `Check_red.out' for details." << endl;
      exit(1);
    }
  }
  in2 >> tmpString;
  if(tmpString == "linearity")
    { 
      in2 >> numOfEqu2; 
      NewIndex = new int[numOfEqu2];
      for(int i = 0; i < numOfEqu2; i++)  NewIndex[i] = 0;
      for(int i = 0; i < numOfEqu2; i++) {
	// NewIndex[i] = 0; 
	in2 >> NewIndex[i];
	//  cerr << NewIndex[i];
      }
      
  }
    in2.close();
  if(flag == 2) strcpy(nonneg, "no");
  tmpString = "a";
  if(numOfEqu2 == 0) equ[0] = 'n';
  ifstream in3("Check_red.out");
  if(!in3){
    cerr << "Missing Check_red.out file..." << endl;
    exit(1);
  }
  while (tmpString!="begin") getline(in3,tmpString);

  in3 >> numOfConsts2 >> numOfDims >> tmpString;

  mat_ZZ newEnt;
  newEnt.SetDims(numOfConsts2, numOfDims);
    for(int i = 0; i < numOfConsts2; i++)
      for(int j = 0; j < numOfDims; j++){
	in3 >> newEnt[i][j]; }

    // Don't mess around in the directory of the input files.
    // Put temporary files in the current directory.
    // Then we can safely parallelize tests simply by changing into fresh
    // directories. (Temporary solution.) --mkoeppe
    //strcat(Filename, ".latte");
    Filename = "latte_nonredundant_input";

    ofstream out2(Filename.c_str());

    out2 << numOfConsts2 << " " << numOfDims << endl;
    for(int i = 0; i < numOfConsts2; i++){
      for(int j = 0; j < numOfDims; j++)
	out2 << newEnt[i][j] << " ";
      out2 << endl;
    }
   
    if(numOfEqu2 != 0)
      out2 << numOfEqu2 << " ";
    for(int i = 0; i < numOfEqu2; i++)
      out2 << NewIndex[i] << " ";
    out2 << endl;
    if(numOfEqu2 != 0) equ[0] = 'y';
    /*    if(max[0] == 'y') 
      for(int i = 0; i < numOfDims - 1; i++) out2 << cost[i] << " ";
      out2 << endl; */
    system_with_error_check("rm -f Check_red*");
    cerr << "done." << endl;
    if (NewIndex) delete[] NewIndex;
}

/* ----------------------------------------------------------------- */
void CheckInputFileCDDRep1(const char *InputFile){
  /*
    Author: Ruriko Yoshida
    Date: January 28th, 2003
    Update: January 28th, 2003
    
  */
  
  ifstream in(InputFile);
  ifstream in2(InputFile);
  ifstream in3(InputFile);

  int flag = 0, counter = 0, flag2 = 0, flag3 = 0;
  char tmp[200];
  in3 >> tmp;
  string tmpString;
  //  if(tmp[0] == '*')
    while(tmpString != "begin") {getline(in2, tmpString); 
    if(tmpString[0] == 'l') flag3 = 1;
    //  cerr << tmpString << endl;
    counter++;}
  //cerr << counter << endl;

  if(tmp[0] == '*')
    for(int i = 0; i < counter-2; i++) getline(in, tmpString);
  in >> tmpString;
  if(tmpString == "begin") flag2 = 1;

  if(flag2 == 0){
    if(flag3 == 0){
      if(tmpString  != "H-representation") flag = 1;
      in >> tmpString;
      if(tmpString  != "begin") flag = 1;
      in >> tmpString;
      in >> tmpString;
      in >> tmpString;
      if(tmpString  != "integer") flag = 1;
    }
    //    else if(tmp[0] == '*') ;
    else {
      ifstream in4(InputFile);
      int tmp_int = 0, num; 
      string tmpString2;
      /*      while(tmpString2 != "linearity"){ in4 >>tmpString2; } 
      in4 >> tmp_int; //cerr << tmpString2 << " " << tmp_int << endl;
      for(int i = 0; i < tmp_int; i++) in4 >> num;
      in4 >> tmpString2;//cerr << tmpString2 << endl; */
      if(tmpString2  != "begin") flag = 1;
      in4 >> tmpString2;//cerr << tmpString2 << endl;
      in4 >> tmpString2;//cerr << tmpString2 << endl;
      in4 >> tmpString2;//cerr << tmpString2 << endl;
      if(tmpString2  != "integer") flag = 1;
      while(tmpString2 != "linearity"){ in4 >>tmpString2; } 
      in4 >> tmp_int; //cerr << tmpString2 << " " << tmp_int << endl;
      for(int i = 0; i < tmp_int; i++) in4 >> num;

    }
  }
 
  else if(flag2 == 1){
    
    in >> tmpString;
    in >> tmpString;
    in >> tmpString;
    if(tmpString  != "integer") flag = 1;
    
  }

  while(!in.eof()) in >> tmpString;
  
  //  if(tmpString  != "end") flag = 1;  
  // cerr << flag3 << endl;

  if(flag == 1){
    ofstream out("Error");
    out << "Your input file CDD version is not correct!" << endl;
    cerr << "Your input file CDD version is not correct!" << endl;
    exit (1);
  }
  //  delete [] tmp;
  return ;
}

/* ----------------------------------------------------------------- */
void CheckInputFileCDDRep(const char *InputFile){
  /*
    Author: Ruriko Yoshida
    Date: January 28th, 2003
    Update: January 28th, 2003
    
  */
  
  ifstream in(InputFile);
  int flag = 0;
  //  char* tmpString = new char[200];

  string tmpString;
 
  while(in >> tmpString)
    {
	  if(tmpString == "end") flag ++;
          else if(tmpString == "begin") flag++;
	  //else if(tmpString == "H-representation") flag++;
	  else if(tmpString == "integer") flag++;
	  /*     else if(tmpString == "rational") {
	    ofstream out("Error");
	    out << "Your input file must be integer!" << endl;
	    cerr << "Your input file must be integer!" << endl;
	    exit (1);
	    }*/
          else ;
     
    }

  // cerr << flag << endl;
  if(flag != 3){
    ofstream out("Error");
    out << "Your input file is not correct!" << endl;
    out << "Must be H-representation with integer!" << endl;
    cerr << "Your input file is not correct!" << endl;
    cerr << "Must be H-representation with integer!" << endl;
    exit (1);
  }
  //  delete [] tmpString;
  return ;
}

/* ----------------------------------------------------------------- */

void CheckInputFileCDDRep3(const char *InputFile){
  /*
    Author: Ruriko Yoshida
    Date: January 28th, 2003
    Update: January 28th, 2003
    
  */
  
  ifstream in(InputFile);
  int counter = 0, dim, numOfConst, flag = 0;
  string tmpString;
 
  while (tmpString!="begin") getline(in,tmpString);

  in >> numOfConst >> dim >> tmpString;

  while (tmpString!="end") 
    {in >> tmpString;// cerr << tmpString << endl;
    counter ++;
    }

  if(counter != numOfConst * dim + 1) flag = 1;

  if(flag == 1){
    ofstream out("Error");
    out << "Your input file has wrong number of elements!" << endl;
    cerr << "Your input file has wrong number of elements!" << endl;
    exit (1);
  }

  return ;
}

/* ----------------------------------------------------------------- */
void CheckInputFileCDDRep4(const char *InputFile){
  /*
    Author: Ruriko Yoshida
    Date: January 28th, 2003
    Update: January 28th, 2003
    
  */
  
  ifstream in(InputFile);
  int i, len, flag = 0;
  int dim, constrains;
  string tmp;
  ZZ zztmp;
  char* tmpString = new char[200];
  while(tmp != "begin") getline(in, tmp);
  in >> zztmp;
  conv(constrains, zztmp);
  in >> zztmp;
  conv(dim, zztmp);
  in >> tmp;
  for(int j = 0; j < constrains; j++){
    for(int k = 0; k < dim; k++){
      in >> tmpString;
      len = strlen(tmpString);
      for (i = 0; i < len; i++) 
	{
	  if(tmpString[i] == '0') ;
          else if(tmpString[i] == '-') ;
	  else if(tmpString[i] == '1') ;
	  else if(tmpString[i] == '2') ;
	  else if(tmpString[i] == '3') ;
	  else if(tmpString[i] == '4') ;
	  else if(tmpString[i] == '5') ;
	  else if(tmpString[i] == '6') ;
	  else if(tmpString[i] == '7') ;
	  else if(tmpString[i] == '8') ;
	  else if(tmpString[i] == '9') ;
          else flag = 1;
	}
    }
  }


  if(flag == 1){
    ofstream out("Error");
    out << "Your input file contains non-number!" << endl;
    cerr << "Your input file contains non-number!" << endl;
    exit (1);
  }
  delete [] tmpString;
  return ;
}


/* ----------------------------------------------------------------- */
void CheckInputFile(const char *InputFile){
  /*
    Author: Ruriko Yoshida
    Date: January 28th, 2003
    Update: January 28th, 2003
    
  */
  
  ifstream in(InputFile);
  int i, len, flag = 0;
  char* tmpString = new char[200];
  while(in >> tmpString)
    {
      len = strlen(tmpString);
      if((tmpString[0] != 'l') && (tmpString[0] != 'N') && (tmpString[0] != 'n')){
	for (i = 0; i < len; i++) 
	  {
	    if(tmpString[i] == '0') ;
	    else if(tmpString[i] == '-') ;
	    else if(tmpString[i] == '1') ;
	    else if(tmpString[i] == '2') ;
	    else if(tmpString[i] == '3') ;
	    else if(tmpString[i] == '4') ;
	    else if(tmpString[i] == '5') ;
	    else if(tmpString[i] == '6') ;
	    else if(tmpString[i] == '7') ;
	    else if(tmpString[i] == '8') ;
	    else if(tmpString[i] == '9') ;
	    else flag = 1;
	}
      }
    }

  if(flag == 1){
    ofstream out("Error");
    out << "Your input file contains non-number!" << endl;
    cerr << "Your input file contains non-number!" << endl;
    exit (1);
  }
  delete [] tmpString;
  return ;
}
/* ----------------------------------------------------------------- */
void CheckInputFileVrep(const char *InputFile){
  /*
    Author: Ruriko Yoshida
    Date: January 28th, 2003
    Update: January 28th, 2003
    
  */
  
  ifstream in(InputFile);
  int i, len, flag = 0;
  char* tmpString = new char[200];
  while(in >> tmpString)
    {
      len = strlen(tmpString);
      for (i = 0; i < len; i++) 
	{
	  if(tmpString[i] == '0') ;
          else if(tmpString[i] == '-') ;
	  else if(tmpString[i] == '1') ;
	  else if(tmpString[i] == '2') ;
	  else if(tmpString[i] == '3') ;
	  else if(tmpString[i] == '4') ;
	  else if(tmpString[i] == '5') ;
	  else if(tmpString[i] == '6') ;
	  else if(tmpString[i] == '7') ;
	  else if(tmpString[i] == '8') ;
	  else if(tmpString[i] == '9') ;
	  else if(tmpString[i] == '/') ;
          else flag = 1;
	}
    }

  if(flag == 1){
    ofstream out("Error");
    out << "Your input file contains non-number!" << endl;
    cerr << "Your input file contains non-number!" << endl;
    exit (1);
  }
  delete [] tmpString;
  return ;
}

/* ---------------------------------------------------------------- */

void CheckLength(const char * filename, char * equ)
{
   ifstream in(filename);
   int numOfConstraints = 0, numOfVars = 0, numOfequ = 0;
   in >> numOfConstraints >> numOfVars;
   int Total = numOfConstraints * numOfVars;
   int counts = 0, tmpint2 = 0;
   char tmpint[2000];
   while(in >> tmpint){
     counts++; 
      if(equ[0] == 'y') 
	if(Total == counts-1){ tmpint2 = atoi(tmpint);
	numOfequ = tmpint2 + 1;}
      }
   counts = counts - numOfequ; 
  

   if(Total > counts){
     ofstream out("Error");
     out << "The wrong number of elements in the file.  The number of elements are less than you indicated" << endl;
     cerr <<"The wrong number of elements in the file.  The number of elements are less than you indicated." << endl;
     exit (1);
   }
   /*   else if(Total < counts){
     ofstream out("Error");
     out << "The number of elements are more than you indicated. " << endl;
     cerr << "The number of elements are more than you indicated. " << endl;
     exit (1);
  
     }*/
}

/* ---------------------------------------------------------------------- */

void CheckLength2(const char * filename, char * equ)
{
   ifstream in(filename);
   int numOfConstraints = 0, numOfVars = 0, numOfequ = 0;
   in >> numOfConstraints >> numOfVars;
   int Total = numOfConstraints * numOfVars;
   int counts = 0, tmpint;
   while( in >> tmpint){
      counts++; 
      if(equ[0] == 'y') 
	if(Total == counts-1){numOfequ = tmpint + 1;}
      }
   counts = counts - numOfequ;
   Total += (numOfVars - 1); 

   if(Total > counts){
     ofstream out("Error");
     out << "The wrong number of elements in the file.  The number of elments are less than you indicated" << endl;
     cerr <<"The wrong number of elements in the file.  The number of elments are less than you indicated." << endl;
     exit (1);
   }
   /*   else if(Total < counts){
     ofstream out("Error");
     out << "The number of elements are more than you indicated. " << endl;
     cerr << "The number of elements are more than you indicated. " << endl;
     exit (1);
  
     }*/
}

/* ---------------------------------------------------------------------- */


ZZ FindBigElt(listVector* equation, int numOfVars){
  ZZ bignum;
  while(equation) {
    for(int i = 0; i < numOfVars; i++)
      if(bignum < equation -> first[i])
	bignum = equation -> first[i];
    equation = equation -> rest;
  }

  return bignum;
}

/* ----------------------------------------------------------------- */
void readLatteProblem(const char *fileName, listVector **equations,
		      listVector **inequalities, 
		      char *equationsPresent,
                      int *numOfVars, char *nonneg, char* dual,
		      char* grobner, char * Vrep) {
  int i,j,eq,ind,numOfVectors,numOfEquations;
  vec_ZZ indexEquations;
  listVector *basis, *endBasis, *tmp, *endEquations, *endInequalities;
  vec_ZZ b;
  ZZ bignum;
  /* Reads numOfVars, matrix A, and rhs b. */


  cerr << "Reading problem.\n";

  //setbuf(stdout,0);

  ifstream in(fileName);
  if(!in){
    cerr << "Cannot open input file " << fileName << " in readLatteProblem." << endl;
    exit(1);
  }

  if(grobner[0] == 'y') strcpy(equationsPresent, "yes");
  in >> numOfVectors;
  in >> (*numOfVars); 
  if(Vrep[0] == 'n'){
  if((dual[0] == 'y')&&(equationsPresent[0] == 'n')) *numOfVars = *numOfVars+1;

  int number = 0, oldNumOfVars = 0;

  if(grobner[0] == 'y'){
    ifstream in2(fileName);
    int dim, row;
    in2 >> dim >> row;
    int Matrix[dim][row];
    for(i = 0; i < dim; i++)
      for(j = 0; j < row; j++) in2 >> Matrix[i][j];
    in2 >> number;
  }
  // if(max[0] == 'y') cost.SetLength(*numOfVars - 1);
  if(grobner[0] == 'y'){
    oldNumOfVars = (*numOfVars);
    (*numOfVars) = 2 *(*numOfVars) + 1; 
  }

  mat_ZZ NonNeg;
  mat_ZZ A, B;
  A.SetDims(numOfVectors, *numOfVars);
  
  if((nonneg[0] == 'y')||(grobner[0] == 'y')) {
   NonNeg.SetDims(*numOfVars, *numOfVars);
   for(i = 0; i < *numOfVars - 1; i++)
     NonNeg[i][i+1] = 1;
    }

  if((dual[0] == 'n')||(equationsPresent[0] == 'y')){
  b=createVector(*numOfVars);
  if(grobner[0] == 'y'){
    for (j=1; j<=oldNumOfVars; j++){ in >> b[j]; b[j] = -b[j]; A[0][j] = b[j]; }
    for (j=oldNumOfVars + 1; j<(*numOfVars); j++){b[j] = - b[j - oldNumOfVars]; A[0][j] = b[j];}
    basis = createListVector(b);
    endBasis = basis;
    
    for (i=1; i<number; i++) {
      b=createVector(*numOfVars);
    for (j=1; j<=oldNumOfVars; j++){ in >> b[j]; b[j] = -b[j]; A[i][j] = b[j]; }
    for (j=oldNumOfVars + 1; j<(*numOfVars); j++){b[j] = - b[j - oldNumOfVars]; A[i][j] = b[j];}

      endBasis = updateBasis(createListVector(b), endBasis);
    }

    for (i=number; i<numOfVectors; i++) {
      b=createVector(*numOfVars); 
      b[0] = -1;
    for (j=1; j<=oldNumOfVars; j++){ in >> b[j]; }
    for (j=oldNumOfVars + 1; j<(*numOfVars); j++){b[j] = - b[j - oldNumOfVars];}

      endBasis = updateBasis(createListVector(b), endBasis);
    }

    B = transpose(A); //cerr << B << endl;
    for(i = 0; i < *numOfVars; i++)
    { ZZ tmp_ZZ;
    tmp_ZZ = B[i]*B[i];
    if(bignum < tmp_ZZ)  bignum = tmp_ZZ;
    }
    bignum = (*numOfVars + 1) * (*numOfVars - number) * power(bignum, number/2);

  }
  else{
    for (j=0; j<(*numOfVars); j++){ in >> b[j]; A[0][j] = b[j];}
    basis = createListVector(b);
    endBasis = basis;
    
    for (i=1; i<numOfVectors; i++) {
      b=createVector(*numOfVars);
      
      for (j=0; j<(*numOfVars); j++) {in >> b[j];
      A[i][j] = b[j];
      }
      endBasis = updateBasis(createListVector(b), endBasis);
    }
    B = transpose(A);
    for(i = 0; i < *numOfVars; i++)
    { ZZ tmp_ZZ;
    tmp_ZZ = B[i]*B[i];
    if(bignum < tmp_ZZ)  bignum = tmp_ZZ;
    }
    bignum = (*numOfVars + 1) * (*numOfVars - numOfVectors) * power(bignum, numOfVectors/2);
  }
 
  if(grobner[0] == 'y'){
  mat_ZZ UB;
  UB.SetDims(*numOfVars, *numOfVars);
  for(i = 0; i < *numOfVars; i++){
    UB[i][0] = bignum;
   
  }

  for(i = 0; i < (*numOfVars - 1); i++){
    UB[i][i + 1] = -1;
  }
   for(i = 0; i < *numOfVars - 1; i++){
     b=createVector(*numOfVars);
     for(j = 0; j < *numOfVars; j++) b[j] = UB[i][j];
     endBasis = updateBasis(createListVector(b), endBasis);
   }
   for(i = 0; i < *numOfVars - 1; i++){
     b=createVector(*numOfVars);
     for(j = 0; j < *numOfVars; j++) b[j] = NonNeg[i][j];
     endBasis = updateBasis(createListVector(b), endBasis);
   }
  }
 

  if(nonneg[0] == 'y'){
   for(i = 0; i < *numOfVars - 1; i++){
     b=createVector(*numOfVars);
     for(j = 0; j < *numOfVars; j++) b[j] = NonNeg[i][j];
     endBasis = updateBasis(createListVector(b), endBasis);
   }
  }
  }
  /*
  in >> numOfVectors;
  in >> (*numOfVars);
  if((dual[0] == 'y')&&(equationsPresent[0] == 'n')) *numOfVars = *numOfVars+1;
  mat_ZZ NonNeg;
  if(nonneg[0] == 'y') {
   NonNeg.SetDims(*numOfVars, *numOfVars);
   for(i = 0; i < *numOfVars - 1; i++)
     NonNeg[i][i+1] = 1;
    }
  if((dual[0] == 'n')||(equationsPresent[0] == 'y')){
  b=createVector(*numOfVars);
  for (j=0; j<(*numOfVars); j++) in >> b[j];
  basis = createListVector(b);
  endBasis = basis;

  for (i=1; i<numOfVectors; i++) {
    b=createVector(*numOfVars);
    for (j=0; j<(*numOfVars); j++) in >> b[j];
    endBasis = updateBasis(createListVector(b), endBasis);
  }

  if(nonneg[0] == 'y'){
   for(i = 0; i < *numOfVars - 1; i++){
     b=createVector(*numOfVars);
     for(j = 0; j < *numOfVars; j++) b[j] = NonNeg[i][j];
     endBasis = updateBasis(createListVector(b), endBasis);
   }
  }
  }
  */
  if((dual[0] == 'y')&&(equationsPresent[0] == 'n')){
 
  b=createVector(*numOfVars);
  ZZ hold;
  in >> hold;// cerr << hold << endl;
  for (j=1; j<(*numOfVars)-1; j++) in >> b[j];
  b[*numOfVars-1] = hold;// cerr << b << endl;
  basis = createListVector(b);
  endBasis = basis;

  for (i=1; i<numOfVectors; i++) {
    b=createVector(*numOfVars);
    in >> hold;//cerr << hold << endl;
    for (j=1; j<(*numOfVars)-1; j++) in >> b[j];
  b[*numOfVars-1] = hold;
    endBasis = updateBasis(createListVector(b), endBasis);
  }

  }
  if (equationsPresent[0]=='n') {
    (*inequalities)=basis;
    (*equations)=0;
  } else {

    /* Read indices of equations and split basis into list of
       equations and inequalities. */

    in >> numOfEquations;
    indexEquations=createVector(numOfEquations);

    for (i=0; i<numOfEquations; i++) in >> indexEquations[i];
#if 0
    cerr << "\nEquation indices: ";
    printVectorToFile(cerr,indexEquations,numOfEquations);
#endif

    (*equations)=createListVector(createVector(*numOfVars));
    (*inequalities)=createListVector(createVector(*numOfVars));
    endEquations=(*equations);
    endInequalities=(*inequalities);

    eq=0;
    ind=1;

    tmp=basis;
    while (tmp) {
      if (ind==indexEquations[eq]) {
	endEquations->rest=createListVector(tmp->first);
	endEquations=endEquations->rest;
	eq++;
	{
	  listVector *t = tmp->rest;
	  delete tmp;
	  tmp=t;
	}
	if (eq==numOfEquations) {
	  endInequalities->rest=tmp;
	  tmp=0;
	}	
      } else {
	endInequalities->rest=createListVector(tmp->first);
	endInequalities=endInequalities->rest;
	{
	  listVector *t = tmp->rest;
	  delete tmp;
	  tmp=t;
	}
      }
      ind++;
    }
    {
      listVector *t = (*equations)->rest;
      delete *equations;
      *equations = t;
    }
    {
      listVector *t = (*inequalities)->rest;
      delete *inequalities;
      *inequalities = t;
    }
  }
  //if(max[0] == 'y') for(i = 0; i < (*numOfVars-1); i++) in >> cost[i];
  if(Vrep[0] == 'n'){
#if 0
  cerr << endl;
  cerr << "Ax <= b, given as (b|-A):\n";
  cerr << "=========================\n";
  printListVectorToFile(cerr, *inequalities,*numOfVars);
  
  cerr << endl;

  cerr << "Ax = b, given as (b|-A):\n";
  cerr << "========================\n";
  printListVectorToFile(cerr, *equations,*numOfVars);

  cerr << endl;
#endif
  }

  else{
  cerr << endl;
  cerr << "The vertex Set, given by (1| V):\n";
  cerr << "=========================\n";
  printListVectorToFile(cerr, *inequalities,*numOfVars);
  
  cerr << endl;
  }
  }
  return;
}

/* ----------------------------------------------------------------- */
int CDDstylereadLatteProblem(const char *fileName, listVector **equations,
		      listVector **inequalities, 
		      char *equationsPresent,
                      int *numOfVars, char *nonneg, char* dual,
                      char* taylor, int & degree, 
                      char* rat, int & cone_output, 
                      char* Memory_Save, char* uni, char* inthull,
		      char* grobner) {
  int i,j,eq,ind, length = 0, f = 0, numOfVectors,numOfEquations;
  vec_ZZ indexEquations;
  listVector *basis, *endBasis, *tmp, *endEquations, *endInequalities;
  vec_ZZ b;
  string tmpString;
  
  if(grobner[0] == 'y') strcpy(equationsPresent, "yes");
  if(rat[0] == 'y') strcpy(rat, "yes");
  /* Reads numOfVars, matrix A, and rhs b. */

  ifstream in2(fileName);
  if(!in2){
    cerr << "Cannot open input file " << fileName << " in CDDstylereadLatteProblem." << endl;
    exit(1);
  }

  while(in2 >> tmpString){ //cerr << tmpString << endl;

    if(tmpString == "dual") {strcpy(dual, "yes");}
    else if(tmpString == "nonneg") strcpy(nonneg, "yes");
    else if(tmpString == "taylor") 
        {strcpy(taylor, "yes");
	in2 >> degree; //cerr << degree << endl;
	}
    else if(tmpString == "rational") 
        {strcpy(rat, "yes");
	in2 >> cone_output;
	}
    else if(tmpString == "unimodular") strcpy(uni, "yes");
    else if(tmpString == "memsave") strcpy(Memory_Save, "yes");
    else if(tmpString == "inthull") strcpy(inthull, "yes");
    else if(tmpString == "grobner") strcpy(grobner, "yes");
    else if(tmpString == "linearity") {strcpy(equationsPresent, "yes");
    in2 >> length;
    indexEquations.SetLength(length);
    for(i = 0; i < length; i++) in2 >> indexEquations[i];
    }
  }

  if(grobner[0] == 'y') return 0;
  else{
    // int length = 0;
  cerr << "Reading problem.\n";
  //  cerr << degree << endl;
  //setbuf(stdout,0);
  ZZ bignum;
  ifstream in(fileName);
  if(!in){
    cerr << "Cannot open input file " << fileName << " in CDDstylereadLatteProblem." << endl;
    exit(1);
  }
    while (tmpString!="begin"){
      getline(in,tmpString);
      if(tmpString[0] == 'l') f = 1;}

    if(f== 1){
      ifstream in2(fileName);
      strcpy(equationsPresent,"yes");
      while (tmpString!="linearity") in2 >> tmpString;
      //  int length;
      in2 >> length;
      //     cerr << length << endl;
      indexEquations.SetLength(length);
      for(i = 0; i < length; i++) in2 >> indexEquations[i];
    }
  in >> numOfVectors;
  in >> (*numOfVars) >> tmpString;
  if((dual[0] == 'y')&&(equationsPresent[0] == 'n')) *numOfVars = *numOfVars+1;
  int oldNumOfVars = 0;
  if(grobner[0] == 'y'){
    oldNumOfVars = (*numOfVars);
    (*numOfVars) = 2 *(*numOfVars) + 1; 
  }

  mat_ZZ NonNeg;
  mat_ZZ A, B;
  A.SetDims(numOfVectors, *numOfVars);

  if((nonneg[0] == 'y')||(grobner[0] == 'y')) {
   NonNeg.SetDims(*numOfVars, *numOfVars);
   for(i = 0; i < *numOfVars - 1; i++)
     NonNeg[i][i+1] = 1;
    }

  if((dual[0] == 'n')||(equationsPresent[0] == 'y')){
  b=createVector(*numOfVars);
  if(grobner[0] == 'y'){
    for (j=1; j<=oldNumOfVars; j++){ in >> b[j]; b[j] = -b[j]; A[0][j] = b[j]; }
    for (j=oldNumOfVars + 1; j<(*numOfVars); j++){b[j] = - b[j - oldNumOfVars]; A[0][j] = b[j];}
    basis = createListVector(b);
    endBasis = basis;
    
    for (i=1; i<length; i++) {
      b=createVector(*numOfVars);
    for (j=1; j<=oldNumOfVars; j++){ in >> b[j]; b[j] = -b[j]; A[i][j] = b[j]; }
    for (j=oldNumOfVars + 1; j<(*numOfVars); j++){b[j] = - b[j - oldNumOfVars]; A[i][j] = b[j];}

      endBasis = updateBasis(createListVector(b), endBasis);
    }

    for (i=length; i<numOfVectors; i++) {
      b=createVector(*numOfVars); 
      b[0] = -1;
    for (j=1; j<=oldNumOfVars; j++){ in >> b[j];  }
    for (j=oldNumOfVars + 1; j<(*numOfVars); j++){b[j] = - b[j - oldNumOfVars];}

      endBasis = updateBasis(createListVector(b), endBasis);
    }

    B = transpose(A); cerr << B << endl;
    for(i = 0; i < *numOfVars; i++)
    { ZZ tmp_ZZ;
    tmp_ZZ = B[i]*B[i];  //cerr << B[i] << endl;
    if(bignum < tmp_ZZ)  bignum = tmp_ZZ;
    }
    bignum = (*numOfVars + 1) * (*numOfVars - length) * power(bignum, length/2);

  }

  else{
    for (j=0; j<(*numOfVars); j++){ in >> b[j]; A[0][j] = b[j];}
    basis = createListVector(b);
    endBasis = basis;
    
    for (i=1; i<numOfVectors; i++) {
      b=createVector(*numOfVars);
      for (j=0; j<(*numOfVars); j++) {in >> b[j];
    A[i][j] = b[j];
      }
      endBasis = updateBasis(createListVector(b), endBasis);
    }
    B = transpose(A);
    for(i = 0; i < *numOfVars; i++)
    { ZZ tmp_ZZ;
    tmp_ZZ = B[i]*B[i];
    if(bignum < tmp_ZZ)  bignum = tmp_ZZ; //cerr << bignum << endl;
    }
    bignum = (*numOfVars + 1) * (*numOfVars - length) * power(bignum, (length)/2);
  }
    if(grobner[0] == 'y'){
      mat_ZZ UB;
      UB.SetDims(*numOfVars, *numOfVars);
      for(i = 0; i < *numOfVars; i++){
	UB[i][0] = bignum;
      }

  for(i = 0; i < (*numOfVars - 1); i++){
    UB[i][i + 1] = -1;
  }
   for(i = 0; i < *numOfVars - 1; i++){
     b=createVector(*numOfVars);
     for(j = 0; j < *numOfVars; j++) b[j] = UB[i][j];
     endBasis = updateBasis(createListVector(b), endBasis);
   }
   for(i = 0; i < *numOfVars - 1; i++){
     b=createVector(*numOfVars);
     for(j = 0; j < *numOfVars; j++) b[j] = NonNeg[i][j];
     endBasis = updateBasis(createListVector(b), endBasis);
   }
  }
 

  if(nonneg[0] == 'y'){
   for(i = 0; i < *numOfVars - 1; i++){
     b=createVector(*numOfVars);
     for(j = 0; j < *numOfVars; j++) b[j] = NonNeg[i][j];
     endBasis = updateBasis(createListVector(b), endBasis);
   }
  }
  }

  if((dual[0] == 'y')&&(equationsPresent[0] == 'n')){
 
  b=createVector(*numOfVars);
  ZZ hold;
  in >> hold;// cerr << hold << endl;
  for (j=1; j<(*numOfVars)-1; j++) in >> b[j];
  b[*numOfVars-1] = hold;// cerr << b << endl;
  basis = createListVector(b);
  endBasis = basis;

  for (i=1; i<numOfVectors; i++) {
    b=createVector(*numOfVars);
    in >> hold;//cerr << hold << endl;
    for (j=1; j<(*numOfVars)-1; j++) in >> b[j];
  b[*numOfVars-1] = hold;
    endBasis = updateBasis(createListVector(b), endBasis);
  }

  }
  if (equationsPresent[0]=='n') {
    (*inequalities)=basis;
    (*equations)=0;
  } else {

    /* Read indices of equations and split basis into list of
       equations and inequalities. */
    /*   if(f == 0)
      in >> numOfEquations;
    indexEquations=createVector(numOfEquations);
    if(f == 0)
    for (i=0; i<numOfEquations; i++) in >> indexEquations[i]; */
    numOfEquations = 0;
    cerr << "\nEquation indices: ";
    printVectorToFile(cerr,indexEquations,numOfEquations);

    (*equations)=createListVector(createVector(*numOfVars));
    (*inequalities)=createListVector(createVector(*numOfVars));
    endEquations=(*equations);
    endInequalities=(*inequalities);

    eq=0;
    ind=1;

    tmp=basis;
    while (tmp) {
      if (ind==indexEquations[eq]) {
	endEquations->rest=createListVector(tmp->first);
	endEquations=endEquations->rest;
	eq++;
	tmp=tmp->rest;
	if (eq==numOfEquations) {
	  endInequalities->rest=tmp;
	  tmp=0;
	}	
      } else {
	endInequalities->rest=createListVector(tmp->first);
	endInequalities=endInequalities->rest;
	tmp=tmp->rest;
      }
      ind++;
    }
    (*equations)=(*equations)->rest;
    (*inequalities)=(*inequalities)->rest;
  }

  cerr << endl;
  cerr << "Ax <= b, given as (b|-A):\n";
  cerr << "=========================\n";
  printListVectorToFile(cerr, *inequalities,*numOfVars);
  
  cerr << endl;

  cerr << "Ax = b, given as (b|-A):\n";
  cerr << "========================\n";
  printListVectorToFile(cerr, *equations,*numOfVars);

  cerr << endl;

  return 0;
  }
}


/* ----------------------------------------------------------------- */
