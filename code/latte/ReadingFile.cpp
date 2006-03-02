/* ----------------------------------------------------------------- */
/*                                                                   */
/* LattE (Lattice Point Enumeration)                                 */
/*                                                                   */
/* Reading file and check the input is correct.                      */
/*                                                                   */
/* Author     : Raymond Hemmecke, Ruriko Yoshida                     */
/*                                                                   */
/* Created    : 07-JUN-02                                            */
/* Last Update: 03-Mar-03                                            */
/*                                                                   */
/* ----------------------------------------------------------------- */

#include <string.h>
#include <stdio.h>

#include "config.h"
#include "myheader.h"
#include "barvinok/dec.h"
#include "barvinok/barvinok.h"
#include "barvinok/ConeDecom.h"
#include "barvinok/Triangulation.h"
#include "vertices/cdd.h"
#include "genFunction/maple.h"
#include "genFunction/piped.h"
#include "cone.h"
#include "dual.h"
#include "RudyResNTL.h"
//  #include "jesus.h"
#include "preprocess.h"
#include "print.h"
#include "ramon.h"
#include "rational.h"
#include "timing.h"
#include "flags.h"
//#include "testing.h"
#include "IntegralHull.h"

void IntVector(listVector* basis, int numOfVars) {

  while(basis) {
    basis->first = -basis->first;
    basis = basis->rest;
  }
/*  printf("\n"); */
  return ;
}
/* ----------------------------------------------------------------- */
void IntCone(listCone* cones, int numOfVars) {
  IntVector(cones->rays,numOfVars);

  return ;
}

/* ----------------------------------------------------------------- */
listCone* IntCone2(listCone* cones, int numOfVars) {
 listCone* cone;
 cone = cones;
  while(cone) {
    IntCone(cone,numOfVars);
    cone = cone->rest;
  }
  //cout << endl;
  //cones = cone;
  return cones;
}

/* ----------------------------------------------------------------- */

void CheckRed(char* Filename, char *equ, char * max, char* nonneg, char* interior, char* dil, int dilation){
  int numOfConsts, numOfDims, numOfEqu = 0, flag = 0;
  string tmpString;
  int numOfNonneg = 0, hold = 0;
  vec_ZZ Index, Index2;
  cout << "Removing redundant inequalities and finding hidden equalities...." << endl;
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
    conv(tmpInt, Index2[i]); //cout << Index2[i] << " ";
    NONNEG[tmpInt-1][tmpInt-1] = 1;
  }
  
  ofstream out("Check_red");
  if(!out){
    cerr << "Cannot write Check_red file..." << endl;
    exit(0);
  }
  out << "H-representation" << endl;
  out << "begin" << endl;
  out << numOfConsts + hold << " " << numOfDims << " rational" << endl;
  for(int i = 0; i < numOfConsts; i++){
    for(int j = 0; j < numOfDims; j++){
      out << entries[i][j] << " "; }//cout << entries[i][j] << " ";}
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
 
  int * NewIndex;
  system(REDCHECK_PATH " Check_red > Check_red.out 2>&1");
  int numOfEqu2 = 0, numOfConsts2;
  ifstream in2("Check_red.out");
  if(!in2){
    cerr << "Missing Check_red.out file..." << endl;
    exit(0);
  }
  while(tmpString != "H-representation")
    in2 >> tmpString;
  in2 >> tmpString;
  if(tmpString == "linearity")
    { 
      in2 >> numOfEqu2; 
      NewIndex = new int[numOfEqu2];
      for(int i = 0; i < numOfEqu2; i++)  NewIndex[i] = 0;
      for(int i = 0; i < numOfEqu2; i++) {
	// NewIndex[i] = 0; 
	in2 >> NewIndex[i];
	//  cout << NewIndex[i];
      }
      
  }
    in2.close();
  if(flag == 2) strcpy(nonneg, "no");
  tmpString = "a";
  if(numOfEqu2 == 0) equ[0] = 'n';
  ifstream in3("Check_red.out");
  if(!in3){
    cerr << "Missing Check_red.out file..." << endl;
    exit(0);
  }
  while (tmpString!="begin") getline(in3,tmpString);

  in3 >> numOfConsts2 >> numOfDims >> tmpString;

   ZZ newEnt[numOfConsts2][numOfDims];
    for(int i = 0; i < numOfConsts2; i++)
      for(int j = 0; j < numOfDims; j++){
	in3 >> newEnt[i][j]; }
   
    strcat(Filename, ".latte");

    ofstream out2(Filename);

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
    system("rm Check_red*");
    cout << "done." << endl;
}

/* ----------------------------------------------------------------- */
void CheckInputFileCDDRep1(char *InputFile){
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
    //  cout << tmpString << endl;
    counter++;}
  //cout << counter << endl;

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
      in4 >> tmp_int; //cout << tmpString2 << " " << tmp_int << endl;
      for(int i = 0; i < tmp_int; i++) in4 >> num;
      in4 >> tmpString2;//cout << tmpString2 << endl; */
      if(tmpString2  != "begin") flag = 1;
      in4 >> tmpString2;//cout << tmpString2 << endl;
      in4 >> tmpString2;//cout << tmpString2 << endl;
      in4 >> tmpString2;//cout << tmpString2 << endl;
      if(tmpString2  != "integer") flag = 1;
      while(tmpString2 != "linearity"){ in4 >>tmpString2; } 
      in4 >> tmp_int; //cout << tmpString2 << " " << tmp_int << endl;
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
  // cout << flag3 << endl;

  if(flag == 1){
    ofstream out("Error");
    out << "Your input file CDD version is not correct!" << endl;
    cerr << "Your input file CDD version is not correct!" << endl;
    exit (0);
  }
  //  delete [] tmp;
  return ;
}


/* ----------------------------------------------------------------- */
listCone* ProjectUp(listCone* cone, int & oldNumOfVars, int & newNumOfVars, 
             listVector *equations){

  listCone *current_cone = cone;
  vec_ZZ newVector;

  newVector.SetLength(oldNumOfVars);
 
  listVector *temp, *temp2, *current_ray, *new_ray;
  int i;

  while(current_cone)
  {

    temp2 = equations;
    //  cout << " Here 1" << endl;
    i = 0;
    while(temp2)
      {	  
	newVector[i] = temp2->first * current_cone->latticePoints->first;
	temp2 = temp2->rest;
	i++;
      }
    //  cout << " Here 2" << endl;
    for(i = oldNumOfVars - newNumOfVars; i < oldNumOfVars; i++)
      {
	newVector[i] = current_cone->latticePoints->first[i - oldNumOfVars + newNumOfVars];
      }
    // cout << " Here 3" << endl;
    delete current_cone->latticePoints;
    current_cone->latticePoints = new listVector;
    current_cone->latticePoints->rest = NULL;

    current_cone->latticePoints->first.SetLength(oldNumOfVars);
    //  cout << " Here 4" << endl;
    for(i = 0; i < oldNumOfVars; i++)
      current_cone->latticePoints->first[i] = newVector[i];

    current_ray = current_cone->rays;
    new_ray = new listVector;
    current_cone->rays = new_ray;
    //   cout << " Here 5" << endl;
    while(current_ray)
      {
	temp2 = equations;

	i = 0;
	while(temp2)
	  {	  
	    newVector[i] = temp2->first * current_ray->first;
	    temp2 = temp2->rest;
	    i++;
	  }

	for(i = oldNumOfVars - newNumOfVars; i < oldNumOfVars; i++)
	  {
	    newVector[i] = current_ray->first[i - oldNumOfVars + newNumOfVars];
	  }
	
	temp = current_ray;
	current_ray = current_ray->rest;
	delete temp;
	
	new_ray->first.SetLength(oldNumOfVars);
	//   cout << " Here 6" << endl;
	for(i = 0; i < oldNumOfVars; i++)
	  new_ray->first[i] = newVector[i];
      

	if(current_ray != NULL)
	  {
	    new_ray->rest = new listVector;
	    new_ray = new_ray->rest;
	  }
	else
	  new_ray->rest = NULL;
      }

    current_cone = current_cone->rest;
  }
  return cone;
}

/* ----------------------------------------------------------------- */
listCone* ProjectUp2(listCone* cone, int & oldNumOfVars, int & newNumOfVars, 
             mat_ZZ AA, vec_ZZ b){

  // d =  oldNumOfVars and k = newNumOfVars
  
  listCone *current_cone = cone;
  vec_ZZ newVector;
  
  newVector.SetLength(oldNumOfVars);
  
  listVector *temp, *current_ray, *new_ray;
  int i;
  
  while(current_cone)
    {
      
      //  cout << " Here 1" << endl;
      i = 0;
      newVector = b;
      
      for(i = 0; i < oldNumOfVars; i++){	  
	newVector[i] += AA[i] * current_cone->latticePoints->first;
      }
      
      //  cout << " Here 2" << endl;
      /*    for(i = oldNumOfVars - newNumOfVars; i < oldNumOfVars; i++)
	    {
	    newVector[i] = current_cone->latticePoints->first[i - oldNumOfVars + newNumOfVars];
	    }*/
	// cout << " Here 3" << endl;
      delete current_cone->latticePoints;
      current_cone->latticePoints = new listVector;
      current_cone->latticePoints->rest = NULL;
      
      current_cone->latticePoints->first.SetLength(oldNumOfVars);
      //  cout << " Here 4" << endl;
      for(i = 0; i < oldNumOfVars; i++)
	current_cone->latticePoints->first[i] = newVector[i];
      
      current_ray = current_cone->rays;
      new_ray = new listVector;
      current_cone->rays = new_ray;
      //   cout << " Here 5" << endl;
      while(current_ray)
	{
	  i = 0;
	  for(i = 0; i < oldNumOfVars; i++)
	    {	  
	      newVector[i] = AA[i] * current_ray->first;
	    }
	  
	  // 	for(i = oldNumOfVars - newNumOfVars; i < oldNumOfVars; i++)
	  // 	  {
	  // 	    newVector[i] = current_ray->first[i - oldNumOfVars + newNumOfVars];
	  // 	  }
	  
	  temp = current_ray;
	  current_ray = current_ray->rest;
	  delete temp;
	  
	  new_ray->first.SetLength(oldNumOfVars);
	  //   cout << " Here 6" << endl;
	  for(i = 0; i < oldNumOfVars; i++)
	    new_ray->first[i] = newVector[i];
	  
	  
	  if(current_ray != NULL)
	    {
	      new_ray->rest = new listVector;
	      new_ray = new_ray->rest;
	    }
	  else
	    new_ray->rest = NULL;
	}
      
      current_cone = current_cone->rest;
    }
  return cone;
}
  
  
/* ----------------------------------------------------------------- */
void CheckInputFileCDDRep(char *InputFile){
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
	    exit (0);
	    }*/
          else ;
     
    }

  // cout << flag << endl;
  if(flag != 3){
    ofstream out("Error");
    out << "Your input file is not correct!" << endl;
    out << "Must be H-representation with integer!" << endl;
    cerr << "Your input file is not correct!" << endl;
    cerr << "Must be H-representation with integer!" << endl;
    exit (0);
  }
  //  delete [] tmpString;
  return ;
}

/* ----------------------------------------------------------------- */

void CheckInputFileCDDRep3(char *InputFile){
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
    {in >> tmpString;// cout << tmpString << endl;
    counter ++;
    }

  if(counter != numOfConst * dim + 1) flag = 1;

  if(flag == 1){
    ofstream out("Error");
    out << "Your input file has wrong number of elements!" << endl;
    cerr << "Your input file has wrong number of elements!" << endl;
    exit (0);
  }

  return ;
}

/* ----------------------------------------------------------------- */
void CheckInputFileCDDRep4(char *InputFile){
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
    exit (0);
  }
  delete [] tmpString;
  return ;
}


/* ----------------------------------------------------------------- */
void CheckInputFile(char *InputFile){
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
    exit (0);
  }
  delete [] tmpString;
  return ;
}
/* ----------------------------------------------------------------- */
void CheckInputFileVrep(char *InputFile){
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
    exit (0);
  }
  delete [] tmpString;
  return ;
}

/* ---------------------------------------------------------------- */

void CheckLength(char * filename, char * equ)
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
     exit (0);
   }
   /*   else if(Total < counts){
     ofstream out("Error");
     out << "The number of elements are more than you indicated. " << endl;
     cerr << "The number of elements are more than you indicated. " << endl;
     exit (0);
  
     }*/
}

/* ---------------------------------------------------------------------- */

void CheckLength2(char * filename, char * equ)
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
     exit (0);
   }
   /*   else if(Total < counts){
     ofstream out("Error");
     out << "The number of elements are more than you indicated. " << endl;
     cerr << "The number of elements are more than you indicated. " << endl;
     exit (0);
  
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
void readLatteProblem(char *fileName, listVector **equations,
		      listVector **inequalities, 
		      char *equationsPresent,
                      int *numOfVars, char *nonneg, char* dual,
		      char* grobner, char* max, vec_ZZ & cost, char * Vrep) {
  int i,j,eq,ind,numOfVectors,numOfEquations;
  vec_ZZ indexEquations, tmpVector;
  listVector *basis, *endBasis, *tmp, *endEquations, *endInequalities;
  vec_ZZ b;
  ZZ bignum;
  /* Reads numOfVars, matrix A, and rhs b. */


  cout << "Reading problem.\n";

  setbuf(stdout,0);

  ifstream in(fileName);
  if(!in){
    cerr << "Cannot open input file in readLatteProblem." << endl;
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

    B = transpose(A); //cout << B << endl;
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
  in >> hold;// cout << hold << endl;
  for (j=1; j<(*numOfVars)-1; j++) in >> b[j];
  b[*numOfVars-1] = hold;// cout << b << endl;
  basis = createListVector(b);
  endBasis = basis;

  for (i=1; i<numOfVectors; i++) {
    b=createVector(*numOfVars);
    in >> hold;//cout << hold << endl;
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
    cout << "\nEquation indices: ";
    printVector(indexEquations,numOfEquations);

    tmpVector=createVector(*numOfVars);
    createListVector(tmpVector);
    (*equations)=createListVector(tmpVector);

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
  //if(max[0] == 'y') for(i = 0; i < (*numOfVars-1); i++) in >> cost[i];
  if(Vrep[0] == 'n'){
  cout << endl;
  cout << "Ax <= b, given as (b|-A):\n";
  cout << "=========================\n";
  printListVector(*inequalities,*numOfVars);
  
  cout << endl;

  cout << "Ax = b, given as (b|-A):\n";
  cout << "========================\n";
  printListVector(*equations,*numOfVars);

  cout << endl;}

  else{
  cout << endl;
  cout << "The vertex Set, given by (1| V):\n";
  cout << "=========================\n";
  printListVector(*inequalities,*numOfVars);
  
  cout << endl;
  }
  }
  return;
}

/* ----------------------------------------------------------------- */
int CDDstylereadLatteProblem(char *fileName, listVector **equations,
		      listVector **inequalities, 
		      char *equationsPresent,
                      int *numOfVars, char *nonneg, char* dual,
                      char* taylor, int & degree, 
                      char* rat, int & cone_output, int & flags,
                      char* Memory_Save, char* uni, char* inthull,
		      char* grobner) {
  int i,j,eq,ind, length = 0, f = 0, numOfVectors,numOfEquations;
  vec_ZZ indexEquations, tmpVector;
  listVector *basis, *endBasis, *tmp, *endEquations, *endInequalities;
  vec_ZZ b;
  string tmpString;
  
  if(grobner[0] == 'y') strcpy(equationsPresent, "yes");
  if(dual[0] == 'y') flags |= DUAL_APPROACH;
  if(rat[0] == 'y') strcpy(rat, "yes");
  /* Reads numOfVars, matrix A, and rhs b. */

  ifstream in2(fileName);
  if(!in2){
    cerr << "Cannot open input file in readLatteProblem." << endl;
    exit(1);
  }

  while(in2 >> tmpString){ //cout << tmpString << endl;

    if(tmpString == "dual") {strcpy(dual, "yes");
    flags |= DUAL_APPROACH;}
    else if(tmpString == "nonneg") strcpy(nonneg, "yes");
    else if(tmpString == "taylor") 
        {strcpy(taylor, "yes");
	in2 >> degree; //cout << degree << endl;
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
  cout << "Reading problem.\n";
  //  cout << degree << endl;
  setbuf(stdout,0);
  ZZ bignum;
  ifstream in(fileName);
  if(!in){
    cerr << "Cannot open input file in readLatteProblem." << endl;
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
      //     cout << length << endl;
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

    B = transpose(A); cout << B << endl;
    for(i = 0; i < *numOfVars; i++)
    { ZZ tmp_ZZ;
    tmp_ZZ = B[i]*B[i];  //cout << B[i] << endl;
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
    if(bignum < tmp_ZZ)  bignum = tmp_ZZ; //cout << bignum << endl;
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
  in >> hold;// cout << hold << endl;
  for (j=1; j<(*numOfVars)-1; j++) in >> b[j];
  b[*numOfVars-1] = hold;// cout << b << endl;
  basis = createListVector(b);
  endBasis = basis;

  for (i=1; i<numOfVectors; i++) {
    b=createVector(*numOfVars);
    in >> hold;//cout << hold << endl;
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
    cout << "\nEquation indices: ";
    printVector(indexEquations,numOfEquations);

    tmpVector=createVector(*numOfVars);
    createListVector(tmpVector);
    (*equations)=createListVector(tmpVector);

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

  cout << endl;
  cout << "Ax <= b, given as (b|-A):\n";
  cout << "=========================\n";
  printListVector(*inequalities,*numOfVars);
  
  cout << endl;

  cout << "Ax = b, given as (b|-A):\n";
  cout << "========================\n";
  printListVector(*equations,*numOfVars);

  cout << endl;

  return 0;
  }
}


/* ----------------------------------------------------------------- */
