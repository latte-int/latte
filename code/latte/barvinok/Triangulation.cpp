/***********************************************************************************
   Author: Ruriko Yoshida
   Date: August 18th, 2002
   Update: November 25th, 2002
   Project for LattE
   This program computes a triangulation of a cone
   Input: the number of vectors, the dimension of the vectors and
          the cone generators.
   Output: the cone triangulation in R^n.
************************************************************************************/
//  #ifndef TRIANGULATION__H
//  #define TRIANGULATION__H
#include <iostream>
#include <fstream>
#include <cctype>
#include <cstring>
#include <cstdlib> // exit()
#include <string>
#include <math.h>
#include "latte_ntl.h"
#include "Triangulation.h"
#include "../flags.h"
#include <string.h>
#include <stdio.h>
#include "config.h"
#include "latte_system.h"

using namespace std;

/*
  This programs uses the following method.  The first step is to set up 
  the following vector configuration. For each p in the vertex set, set 
  0 p_1 p_2 ... p_d r_i where r_i is a random rational number for each 
  row of the input matrix, then excute CDD+.  Read the output file *.ine, 
  if there is negative component at the end of each inequality.  If it is, 
  consider the corresponding plane. Then, list up all vectors in the plane.  
  The family of correction of vectors for each such a plane is 
  triangulation of the input cone.
*/

/* ----------------------------------------------------------------- */
void writeCDDextFileRudy(const int & m, const int & n, const mat_ZZ & Mat) {
   /*
     Create an input file filename.ext for CDD+
   */
   ofstream OUT;
   OUT.open("tri.ext");
   OUT << "V-representation" << endl;
   OUT << "begin " << endl;
   OUT << m << " " << n + 2 << " rational" << endl;
   for(int i = 0; i < m; i++){
      OUT << 0 << " " ;
      for(int j = 0; j < n; j++){
         OUT << Mat[i][j] << " ";
         }
      OUT << rand() % 100 + 1 << "/" << rand() % 2000 + 10 << "\n" ;

      }
   OUT << "end" << endl;
   OUT << "incidence" << endl;
   OUT << "input_adjacency" << endl;
   OUT.close();

  return ;
}
/* ----------------------------------------------------------------- */
vec_ZZ readCDDineFileRudy(int& face3, int& dim3) {
   ifstream IN;
   IN.open("tri.ine");

   if(!IN)
    {
     cerr<<"File could not be opened in readCDDineFileRudy. "<<endl;
      exit(2);
    }

   string trash2;
   while(trash2 != "begin")
     getline(IN, trash2);

   IN >> face3 >> dim3;

   vec_ZZ neg;
   neg.SetLength(face3);
   for( int i = 0; i < face3; i++) neg[i] = 0;
   /*
     Reading filename.ine file.  Check each plane and if the last 
     coordinate of inequality is negative, we keep track and assign 
     neg[i] = 1.  We need this process in order to pick right 
     triangulation later.
   */
   IN >> trash2;
   for(int i = 0; i < face3; i++){
     for(int j = 0; j < dim3; j++)
       IN >> trash2;
       if(trash2[0] == '-')
         neg[i] = 1;
       if(trash2[0] == '0')
         neg[i] = 2;
   }

  return(neg); 
}
/* ----------------------------------------------------------------- */
void readCDDicdFileRudy(int & face2, vec_ZZ & numOfPoints, mat_ZZ & Result) {
  ifstream File3;
  File3.open("tri.icd");
  /*
    Reading filename.icd file.  This file has the list of points 
    in each plane.
  */
  if(!File3)
    {
      cerr<<"File could not be opened in readCDDicdFileRudy."<<endl;
      exit(5);
    }
  int num, art;
  string tmp;
  
  while(tmp != "begin")
    getline(File3, tmp);
  
  File3 >> face2 >> num >> art; 
  Result.SetDims(face2, num);
  numOfPoints.SetLength(face2);
  char trash;
 
  int number = 0;
  
  /*
    If numOfPoints[i] is negative, then, numOfPoints[i]
    is the number of points which are NOT on the plane and
    filename.icd has a list of points which
    are NOT in the plane.  
    In that case, we want to have points which are NOT in the 
    list.
  */
  for(int i = 0; i < face2; i++){
    File3 >> numOfPoints[i] >> trash;
    if(numOfPoints[i] > 0)
      for(int j = 0; j < numOfPoints[i]; j++)
	File3 >> Result[i][j];
    else { 
      number = 0;
      numOfPoints[i] = - numOfPoints[i];
      for(int j = 0; j < numOfPoints[i]; j++)
	File3 >> Result[i][j];
      
      numOfPoints[i] = num - numOfPoints[i];
      vec_ZZ Temp;
      Temp.SetLength(num + 1);
      
      for(int j = 0; j <= num; j++) conv(Temp[j], j+1);
      for(int j = 0; j < num - numOfPoints[i]; j++) {
	for(int k = 0; k <= num; k++) {
	  if(Temp[k] == Result[i][j]) {
	    Temp[k] = 0;   
	  }
	}
      }
      for(int j = 0; j <= num; j++) {
	if(Temp[j] != 0)
	  Result[i][j - number] = Temp[j];
	else
	  number++;
      }
    }
  }
  
  File3.close();
  
  return ; 
}
/* ----------------------------------------------------------------- */
int Triangulation(const mat_ZZ & Mat, const int & m, const int & n, 
		  char* a, list< int >& List)
{
  vec_ZZ neg;
  
  writeCDDextFileRudy(m,n,Mat);
  system_with_error_check(CDD_PATH " tri.ext > tri.out");
  system_with_error_check("rm -f tri.out");
  
  int face3 = 0, dim3 = 0;
  neg = readCDDineFileRudy(face3, dim3);
  
  vec_ZZ numOfPoints;
  mat_ZZ Result;
  int face2 = 0;

  readCDDicdFileRudy(face2,numOfPoints,Result);

  int count = 0;
  int ZZToInt = 0;
  for(int i = 0; i < face2; i++){
    if(neg[i] == 1 ){
      for(int j = 0; j < numOfPoints[i]; j++)
        if(Result[i][j] != 0){
          conv(ZZToInt, Result[i][j]);
	  List.push_front(ZZToInt);
	}
      count++;
    }
  }
  if( count == 0)
    for(int i = 0; i < face2; i++){
      if(neg[i] == 0 ){
	for(int j = 0; j < numOfPoints[i]; j++)
	  if(Result[i][j] != 0){
            conv(ZZToInt, Result[i][j]);
            List.push_front(ZZToInt);
	  }
	count++;
      }
    }

  //system_with_error_check("rm -f tri.*");

  return count;
}
/* ----------------------------------------------------------------- */



// NEW TRIANGULATION FUNCTION TO REUSE PREVIOUSLY CALCULATED RESULTS

/* ----------------------------------------------------------------- */
int Triangulation_Load_Save (const mat_ZZ & Mat, const int & m, const int & n, char* a, list< int >& List, char *File_Name, int Cone_Index, unsigned int Flags) 
{
//	cout << "Triangulation_Load_Save: Cone_Index: " << Cone_Index << "  Flags: " << Flags << endl;
  vec_ZZ neg; 
 	int	File_Present = 1;	 
	char		File_Path [256];
	char		System_Command[256];
	char		*Integer_String;
	char 		*Integer_String_Original;  
	
	Integer_String = new char [100];
	Integer_String_Original = Integer_String;
	
	sprintf (Integer_String, " %d", Cone_Index ); 
	
	while (Integer_String[0] == ' ')
		Integer_String += 1;		
	
	strcpy		(File_Path, "triangulations/");
	strcat		(File_Path, File_Name );
	strcat		(File_Path, Integer_String ); 
	strcat		(File_Path, ".tar.gz");
	
	if ( Flags & LOAD)
	{
		//check to see if file is present   [File_Name].tar.gz
		
		ifstream	checkforfile;
		

		cout << "Triangulation_Load_Save: Checking for " << File_Path << endl;

		checkforfile.open (File_Path);

		if (checkforfile.fail ())
		{
			//file not present
			File_Present = 0;

			cout << "Triangulation_Load_Save: File not present.  Calculation triangulation." << endl;
  			writeCDDextFileRudy(m,n,Mat);
  			system_with_error_check(CDD_PATH "tri.ext > tri.out");
  			system_with_error_check("rm -f tri.out");
		}	
		else
		{
			//unzip file, rename [File_name][Cone_Index] to tri.ext
			strcpy (System_Command, "tar -zxf ");
			strcat (System_Command, File_Path );

			cout << "Triangulation_Load_Save: File present.  Unziping using command: " << System_Command << endl;
	
			system_with_error_check (System_Command);
		}		
	}
	else // NOT LOADING
	{

  			writeCDDextFileRudy(m,n,Mat);
  			system_with_error_check(CDD_PATH " tri.ext > tri.out");
  			system_with_error_check("rm -f tri.out");
	}	
  
  int face3 = 0, dim3 = 0;
  
  /////  read stuff from file
  neg = readCDDineFileRudy(face3, dim3);
  
  vec_ZZ numOfPoints;
  mat_ZZ Result;
  int face2 = 0;
  /////// read stuff from file
  readCDDicdFileRudy(face2,numOfPoints,Result);

  int count = 0;
  int ZZToInt = 0;
  for(int i = 0; i < face2; i++){
    if(neg[i] == 1 ){
      for(int j = 0; j < numOfPoints[i]; j++)
        if(Result[i][j] != 0){
          conv(ZZToInt, Result[i][j]);
	  List.push_front(ZZToInt);
	}
      count++;
    }
  }
  if( count == 0)
    for(int i = 0; i < face2; i++){
      if(neg[i] == 0 ){
	for(int j = 0; j < numOfPoints[i]; j++)
	  if(Result[i][j] != 0){
            conv(ZZToInt, Result[i][j]);
            List.push_front(ZZToInt);
	  }
	count++;
      }
    }
	
  	if ( (Flags & SAVE)  && ( !(Flags & LOAD) || File_Present == 0)  )
	{
		//save files to triangulation/[File_Name][Cone_Index].tar.gz 	

		strcpy (System_Command, "tar -zcf ");
		strcat (System_Command, File_Path);
		strcat (System_Command, " tri.ine tri.icd");

		cout << "Triangulation_Load_Save: Save mode.  Creating archive of tri.ine tri.icd with command: " << System_Command << endl;

		system_with_error_check (System_Command);
		
	}
  
  system_with_error_check("rm -f tri.*");


  delete [] Integer_String_Original;
  return count;
}
