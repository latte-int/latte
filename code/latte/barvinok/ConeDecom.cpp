/************************************************************
  Author: Ruriko Yoshida
  Date: August 24th, 2002
  Update: Febrary 3rd, 2003
  This is a main file for ConeDecom.
  Input: a matrix whose columns represent the generators of a cone.
  Output: Unimodular cones.
**************************************************************/

#include <list>
#include <fstream.h>
#include <cstdlib>
#include <cstring>
#include <string>

#include "../myheader.h"
#include "../ramon.h"
#include "../print.h"
#include "../cone.h"
#include "Cone.h"
#include "barvinok.h"
#include "Triangulation.h"
#include <NTL/mat_ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_RR.h>
#include <NTL/vec_RR.h>
#include <NTL/config.h>
#include <NTL/LLL.h>
#include <NTL/HNF.h>
#include <NTL/ZZ.h>

/* ----------------------------------------------------------------- */
listVector* transformArrayBigVectorToListVector(mat_ZZ A, int numOfVectors, 
						int numOfVars) {
  int i;
  vec_ZZ v;
  listVector *L, *endL;

  v=createVector(numOfVars);
  L=createListVector(v);
  endL=L;

  for (i=0; i<numOfVectors; i++) {
    v=A[i];
    endL->rest = createListVector(v);
    endL = endL->rest;
  }

  return (L->rest);
}

/* ----------------------------------------------------------------- */
listCone* transformRudyListConeIntoRamonListCone(list< PtrCone > RudyCones,
						 int numOfVars) {
  int s;
  listCone *cones, *endCones, *newCone;
  PtrCone tmp;

  cones=createListCone();
  endCones=cones;

  while (!RudyCones.empty()) {
    newCone=createListCone();

    tmp = RudyCones.back();
    RudyCones.pop_back();
    s = tmp.sign;
    if (s==0) s=-1;

    newCone->coefficient=s;
    newCone->rays = tmp.Generator;
    endCones->rest=newCone;
    endCones=endCones->rest;
   // tmp.generator.kill();
  }
//    printListCone(cones->rest,numOfVars);
  return (cones->rest);
}
/* ----------------------------------------------------------------- */

/*
  The first step is to triangulate a cone into simplicial cones.
  Then, by using Barvinok's decomposition, we decompose each
  simplicial cone into unimodular cones.
*/

/* ----------------------------------------------------------------- */
listCone* barvinokDecomposition(mat_ZZ Mat, int m, int n, int & numOfUniCones, char *File_Name, unsigned int Flags, int Cone_Index) {

  list< PtrCone > Uni;

  /* m is the number of vectors and n is the number of dims. */
  if((m == 0) || (n == 0)){
    cerr << "The polytope is empty!" << endl;
    //    system("rm latte_dec");
    //system("rm core");
    exit(2);
  }
//      mat_ZZ Mat;

//      Mat.SetDims(m, n);

//      for(int i= 0; i< m; i++)
//        for(int j = 0; j < n; j++)
//          conv(Mat[i][j], MatRays[i][j]);
    
  int Face = 1, Faces = 10000;
  char* s1 = "latte_dec";
  list< int > List;
  if(m != n){
    Face = Triangulation_Load_Save(Mat, m, n, s1, List, File_Name, Cone_Index, Flags);
  } /*Call triangulation fun.*/

   /*
     In this fun, a cone is decomposed into simplicial cones.
   */
  Faces = Face;
  mat_ZZ B[Faces];

  for(int i = 0; i < Faces; i++)
    B[i].SetDims(n, n);

  if(m != n){
    long tmp = 0;
    int counter = 0, index = 0;

    while(!List.empty())
     {
       tmp = List.back();
       List.pop_back();
       B[index]((counter % n) + 1) = Mat(tmp);
       counter++;
       if((counter % n == 0))
         index++;
       
     } 
  }
  if(m == n) 
  {
	  B[0] = Mat;
  }
   /*
     Call barvinok fun.  barvinok function decompose each
     simplicial cone into unimodular cones.

   */
   int total = 0;
    for(int i = 0; i < Faces; i++){
      if(IsZero(B[i]) != 1){
        total += barvinok(B[i],  Uni, numOfUniCones);
        }
      }

  for(int i = 0; i < Faces; i++)
    B[i].kill();

   return (transformRudyListConeIntoRamonListCone(Uni,n));
}



/* ----------------------------------------------------------------- */
int barvinokDecomposition_Single(const mat_ZZ Mat, int m, int & numOfUniCones, rationalVector *vertex, Single_Cone_Parameters *Parameters, Node_Controller *Controller, char *File_Name, int Cone_Index) 
{

  list< PtrCone > Uni;

  /* m is the number of vectors and n is the number of dims. */

  int n = Parameters->Number_of_Variables;
  
  //cout << "BarvinokDecomposition_Single Called:  " << endl;
 
  if((m == 0) || (n == 0)){
    cerr << "The polytope is empty!" << endl;
    //    system("rm latte_dec");
    //system("rm core");
    exit(2);
  }
//      mat_ZZ Mat;

//      Mat.SetDims(m, n);

//      for(int i= 0; i< m; i++)
//        for(int j = 0; j < n; j++)
//          conv(Mat[i][j], MatRays[i][j]);
    
  int Face = 1, Faces = 10000;
  char* s1 = "latte_dec";
  list< int > List;
  if(m != n){
    Face = Triangulation_Load_Save(Mat, m, n, s1, List, File_Name, Cone_Index, Parameters->Flags);
  } /*Call triangulation fun.*/

   /*
     In this fun, a cone is decomposed into simplicial cones.
   */
  Faces = Face;
  mat_ZZ B[Faces];

  for(int i = 0; i < Faces; i++)
    B[i].SetDims(n, n);

  if(m != n){
    long tmp = 0;
    int counter = 0, index = 0;

    while(!List.empty())
     {
       tmp = List.back();
       List.pop_back();
       B[index]((counter % n) + 1) = Mat(tmp);
       counter++;
       if((counter % n == 0))
         index++;
       
     } 
  }
  if(m == n) B[0] = Mat;

   /*
     Call barvinok fun.  barvinok function decompose each
     simplicial cone into unimodular cones.

   */
   //cout << "BarvinokDecomposition: Number of simplicies = " << Faces << endl;
    for(int i = 0; i < Faces; i++){
      if(IsZero(B[i]) != 1){
        if ( barvinok_Single(B[i], numOfUniCones, Parameters, Controller, vertex) == -1)
	{
  		for(int i = 0; i < Faces; i++)
			B[i].kill ();
		return -1;
	}
	//cout << "BarvinokDecomposition: Current total lattice points (not div by ten pow): " << *(Parameters->Total_Lattice_Points) << endl;
	
      	//cout << "BarvinokDecomposition: Taylor Expansion: ";  
	
	//for (int j = 0; j <= Parameters->Degree_of_Taylor_Expansion; j++)
	//	cout << Parameters->Taylor_Expansion_Result[j] << "t^" << j << " ";
        //cout << endl;	
      	}
      }

  for(int i = 0; i < Faces; i++)
    B[i].kill();

  return 1;
}













