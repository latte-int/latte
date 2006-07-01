/************************************************************
  Author: Ruriko Yoshida
  Date: August 24th, 2002
  Update: Febrary 3rd, 2003
  This is a main file for ConeDecom.
  Input: a matrix whose columns represent the generators of a cone.
  Output: Unimodular cones.
**************************************************************/

#include <list>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <string>

#include "../myheader.h"
#include "../ramon.h"
#include "../print.h"
#include "../cone.h"
#include "convert.h"
#include "barvinok.h"
#include "Triangulation.h"

/*
  The first step is to triangulate a cone into simplicial cones.
  Then, by using Barvinok's decomposition, we decompose each
  simplicial cone into unimodular cones.
*/

/* ----------------------------------------------------------------- */
int barvinokDecomposition_Single(const mat_ZZ &Mat, rationalVector *vertex,
				 Single_Cone_Parameters *Parameters) 
{
  /* m is the number of vectors and n is the number of dims. */

  int m = Mat.NumRows();
  int n = Parameters->Number_of_Variables;
  assert(Mat.NumCols() == n);
  
  if((m == 0) || (n == 0)){
    cerr << "The polytope is empty!" << endl;
    exit(2);
  }
  int Face = 1, Faces = 10000;
  char* s1 = "latte_dec";
  list< int > List;
  if(m != n){
    Parameters->triangulate_time.start();
    Face = Triangulation_Load_Save(Mat, m, n, s1, List, Parameters->File_Name, Parameters->Cone_Index, Parameters->Flags);
    Parameters->triangulate_time.stop();
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
    for(int i = 0; i < Faces; i++){
      if(IsZero(B[i]) != 1){
        if ( barvinok_Single(B[i], Parameters, vertex) == -1)
	{
  		for(int i = 0; i < Faces; i++)
			B[i].kill ();
		return -1;
	}
      	}
      }

  for(int i = 0; i < Faces; i++)
    B[i].kill();

  return 1;
}

int
barvinokDecomposition_Single (listCone *cone,
			      Single_Cone_Parameters *Parameters)
{
  int numOfVars = Parameters->Number_of_Variables;
  int numOfRays = lengthListVector(cone->rays);
  mat_ZZ mat = createConeDecMatrix(cone,numOfRays,numOfVars);
  return barvinokDecomposition_Single(mat, cone->vertex,
				      Parameters);
}

