/*******************************************************************
   Author: Ruriko Yoshida
   Date: July 25th, 2002
   Update: Febrary 3rd, 2003
   This program computes Barvinok's decomposition of a cone.
   This program is for the project "LattE."

*********************************************************************/
#include <list>

#include <fstream>
#include <cstdlib>
#include <cstring>
#include <math.h>
#include <algorithm>
#include <time.h>

#include "Cone.h"
#include "barvinok.h"
#include "../myheader.h"
#include "../ramon.h"
#include "../dual.h"
#include "../RudyResNTL.h"
#include "../print.h"
#include "../genFunction/piped.h"

#include <NTL/mat_ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_RR.h>
#include <NTL/vec_RR.h>
#include <NTL/config.h>
#include <NTL/LLL.h>
#include <NTL/HNF.h>
#include <NTL/ZZ.h>
using namespace std;

 /* Note:  We are dealing with the "Row space" of the
    input matrix due to NTL. */

/**********************************************************************/
vec_ZZ CheckOmega( const mat_ZZ & U, vec_ZZ & Z){

  int m;
  m = U.NumCols();
  mat_ZZ Dummy;
  Dummy.SetDims(m + 1, m);
  Dummy[0] = Z;
  ZZ d;

  for(int i = 0; i < m; i++)
    Dummy[i + 1] = U[i];

  mat_ZZ dummy;
  image(d, Dummy, dummy);

  int flag = 1, number = 0; 

  for(int i = 0; i <= m; i++)
      if(dummy[0][i] >= 0) number++;
  if(number == (m + 1))  flag = 0;

  if(flag != 0){
  number = 0;
  for(int i = 0; i <= m; i++)
     if(dummy[0][i] <= 0) number++;
  if(number == (m + 1)) flag = 0;
  }     
  if(flag == 0){
    Z = - Z; 
  }
  Dummy.kill();
  dummy.kill();
  return Z;

}
/**********************************************************************/
 
void MatrixGCD(mat_ZZ & B, long & m){
  ZZ gcds[m];
  for(int i = 1; i <= m; i++)
    for(int j = 1; j <= m; j++)
      if(B(i, j) != 0)
	gcds[i-1] = GCD(gcds[i-1], B(i, j));
  for(int i = 1; i <= m; i++)
    for(int j = 1; j <= m; j++)
      if(B(i, j) != 0)
	B(i, j) = B(i, j) / gcds[i-1];

}
/**********************************************************************/
void AssignSign(const Cone& tmp, Cone & cones){ 

  ZZ Det = determinant(tmp.generator);
  if ((tmp.sign)>0) {
    if((Det * determinant(cones.generator)) >= 0)
      cones.sign = 1;
    else
      cones.sign = 0;
  } else {
    if((Det * determinant(cones.generator)) >= 0)
      cones.sign = 0;
    else
      cones.sign = 1;
  }

}

void AssignSign_Single( Cone *tmp, Cone *cones){ 

  ZZ Det = determinant(tmp->generator);
  if ((tmp->sign)>0) {
    if((Det * determinant(cones->generator)) >= 0)
      cones->sign = 1;
    else
      cones->sign = 0;
  } else {
    if((Det * determinant(cones->generator)) >= 0)
      cones->sign = 0;
    else
      cones->sign = 1;
  }
}

/**********************************************************************/
int barvinok(mat_ZZ & B, list< PtrCone > & Uni, int & numOfUniCones,
	     int max_determinant)
{
   long m, n;
  m = B.NumRows();


  n = B.NumCols();

  vec_ZZ v;
   v.SetLength(m);

   PtrCone tmpPtrCone;

   if( m != n)
   {
       cerr << "Input must be square. " << endl;
       exit(2);
   }

   ZZ D = determinant(B);

         if( D == 0)
   {
       cerr << "Input must be linearly independent. " << endl;
       exit(3);
   }

   vec_ZZ Z, Z2;
   ZZ Det;

   /* The following routine is to get the minimal
      integral generators for the cone.  */

   MatrixGCD(B, m);

   int  width = 0;
   list< Cone > QNonUni;

   Cone cones1[m];
   Cone dummy;
   dummy.generator = B;
   dummy.sign = 1;

   ZZ dummy_determinant = abs(determinant(dummy.generator));
   if(dummy_determinant == 1){
     listVector *L, *endL;

     L=createListVector(dummy.generator[0]);
     endL=L;

     for (int k=1; k<m; k++) {
       v=dummy.generator[k];
       endL->rest = createListVector(v);
       endL = endL->rest;
     }
     tmpPtrCone.Generator = L;
     tmpPtrCone.sign = 1;
     tmpPtrCone.determinant = dummy_determinant;
     Uni.push_back(tmpPtrCone);
   }

   else
     QNonUni.push_back(dummy);
   int counter = 0, count = 0;
   Cone tmp;
   
   while(!QNonUni.empty()){
     tmp = Equal(QNonUni.back());  
     QNonUni.pop_back();
     
     /* ComputeOmega(const mat_ZZ &, long& ) computes
	an integral vector in the parallelogram. */
     
     Z = ComputeOmega(tmp.generator, m, 0, 0);
     
     Z = CheckOmega(tmp.generator, Z);
     /* Copy the original matrix into m matrices. */
     
     for(int i = 0; i < m; i++){
       cones1[i].generator = tmp.generator;
     }
     /* Copy the vector which calculated into each ith row
	of each matrix (we are dealing with the row space). */
     
     for(int i = 1; i <= m; i++)
       for(int j = 1; j <= m; j++)
         cones1[i-1].generator(i, j) = Z(j);
     
     /* Compute a determinant of each matrices of an array.
	If the product of the determinant of the tmp and a determinant
	if each matrix is negative, we assign negative sign.  Otherwise,
	we assign positive to each matrix. */
    
     Det = determinant(tmp.generator);

     ZZ *Dets = new ZZ[m];
     for(int i = 0; i < m; i++)
       Dets[i] = determinant(cones1[i].generator);

     for(int i = 0; i < m; i++){
       if(abs(Dets[i]) >= abs(Det)){
         cout << "Second loop... " << endl;
	 Z = ComputeOmega(tmp.generator, m, 2, 2);
	 Z = CheckOmega(tmp.generator, Z);
	 for(int k = 1; k <= m; k++)
	   for(int j = 1; j <= m; j++)
	     cones1[k-1].generator(k, j) = Z(j);
       }
     }

     for(int i = 0; i < m; i++) {
       //This is AssignSign(tmp, cones1[i]): 
       if ((tmp.sign)>0) {
	 if((Det * Dets[i]) >= 0)
	   cones1[i].sign = 1;
	 else
	   cones1[i].sign = 0;
       } else {
	 if((Det * Dets[i]) >= 0)
	   cones1[i].sign = 0;
	 else
	   cones1[i].sign = 1;
       }
     }

#ifdef SHOWDETS
     cout << "Determinant: " << Det << " -> ";
     for(int i = 0; i < m; i++)
       cout << Dets[i] << " ";
     cout << endl;
#endif
     
     /* If a matrix is unimodular, then we put a matrix into QUni queue.
	Otherwise, we put it inot QNonUni queue. */
     
     for(int i = 0; i < m; i++) {
       ZZ Det_i = Dets[i];
       if (Det_i == 0) {
	 /* do nothing */
       }
       else if(abs(Det_i) <= max_determinant) {
	 listVector *L, *endL;

	 L=createListVector(cones1[i].generator[0]);
	 endL=L;

	 for (int k=1; k<m; k++) {
	   endL->rest = createListVector(cones1[i].generator[k]);
	   endL = endL->rest;
	 }
	 tmpPtrCone.Generator = L;
	 tmpPtrCone.sign = cones1[i].sign;
	 tmpPtrCone.determinant = abs(Det_i);
	 Uni.push_back(tmpPtrCone);
	 numOfUniCones++;
	 if((numOfUniCones % 1000) == 0)
	   cout << numOfUniCones << " unimodular cones are done. " << endl;
	 width++;
       }
       else if(abs(Det_i) < abs(Det))
	 QNonUni.push_back(cones1[i]);
     
       else
	 {
	   cerr << "Error!  We cannot have smaller determinant!" << endl;
	   exit(5);}
     }
     
      for(int i = 0; i < m; i++) cones1[i].generator.kill();
      delete[] Dets;
     Z.kill();
     Z2.kill();
     tmp.generator.kill();
     counter++;
     count = 0;  
   }
   for(int i = 0; i < m; i++) cones1[i].generator.kill();
   return width;
 }




int barvinok_Single(mat_ZZ & B, int & numOfUniCones, Single_Cone_Parameters *Parameters, Node_Controller *Controller, rationalVector *vertex)
{
	//cout << "barvinok_Single Called." << endl;;
	
	long m, n;
  	m = B.NumRows();


  	n = B.NumCols();

  	vec_ZZ v;
   	v.SetLength(m);


   	if( m != n)
   	{	
       		cerr << "Input must be square. " << endl;
       		exit(2);
   	}

   	ZZ D = determinant(B);

         if( D == 0)
   	{
       		cerr << "Input must be linearly independent. " << endl;
       		exit(3);
   	}

   	vec_ZZ Z, Z2;
   	ZZ Det;

   	/* The following routine is to get the minimal
      	integral generators for the cone.  */

   	MatrixGCD(B, m);


   	Cone *dummy = new Cone;
   	dummy->generator = B;
   	dummy->sign = 1;

	Barvinok_DFS_Parameters *DFS_Parameters = new Barvinok_DFS_Parameters;

	DFS_Parameters->vertex = vertex;
	DFS_Parameters->Controller = Controller;
	DFS_Parameters->Number_of_Variables = Parameters->Number_of_Variables;
	DFS_Parameters->Degree_of_Taylor_Expansion = Parameters->Degree_of_Taylor_Expansion;
	DFS_Parameters->Flags = Parameters->Flags;
	DFS_Parameters->Taylor_Expansion_Result = Parameters->Taylor_Expansion_Result;
	DFS_Parameters->Random_Lambda = Parameters->Random_Lambda;
	DFS_Parameters->Ten_Power = Parameters->Ten_Power;
	DFS_Parameters->Total_Lattice_Points = Parameters->Total_Lattice_Points;	
	DFS_Parameters->Total_Uni_Cones = Parameters->Total_Uni_Cones;    	
	DFS_Parameters->Current_Simplicial_Cones_Total = Parameters->Current_Simplicial_Cones_Total;
	DFS_Parameters->Max_Simplicial_Cones_Total = Parameters->Max_Simplicial_Cones_Total;
	
	
	int result;
	
	result = barvinok_DFS(dummy, DFS_Parameters);

	delete DFS_Parameters;
	dummy->generator.kill ();
	delete dummy;
	
	return result;
}
	
listCone* transformRudyListConeIntoRamonListCone_Single( PtrCone RudyCone,
						 int numOfVars) 
{
  int s;
  listCone *cones;

  	cones=createListCone();
  	
  	s = RudyCone.sign;
	if (s==0) s=-1;

	cones->coefficient=s;
	cones->determinant = RudyCone.determinant;
	cones->rays = RudyCone.Generator;
  return (cones);
}


int barvinok_DFS(Cone *C, Barvinok_DFS_Parameters *Parameters)
{
       	
    	ZZ Det = abs(determinant(C->generator));
    	int result = 1;
       	
     	if(Det == 0)
     	{
	 	//cout << "barvinok_DFS: Det = 0." << endl;
		return 1;	
   	  }		     
    	 else if(Det == 1)
    	 {
	 	//cout << "barvinok_DFS: Cone is unimod " << endl;
		
		*(Parameters->Total_Uni_Cones) += 1;

		if ( *(Parameters->Total_Uni_Cones) % 1000 == 0)
			cout << *(Parameters->Total_Uni_Cones) << " unimodular cones dones." << endl;
		 
		listVector *L, *endL;
   		vec_ZZ v;
		PtrCone tmpPtrCone;
	
		//cout << "barvinok_DFS: Setting NumRows...";	
		int m = C->generator.NumRows(); 
		//cout << "done." << endl;
	
		v.SetLength(m);
	
		//cout << "barvinok_DFS: Creating list... ";
		L=createListVector(C->generator[0]);
		//cout << "done." << endl;
  		endL=L;

 	 	for (int k=1; k<m; k++) 
		{
 	   		v=C->generator[k];
 	   		endL->rest = createListVector(v);
 	   		endL = endL->rest;
 	 	}
      	
		tmpPtrCone.Generator = L;
		tmpPtrCone.determinant = Det;
		
		// if something don't work, check this!
      		tmpPtrCone.sign = C->sign;

		listCone *Cone;

		//cout << "barvinok_DFS: Transforming Rudy list to Raymond list...";
		Cone = transformRudyListConeIntoRamonListCone_Single (tmpPtrCone, Parameters->Number_of_Variables );
		//cout << "done." << endl;
	
		Cone = dualizeBackCones (Cone, Parameters->Number_of_Variables );   	

		//cout << "barvinok_DFS: Calculating points in Parallelepiped" << endl;
		if ( (Parameters->Flags & DUAL_APPROACH) == 0)
			Cone->latticePoints = pointsInParallelepipedOfUnimodularCone (Parameters->vertex, Cone->rays, Parameters->Number_of_Variables);	
		
		Single_Cone_Parameters *Residue_Parameters = new Single_Cone_Parameters;

		Residue_Parameters->Cone = Cone;
		Residue_Parameters->Number_of_Variables = Parameters->Number_of_Variables;
		Residue_Parameters->Degree_of_Taylor_Expansion = Parameters->Degree_of_Taylor_Expansion;
		Residue_Parameters->Flags = Parameters->Flags;
		Residue_Parameters->Ten_Power = Parameters->Ten_Power;
		Residue_Parameters->Random_Lambda = Parameters->Random_Lambda;
		Residue_Parameters->Taylor_Expansion_Result = Parameters->Taylor_Expansion_Result;
	
		if (Parameters->Flags & DUAL_APPROACH)	
		{
			//cout << "barvinok_DFS: Calling ResidueFunction_Single_Cone" << endl;
			result =  ResidueFunction_Single_Cone ( Residue_Parameters, Parameters->Controller);
			
		}	
		else
		{
			result = Residue_Single_Cone (Cone, Parameters->Number_of_Variables, Parameters->Random_Lambda, Parameters->Total_Lattice_Points, Parameters->Ten_Power);
		
		}	

		//clean up
		//
		
		/*listVector *tempvector;
		while (Cone->rays)
		{
			tempvector = Cone->rays;
			Cone->rays = Cone->rays->rest;
			delete tempvector;
			
		}*/

		delete Residue_Parameters;
		
		
		return result;
  	} 		     
     
    
	//cout << "barvinok_DFS: non-uni cone." << endl;

     
     	vec_ZZ Z;    
     	long m = C->generator.NumRows();
     	//long n = C->generator.NumCols();

     	ZZ Dets[m];	     
     	Cone  *cones1 [m];
     
     	Z = ComputeOmega(C->generator, m, 0, 0);
     
     	Z = CheckOmega(C->generator, Z);
     
     	/* Copy the original matrix into m matrices. */
     
     	for(int i = 0; i < m; i++)
     	{
		cones1[i] = new Cone;	
      	 	cones1[i]->generator = C->generator;
     	}
     
     	/* Copy the vector which calculated into each ith row
	of each matrix (we are dealing with the row space). */
     
     	for(int i = 1; i <= m; i++)
       		for(int j = 1; j <= m; j++)
       		{       
         		cones1[i-1]->generator(i, j) = Z(j);
       		} 
     
     	/* Compute a determinant of each matrices of an array.
	If the product of the determinant of the tmp and a determinant
	if each matrix is negative, we assign negative sign.  Otherwise,
	we assign positive to each matrix. */
    

     
     
     	for(int i = 0; i < m; i++) AssignSign_Single(C, (cones1[i]));
     
     	for(int i = 0; i < m; i++)
     	{
       		if(abs(determinant(cones1[i]->generator)) >= Det)
		{
         		cout << "Second loop... " << endl;
	 		Z = ComputeOmega(C->generator, m, 2, 2);
	 		Z = CheckOmega(C->generator, Z);
	 		for(int k = 1; k <= m; k++)
	   			for(int j = 1; j <= m; j++)
	     				cones1[k-1]->generator(k, j) = Z(j);

	 		for(int s = 0; s < m; s++)  
				AssignSign_Single(C, (cones1[s]));
       		}
       
     	}

     	ZZ max;
     	max = -1;

#ifdef SHOWDETS
	cout << "Determinant " << Det << " -> ";
#endif
     	for(int i = 0; i < m; i++)
     	{
		Dets[i] = abs(determinant(cones1[i]->generator));
#ifdef SHOWDETS
		cout << Dets[i] << ", ";
#endif
	     	if(Dets[i] > max)
 			max = Dets[i];
     	}
#ifdef SHOWDETS
	cout << endl;
#endif
     
     	int current;
     	ZZ min;


	*(Parameters->Current_Simplicial_Cones_Total) = *(Parameters->Current_Simplicial_Cones_Total) + m;

	
	
	if(*(Parameters->Current_Simplicial_Cones_Total) > *(Parameters->Max_Simplicial_Cones_Total))
		*(Parameters->Max_Simplicial_Cones_Total) = *(Parameters->Current_Simplicial_Cones_Total);
	
	
     	//cout << "barvinok_DFS: max = " << max << endl;
     	for(int i = 0; i < m; i++)
     	{
	
		min = max + 1;
	     
		for(int j = 0; j < m; j++)
		{
			if(Dets[j] < min && Dets[j] != -1)
			{
				current = j;
				min = Dets[j];
			}
		}

	
	//cout << "barvinok_DFS: current = " << current << "   Dets[current] = " << Dets[current] << endl;	
   	if(barvinok_DFS(cones1[current], Parameters) == -1)
		result = -1;
	
      	cones1[current]->generator.kill();
	Dets[current] = -1;
	delete cones1[current];

	*(Parameters->Current_Simplicial_Cones_Total) = *(Parameters->Current_Simplicial_Cones_Total) - 1;
	//cout << " - " ;
	
     	}
  
    	//cout << "barvinok_DFS: Done with DFS on " << m << " items." << endl; 
      
     	Z.kill();
     
     	//delete [] cones1;

     	return result;
   
 }








