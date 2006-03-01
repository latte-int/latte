/* ----------------------------------------------------------------- */
/*                                                                   */
/* LattE (Lattice Point Enumeration)                                 */
/*                                                                   */
/* Calling Rudy's Decomposition Procedure                            */
/*                                                                   */
/* Author     : Raymond Hemmecke, Ruriko Yoshida                     */
/*                                                                   */
/* Created    : 13-SEP-02                                            */
/* Last Update: 03-Feb-03 by Rudy                                    */
/*                                                                   */
/* ----------------------------------------------------------------- */

#include <cassert>
#include <fstream>
#include <stdlib.h>
#include <time.h>

#include "../myheader.h"
#include "../cone.h"
#include "../print.h"
#include "../ramon.h"
#include "ConeDecom.h"
#include "../RudyResNTL.h"
#include "../PolyTree.h"
#include "../flags.h"
#include "dual.h"
#include "genFunction/piped.h"

#define MODULUS 1000000000
#define Exponent_Ten_Power 10


/* ----------------------------------------------------------------- */
listCone* readListCone(rationalVector *vertex, int numOfVars) {
  int i,j,k,numOfCones,coeff;
  vec_ZZ v;
  listVector *tmp, *endtmp;
  listCone *cones, *endCones;
  char fileName[127];

  cones=createListCone();
  endCones=cones;

  strcpy(fileName,"latte_dec.latte");

  ifstream in(fileName);
  if(!in){
    cerr << "Cannot open input file in readListCone." << endl;
    exit(1);
  }

  in >> numOfCones;

  for (i=0; i<numOfCones; i++) {
    in >> coeff;
    v=createVector(numOfVars);
    tmp=createListVector(v);
    endtmp=tmp;
    for (j=0; j<numOfVars; j++) {
      v=createVector(numOfVars);
      for (k=0; k<numOfVars; k++) {
	in >> v[k];
      }
      endtmp->rest=createListVector(v);
      endtmp=endtmp->rest;
    }
    endCones->rest=createListCone();
    endCones=endCones->rest;
    if (coeff==1) endCones->coefficient=1;
    else endCones->coefficient=-1;
    endCones->vertex=vertex;
    endCones->rays=tmp->rest;
  }

  in.close();

  return (cones->rest);
}
/* ----------------------------------------------------------------- */
mat_ZZ createConeDecMatrix(listCone *cones, int numOfRays, int numOfVars) {
  int i;
  mat_ZZ mat;
  listVector *tmp;

  mat.SetDims(numOfRays, numOfVars);

  tmp=cones->rays;
  for (i=0; i<numOfRays; i++) {
    mat[i]=copyVector(tmp->first,numOfVars);
    tmp=tmp->rest;
  }
  //removeListVector(cones->rays);
  return (mat);
}
/* ----------------------------------------------------------------- */

Collecting_Single_Cone_Parameters::Collecting_Single_Cone_Parameters()
  : Decomposed_Cones(NULL)
{
}

int Collecting_Single_Cone_Parameters::ConsumeCone(listCone *cone)
{
  assert(cone->rest == NULL);
  cone->rest = Decomposed_Cones;
  Decomposed_Cones = cone;
  return 1; // means "success, please continue"
}

listCone*
decomposeCones(listCone *cones, int numOfVars, unsigned int Flags,
	       char *File_Name, int max_determinant)
{
  int numOfConesDecomposed,numOfAllCones,numOfRays, numOfUniCones = 0;
  mat_ZZ mat;
  listCone *tmp;

  Collecting_Single_Cone_Parameters parameters;
  parameters.Flags = Flags;
  parameters.Number_of_Variables = numOfVars;
  parameters.max_determinant = max_determinant;

  cout << "Decomposing all cones.\n";
  numOfConesDecomposed=0;
  numOfAllCones=lengthListCone(cones);

  int Cone_Index = 0;
  
  tmp=cones;
  while (tmp) {
    numOfRays=lengthListVector(tmp->rays);
    mat=createConeDecMatrix(tmp,numOfRays,numOfVars);

    parameters.Cone = tmp;

    
    int result = barvinokDecomposition_Single(mat, numOfRays,
					      numOfUniCones, tmp->vertex,
					      &parameters, File_Name,
					      Cone_Index);
    assert(result >= 0);
    
    tmp=tmp->rest;
    numOfConesDecomposed++;
    if (numOfConesDecomposed==50*(numOfConesDecomposed/50)) {
      cout << numOfConesDecomposed << " / " << numOfAllCones << " done.\n";
    }
    Cone_Index++;
  }

  cout << "All cones have been decomposed.\n";
  cout << lengthListCone(parameters.Decomposed_Cones) << " cones in total.\n";
  
  return parameters.Decomposed_Cones;
}

/* ----------------------------------------------------------------- */

int Standard_Single_Cone_Parameters::ConsumeCone(listCone *Cone)
{
  Cone = dualizeBackCones (Cone, this->Number_of_Variables );   	
  //cout << "barvinok_DFS: Calculating points in Parallelepiped" << endl;
  if ( (this->Flags & DUAL_APPROACH) == 0)
    computePointsInParallelepiped(Cone, this->Number_of_Variables);
		
  if (this->Flags & DUAL_APPROACH)	
    {
      //cout << "barvinok_DFS: Calling ResidueFunction_Single_Cone" << endl;
      return ResidueFunction_Single_Cone (Cone, (Standard_Single_Cone_Parameters *) this);
    }	
  else
    {
      return Residue_Single_Cone (Cone, this->Number_of_Variables, this->Random_Lambda, &this->Total_Lattice_Points, &this->Ten_Power);
    }	
}

void decomposeCones_Single (listCone *cones, int numOfVars, int degree, unsigned int flags, char *File_Name) 
{
	int numOfConesDecomposed,numOfAllCones,numOfRays, numOfUniCones = 0;
  	mat_ZZ mat;
  	listCone *tmp;
	/*  char command[127]; */
	Standard_Single_Cone_Parameters	*Barvinok_Parameters = new Standard_Single_Cone_Parameters;
	
	int	Success = 0;
	
  	cout << "decomposeCones_Single: Decomposing all cones. (Memory Save on)\n";
  	numOfAllCones=lengthListCone(cones);
	cout << numOfAllCones << " cones total to be done!";	

	
	//Set Ten_Power to 100 billion
	Barvinok_Parameters->Ten_Power = 1;
	for (int i = 0; i < Exponent_Ten_Power; i++)	
	  Barvinok_Parameters->Ten_Power *= 10;

	cout << "decomposeCones_Single: degree = " << degree << endl;
	Barvinok_Parameters->Taylor_Expansion_Result = new ZZ [degree + 1 ];
	
	
	Barvinok_Parameters->Degree_of_Taylor_Expansion = degree;
	Barvinok_Parameters->Flags = flags;
	Barvinok_Parameters->Number_of_Variables = numOfVars;
	Barvinok_Parameters->max_determinant = 1;
	Barvinok_Parameters->Random_Lambda = new ZZ [numOfVars];

	Barvinok_Parameters->Controller = new Node_Controller(numOfVars + 1, degree);			
	
	cout << "Number of cones: " << numOfAllCones << endl;
	
	while (Success == 0)
	{
		Success = 1;
  		numOfConesDecomposed = 0;
		tmp=cones;

		Barvinok_Parameters->Current_Simplicial_Cones_Total = numOfAllCones;
		Barvinok_Parameters->Max_Simplicial_Cones_Total = 0;

		//Compute Random lambda
		cout << "decomposeCone_Single: Random Lambda = ";
		for (int i = 0;  i < numOfVars; i++)
		{
			Barvinok_Parameters->Random_Lambda[i] = (rand () % MODULUS) * ((rand() % 2) * 2 - 1);
			cout << Barvinok_Parameters->Random_Lambda[i] << " ";
		}
		cout << endl;
			
		for (int i = 0; i <= degree; i++)
			Barvinok_Parameters->Taylor_Expansion_Result[i] = 0;
		
		Barvinok_Parameters->Total_Lattice_Points = 0;
		Barvinok_Parameters->Total_Uni_Cones = 0;
	
		int Cone_Index = 0;
		
  		while (tmp) 
		{
    			numOfRays=lengthListVector(tmp->rays);
    			mat=createConeDecMatrix(tmp,numOfRays,numOfVars);
			Barvinok_Parameters->Cone = tmp;
				
			// reminder, set vertex, pass vertex
			if(barvinokDecomposition_Single(mat,numOfRays, numOfUniCones, tmp->vertex,Barvinok_Parameters, File_Name, Cone_Index) == -1)
			{
				Success = 0;		
				break;
			}
    		
			
    			numOfConesDecomposed++;
    			if (numOfConesDecomposed==50*(numOfConesDecomposed/50)) 
			{	
      				cout << numOfConesDecomposed << " / " << numOfAllCones << " done.\n";
    			}
			
  			tmp = tmp->rest;
			
			Cone_Index++;
		}
		
		if (Success == 0)
		{
			cout << "Lambda Choosen unsuccessful, trying again." << endl;
		}
	}

	cout << endl << "Total Unimodular Cones: " << Barvinok_Parameters->Total_Uni_Cones << endl;
	ofstream UniOut("numOfUnimodularCones");
	UniOut << Barvinok_Parameters->Total_Uni_Cones << endl;
	cout << "Maximum number of simplicial cones in memory at once: " << Barvinok_Parameters->Max_Simplicial_Cones_Total << endl;

	
	if ( flags & DUAL_APPROACH)
	{
	  cout << "Memory Save Mode: Taylor Expansion:" << endl;
	  if(degree > 1){		
	    for (int i = 0; i<= degree; i++)
	      {
		cout <<	( Barvinok_Parameters->Taylor_Expansion_Result[i] + Barvinok_Parameters->Ten_Power/2)/(Barvinok_Parameters->Ten_Power) ;
		if (i != 0)
		  {
		    cout << "t^" << i;
		  }
		cout << endl;
	      }
	  }
	  else if(degree == 1){
	    cout << "\n****  Total number of lattice points is: " << ( Barvinok_Parameters->Taylor_Expansion_Result[1] + Barvinok_Parameters->Ten_Power/2)/Barvinok_Parameters->Ten_Power  << "  ****" << endl << endl;
	  }
	  
	}
	else 
	{
	  cout << "\n*****  Total number of lattice points: ";
	  ofstream out("numOfLatticePoints");
	  Barvinok_Parameters->Total_Lattice_Points = abs( Barvinok_Parameters->Total_Lattice_Points);
		
		cout <<	( Barvinok_Parameters->Total_Lattice_Points + Barvinok_Parameters->Ten_Power/2)/Barvinok_Parameters->Ten_Power;
		out <<	( Barvinok_Parameters->Total_Lattice_Points + Barvinok_Parameters->Ten_Power/2)/Barvinok_Parameters->Ten_Power << endl;
		cout << "  ****" << endl << endl;
	}
	
  	//cout << lengthListCone(newCones->rest) << " cones in total.\n";

	delete Barvinok_Parameters->Controller;
	delete [] Barvinok_Parameters->Random_Lambda;
	delete [] Barvinok_Parameters->Taylor_Expansion_Result;
	delete Barvinok_Parameters;
	
}
/* ----------------------------------------------------------------- */
