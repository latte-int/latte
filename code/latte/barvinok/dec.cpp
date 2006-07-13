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
#include "convert.h"
#include "dual.h"
#include "genFunction/piped.h"
#include "Residue.h"

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
	       char *File_Name, int max_determinant,
	       bool dualize,
	       BarvinokParameters::DecompositionType decomposition)
{
  int numOfConesDecomposed,numOfAllCones;
  listCone *tmp;

  Collecting_Single_Cone_Parameters parameters;
  parameters.Flags = Flags;
  parameters.Number_of_Variables = numOfVars;
  parameters.max_determinant = max_determinant;
  parameters.File_Name = File_Name;
  parameters.decomposition = decomposition;

  if (dualize) {
    parameters.dualize_time.start();
    cones = dualizeCones(cones, numOfVars);
    parameters.dualize_time.stop();
    cout << parameters.dualize_time;
  }
    
  cout << "Decomposing all cones.\n";
  numOfConesDecomposed=0;
  numOfAllCones=lengthListCone(cones);

  parameters.Cone_Index = 0;
  tmp=cones;
  while (tmp) {
    //parameters.Cone = tmp;
    int result = barvinokDecomposition_Single(tmp,
					      &parameters);
    assert(result >= 0);
    
    tmp=tmp->rest;
    numOfConesDecomposed++;
    if (numOfConesDecomposed==50*(numOfConesDecomposed/50)) {
      cout << numOfConesDecomposed << " / " << numOfAllCones << " done.\n";
    }
    parameters.Cone_Index++;
  }

  cout << "All cones have been decomposed.\n";
  cout << lengthListCone(parameters.Decomposed_Cones) << " cones in total.\n";
  
  return parameters.Decomposed_Cones;
}

/* ----------------------------------------------------------------- */



vec_ZZ
guess_generic_vector(int numOfVars)
{
  vec_ZZ result;
  result.SetLength(numOfVars);
  int i;
  for (i = 0;  i < numOfVars; i++)
    result[i] = (rand () % MODULUS) * ((rand() % 2) * 2 - 1);
  return result;
}

void
Generic_Vector_Single_Cone_Parameters::InitializeComputation()
{
  generic_vector = guess_generic_vector(Number_of_Variables);
}

void
barvinokDecomposition_List(listCone *cones,
			   Generic_Vector_Single_Cone_Parameters &Parameters)
{
  do {
    Parameters.InitializeComputation();
    try {
      listCone *cone;
      for (cone = cones; cone != NULL; cone = cone->rest) {
	int status;
	status = barvinokDecomposition_Single(cone, &Parameters);
	if (status < 0) {
	  static NotGenericException not_generic;
	  throw not_generic;  // FIXME: Later replace this return
			      // value handling by exception.
	}
      }
      return;
    }
    catch (NotGenericException) {
      cout << "Generic vector chosen unsuccessfully, trying again." << endl;
    };
  } while (1);
}



void
Standard_Single_Cone_Parameters::InitializeComputation()
{
  Generic_Vector_Single_Cone_Parameters::InitializeComputation();
  int i;
  for (i = 0; i <= Degree_of_Taylor_Expansion; i++)
    Taylor_Expansion_Result[i] = 0;
  Total_Lattice_Points = 0;
  Total_Uni_Cones = 0;
  Cone_Index = 0; 
}

int
Standard_Single_Cone_Parameters::ConsumeCone(listCone *Cone)
{
  //cout << "barvinok_DFS: Calculating points in Parallelepiped" << endl;
  if ( (Flags & DUAL_APPROACH) == 0)
    computePointsInParallelepiped(Cone, Number_of_Variables);
  //printListCone(Cone, Number_of_Variables);
		
  if (Flags & DUAL_APPROACH)	
    {
      //cout << "barvinok_DFS: Calling ResidueFunction_Single_Cone" << endl;
      return ResidueFunction_Single_Cone (Cone, this);
    }	
  else
    {
      return Residue_Single_Cone (Cone, Number_of_Variables, generic_vector, &Total_Lattice_Points, &Ten_Power);
    }	
}

void
decomposeAndComputeResidue(listCone *cones, int degree, bool dualize, 
			   Standard_Single_Cone_Parameters &param)
{
  int numOfConesDecomposed,numOfAllCones,numOfRays;
  mat_ZZ mat;
  listCone *tmp;
  int	Success = 0;

  if (dualize) {
    param.dualize_time.start();
    cones = dualizeCones(cones, param.Number_of_Variables);
    param.dualize_time.stop();
    cout << param.dualize_time;
  }
	
  cout << "decomposeCones_Single: Decomposing all cones. (Memory Save on)\n";
  numOfAllCones=lengthListCone(cones);
  cout << numOfAllCones << " cones total to be done!";	
	
  //Set Ten_Power to 100 billion
  // FIXME: What is magic about this number? --mkoeppe, Sat Mar  4 21:21:45 PST 2006
  param.Ten_Power = 1;
  for (int i = 0; i < Exponent_Ten_Power; i++)	
    param.Ten_Power *= 10;

  cout << "decomposeCones_Single: degree = " << degree << endl;
  param.Taylor_Expansion_Result = new ZZ [degree + 1 ];
	
  param.Degree_of_Taylor_Expansion = degree;
  param.Controller = new Node_Controller(param.Number_of_Variables + 1, degree);			

  cout << "Number of cones: " << numOfAllCones << endl;

  param.decompose_time.start();
  barvinokDecomposition_List(cones, param);
  param.decompose_time.stop();

  cout << endl << "Total Unimodular Cones: " << param.Total_Uni_Cones << endl;
  ofstream UniOut("numOfUnimodularCones");
  UniOut << param.Total_Uni_Cones << endl;
  cout << "Maximum number of simplicial cones in memory at once: " << param.Max_Simplicial_Cones_Total << endl;

	
  if ( param.Flags & DUAL_APPROACH)
    {
      cout << "Memory Save Mode: Taylor Expansion:" << endl;
      if(degree > 1){		
	for (int i = 0; i<= degree; i++)
	  {
	    cout <<	( param.Taylor_Expansion_Result[i] + param.Ten_Power/2)/(param.Ten_Power) ;
	    if (i != 0)
	      {
		cout << "t^" << i;
	      }
	    cout << endl;
	  }
      }
      else if(degree == 1){
	Integer numOfLatticePoints
	  = ( param.Taylor_Expansion_Result[1] + param.Ten_Power/2)/param.Ten_Power;
	cout << "\n****  Total number of lattice points is: " << numOfLatticePoints << "  ****" << endl << endl;
	ofstream out("numOfLatticePoints");
	out << numOfLatticePoints;
      }
    }
  else 
    {
#if 0
      cout << "Result: " << param.Total_Lattice_Points << " / "
	   << param.Ten_Power << endl;
#endif
      cout << "\n*****  Total number of lattice points: ";
      ofstream out("numOfLatticePoints");
      param.Total_Lattice_Points = abs( param.Total_Lattice_Points);
		
      cout <<	( param.Total_Lattice_Points + param.Ten_Power/2)/param.Ten_Power;
      out <<	( param.Total_Lattice_Points + param.Ten_Power/2)/param.Ten_Power << endl;
      cout << "  ****" << endl << endl;
    }
	
  //cout << lengthListCone(newCones->rest) << " cones in total.\n";

  delete param.Controller;
  delete [] param.Taylor_Expansion_Result;
	
}

void decomposeCones_Single (listCone *cones, int numOfVars, int degree,
			    unsigned int flags, char *File_Name, int max_determinant,
			    bool dualize,
			    BarvinokParameters::DecompositionType decomposition) 
{
  Standard_Single_Cone_Parameters *Barvinok_Parameters
    = new Standard_Single_Cone_Parameters;
	
  Barvinok_Parameters->Flags = flags;
  Barvinok_Parameters->Number_of_Variables = numOfVars;
  Barvinok_Parameters->max_determinant = max_determinant;
  Barvinok_Parameters->File_Name = File_Name;
  Barvinok_Parameters->decomposition = decomposition;

  decomposeAndComputeResidue(cones, degree, dualize, *Barvinok_Parameters);

  delete Barvinok_Parameters;
}
/* ----------------------------------------------------------------- */

