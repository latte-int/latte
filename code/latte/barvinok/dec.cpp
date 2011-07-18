/* dec.cpp -- Calling Rudy's Decomposition Procedure

   Copyright 2002, 2003 Raymond Hemmecke, Ruriko Yoshida
   Copyright 2006 Matthias Koeppe

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

#include <cassert>
#include <fstream>
#include <stdlib.h>
#include <time.h>

#include "cone.h"
#include "print.h"
#include "ramon.h"
#include "latte_ntl_integer.h"
#include "RudyResNTL.h"
#include "PolyTree.h"
#include "flags.h"
#include "convert.h"
#include "dual.h"
#include "genFunction/piped.h"
#include "Residue.h"

#define MODULUS 1000000000
#define Exponent_Ten_Power 10

Collecting_Single_Cone_Parameters::Collecting_Single_Cone_Parameters()
  : Decomposed_Cones(NULL)
{
}

Collecting_Single_Cone_Parameters::Collecting_Single_Cone_Parameters(const BarvinokParameters &params)
  : Single_Cone_Parameters(params), Decomposed_Cones(NULL)
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
decomposeCones(listCone *cones, bool dualize,
	       BarvinokParameters &param)
{
  int numOfConesDecomposed,numOfAllCones;
  listCone *tmp;

  Collecting_Single_Cone_Parameters parameters(param);

  if (dualize) {
    dualizeCones(cones, param.Number_of_Variables, &param);
  }
    
  cerr << "Decomposing all cones.\n";
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
      cerr << numOfConesDecomposed << " / " << numOfAllCones << " done.\n";
    }
    parameters.Cone_Index++;
  }

  cerr << "All cones have been decomposed.\n";
  cerr << lengthListCone(parameters.Decomposed_Cones) << " cones in total.\n";
  
  return parameters.Decomposed_Cones;
}

listCone*
decomposeCones(listCone *cones, int numOfVars, unsigned int Flags,
	       char *File_Name, int max_determinant,
	       bool dualize,
	       BarvinokParameters::DecompositionType decomposition,
	       bool debug_triangulation)
{
  Collecting_Single_Cone_Parameters parameters;
  parameters.Flags = Flags;
  parameters.Number_of_Variables = numOfVars;
  parameters.max_determinant = max_determinant;
  parameters.File_Name = File_Name;
  parameters.decomposition = decomposition;
  parameters.debug_triangulation = debug_triangulation;
  return decomposeCones(cones, dualize, parameters);
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
      int index;
      for (cone = cones, index = 0; cone != NULL; cone = cone->rest, index++) {
	int status;
	status = barvinokDecomposition_Single(cone, &Parameters);
	if (status < 0) {
	  static NotGenericException not_generic;
	  throw not_generic;  // FIXME: Later replace this return
			      // value handling by exception.
	}
	if (index%1 == 0)
	  cerr << index << " vertex cones done. " << endl;
      }
      return;
    }
    catch (NotGenericException) {
      cerr << "Generic vector chosen unsuccessfully, trying again." << endl;
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
  Max_Depth = 0;
  Current_Depth = 0;
}

int
Standard_Single_Cone_Parameters::ConsumeCone(listCone *Cone)
{
  //cerr << "barvinok_DFS: Calculating points in Parallelepiped" << endl;
  if ( (Flags & DUAL_APPROACH) == 0)
    computePointsInParallelepiped(Cone, Number_of_Variables, this);
  //printListCone(Cone, Number_of_Variables);
		
  if (Flags & DUAL_APPROACH)	
    {
      //cerr << "barvinok_DFS: Calling ResidueFunction_Single_Cone" << endl;
      return ResidueFunction_Single_Cone (Cone, this);
    }	
  else
    {
      return Residue_Single_Cone (Cone, Number_of_Variables, generic_vector, &Total_Lattice_Points, &Ten_Power);
    }	
}

vec_ZZ
decomposeAndComputeResidue(listCone *cones, int degree, bool dualize, 
			   Standard_Single_Cone_Parameters &param)
{
	vec_ZZ polynomialAnswer;
  int numOfAllCones;
  mat_ZZ mat;

  if (dualize) {
    dualizeCones(cones, param.Number_of_Variables, &param);
  }
	
  cerr << "decomposeCones_Single: Decomposing all cones. (Memory Save on)\n";
  numOfAllCones=lengthListCone(cones);
  cerr << numOfAllCones << " cones total to be done!";	
	
  //Set Ten_Power to 100 billion
  // FIXME: What is magic about this number? --mkoeppe, Sat Mar  4 21:21:45 PST 2006
  param.Ten_Power = 1;
  for (int i = 0; i < Exponent_Ten_Power; i++)	
    param.Ten_Power *= 10;

  cerr << "decomposeCones_Single: degree = " << degree << endl;
  param.Taylor_Expansion_Result = new ZZ [degree + 1 ];
  polynomialAnswer.SetLength(degree+1);
	
  param.Degree_of_Taylor_Expansion = degree;
  param.Controller = new Node_Controller(param.Number_of_Variables + 1, degree);			

  cerr << "Number of cones: " << numOfAllCones << endl;

  barvinokDecomposition_List(cones, param);

  cerr << endl << "Total Unimodular Cones: " << param.Total_Uni_Cones << endl;
  ofstream UniOut("numOfUnimodularCones");
  UniOut << param.Total_Uni_Cones << endl;
  cerr << "Maximum number of simplicial cones in memory at once: " << param.Max_Simplicial_Cones_Total << endl;

	
  if ( param.Flags & DUAL_APPROACH)
    {
      cerr << "Memory Save Mode: Taylor Expansion:" << endl;
      if(degree > 1){		
	for (int i = 0; i<= degree; i++)
	  {
		ZZ number;
		number = ( param.Taylor_Expansion_Result[i] + param.Ten_Power/2)/(param.Ten_Power) ;
	    cout <<	number;
	    polynomialAnswer[i] = number;

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
	param.deliver_number_of_lattice_points(numOfLatticePoints);
	polynomialAnswer.SetLength(2);
	polynomialAnswer[0] = 0;
	polynomialAnswer[1] = numOfLatticePoints;
      }
    }
  else {
    param.Total_Lattice_Points = abs(param.Total_Lattice_Points);
    ZZ number;
    number = (param.Total_Lattice_Points + param.Ten_Power/2)/param.Ten_Power;
    param.deliver_number_of_lattice_points(number);

    polynomialAnswer.SetLength(2);
    polynomialAnswer[0] = 0;
    polynomialAnswer[1] = number;

  }
	
  //cerr << lengthListCone(newCones->rest) << " cones in total.\n";

  delete param.Controller;
  delete [] param.Taylor_Expansion_Result;
  return polynomialAnswer;
	
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

