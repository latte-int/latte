/* ReadPolyhedron.cpp -- Handle command-line args to read a polyhedron

   Copyright 2007 Matthias Koeppe
   Derived from count.cpp, which is:
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

#include <cstring>

#include "CheckEmpty.h"
#include "ReadingFile.h"
#include "ReadLatteStyle.h"
#include "ReadPolyhedron.h"
#include "vertices/cdd.h"
#include "convert.h"
#include "preprocess.h"
#include "ramon.h"

ReadPolyhedronData::ReadPolyhedronData()
{
  strcpy(Vrepresentation,"no");
  strcpy(interior,"no");
  strcpy(dilation,"no");
  strcpy(dualApproach,"no");
  strcpy(nonneg, "no");
  strcpy(cddstyle, "no");
  strcpy(equationsPresent,"no");
  strcpy(Singlecone,"no");
  strcpy(grobner,"no");
  strcpy(maximum,"no");
  strcpy(minimize,"no");
  strcpy(assumeUnimodularCones,"no");
  strcpy(taylor,"no");
  strcpy(rationalCone,"no");
  strcpy(inthull, "no");
  strcpy(LRS, "no");
  strcpy(Memory_Save, "yes");

  expect_dilation_factor = false;
  dilation_const = 1;
  expect_filename = true;
  degree = 1;

  matrix = NULL;
  templistVec = NULL;
}

bool ReadPolyhedronData::parse_option(const char *arg)
{
  if (strncmp(arg,"vrep",3)==0) strcpy(Vrepresentation,"yes"); 
  else if(strncmp(arg,"int",3)==0) strcpy(interior,"yes");
  else if(strncmp(arg,"homog",3)==0) strcpy(dualApproach,"yes");
  else if(strncmp(arg,"equ",3)==0) strcpy(equationsPresent,"yes");
  else if(strncmp(arg,"+", 1) ==0) strcpy(nonneg,"yes");
  else if(strncmp(arg,"cdd",3)==0) strcpy (cddstyle, "yes");
  else if(strncmp(arg,"dil",3)==0) {
    strcpy (dilation, "yes");
    expect_dilation_factor = true;
  }
  else if(strncmp(arg,"--", 2)!=0) {
    // Regular argument, see if we expect one
    if (expect_dilation_factor) {
      dilation_const = atoi(arg);
      expect_dilation_factor = false;
    }
    else if (expect_filename) {
      filename = arg;
      expect_filename = false;
    }
    else
      return false;
  }
  else return false;
  return true;
}

Polyhedron *ReadPolyhedronData::read_polyhedron(BarvinokParameters *params)
{
  Polyhedron *Poly = NULL;

  if (expect_filename) {
    cerr << "The input file name is missing." << endl;
    exit(2);
  }

  listVector *equations = NULL, *inequalities = NULL;
  
  if (Vrepresentation[0] == 'y') {
    /* The polyhedron is given by its V-representation in a
       LattE-style input format. */
    if (cddstyle[0] == 'y') {
      cerr << "The command-line keyword `vrep' denotes the use of a LattE-style " << endl
	   << "input format.  It is not compatible with using a CDD-style input format." << endl;
      exit(2);
    }
    if(equationsPresent[0] == 'y') {
      cerr<<"You cannot use vrep and equ at the same time." << endl;
      exit(4);
    }
    if (dilation_const != 1) {
      cerr << "Dilation unimplemented for `vrep' input" << endl;
      exit(1);
    }
    Poly = ReadLatteStyleVrep(filename.c_str(), dualApproach[0] == 'y');
  }
  else {
    /* Not VREP. */
    Poly = new Polyhedron;
    int numOfVars;
  
    if((dualApproach[0] == 'y') && (nonneg[0] == 'y')&&(equationsPresent[0] == 'n')){
      cerr<<"You cannot use + and dua at the same time." << endl;
      exit(2);
    }
  
    if((Memory_Save[0] == 'y') && (inthull[0] == 'y')){
      cerr<<"You cannot use int and memsave at the same time." << endl;
      exit(3);
    }

  
    /* Check input file. */
    if(Vrepresentation[0] == 'n'){
      if((cddstyle[0] == 'n') && (maximum[0] == 'n')&& (minimize[0] == 'n')){
	CheckInputFile(filename.c_str());
	CheckLength(filename.c_str(), equationsPresent);
      }
      if(minimize[0] == 'y')  strcpy (maximum, "yes");   
 
      if((cddstyle[0] == 'n') && (maximum[0] == 'y')){
	CheckInputFile(filename.c_str());
	CheckLength(filename.c_str(),equationsPresent);
      }

      if(cddstyle[0] == 'y')
	{ CheckInputFileCDDRep(filename.c_str());
	  CheckInputFileCDDRep1(filename.c_str());
	  CheckInputFileCDDRep3(filename.c_str());
	  CheckInputFileCDDRep4(filename.c_str());
	}
    }else CheckInputFileVrep(filename.c_str());
    if (cddstyle[0] == 'y') {
      cerr << "Warning: Not performing check for empty polytope, "
	   << "because it is unimplemented for the CDD-style input format. "
	   << endl;
    }
    else {
      CheckEmpty(filename.c_str());
    }
    //vec_ZZ cost;
    /* Read problem data. */
    params->read_time.start();
  
    if((cddstyle[0] == 'n') && (Vrepresentation[0] == 'n'))
      CheckRed(filename, equationsPresent, maximum, nonneg, interior, dilation, dilation_const); 

    dilation_const = 1;
    if((cddstyle[0] == 'n'))
      readLatteProblem(filename.c_str(),&equations,&inequalities,equationsPresent,
		       &numOfVars, nonneg, dualApproach, grobner, Vrepresentation);
    if(cddstyle[0] == 'y'){
      int tmpoutput;
      CDDstylereadLatteProblem(filename.c_str(),&equations,&inequalities,equationsPresent,
			       &numOfVars, nonneg, dualApproach, taylor, degree,
			       rationalCone, tmpoutput, Memory_Save,
			       assumeUnimodularCones, inthull, grobner);
    }
  
    if((dualApproach[0] == 'y') && (nonneg[0] == 'y')&&(equationsPresent[0] == 'n')){
      cerr<<"You cannot use + and dua at the same time." << endl;
      exit(2);
    }
  
    if((Memory_Save[0] == 'y') && (inthull[0] == 'y')){
      cerr<<"You cannot use int and memsave at the same time." << endl;
      exit(3);
    }
  
    numOfVars--;

    mat_ZZ ProjU, ProjU2;
    ProjU.SetDims(numOfVars, numOfVars);
    ProjU2.SetDims(numOfVars, numOfVars);
    oldnumofvars = numOfVars;
    {
      listVector *matrixTmp;
      if (equationsPresent[0]=='y') {
	{
	  vec_ZZ *generators = NULL;
	  matrixTmp=preprocessProblem(equations,inequalities,&generators,&numOfVars, cost, ProjU, interior, dilation_const);
	  if (generators) delete[] generators;
	}
	freeListVector(equations);
	freeListVector(inequalities);
	ProjU2 = transpose(ProjU);
	bb = ProjU2[0];
	mat_ZZ AAA;
	AAA.SetDims(ProjU2.NumRows() - 1, ProjU2.NumCols());
	int i;
	for(i = 1; i <= numOfVars; i++){
	  AAA[i - 1] = ProjU2[i];
	}
	AA = transpose(AAA);
	// cout << ProjU << determinant(transpose(AAA)*AAA) <<  endl;
	templistVec = transformArrayBigVectorToListVector(ProjU, ProjU.NumCols(), ProjU.NumRows()); 
      } else {
	dilateListVector(inequalities, numOfVars, dilation_const);
	matrixTmp=inequalities;
      }
      if((dualApproach[0] == 'y')&&(equationsPresent[0]=='y')){
	matrix = TransformToDualCone(matrixTmp,numOfVars);
	freeListVector(matrixTmp);
      }
      else {
	matrix = matrixTmp;
      }
    }
    /* Now matrix contains the new inequalities. */
    params->read_time.stop();
    cout << params->read_time;
  
    //   cout << "Project down cost function: " << cost << endl;
    vec_RR Rat_solution, tmp_den, tmp_num;
    mat_RR ProjU_RR;
    ProjU_RR.SetDims(ProjU.NumRows(), ProjU.NumCols());
    int i;
    for(i = 0; i < ProjU.NumRows(); i++)
      for(int j = 0; j < ProjU.NumCols(); j++) conv(ProjU_RR[i][j], ProjU[i][j]);
    //cout << ProjU << ProjU_RR << endl;
    Rat_solution.SetLength(numOfVars);
    tmp_den.SetLength(numOfVars);
    tmp_num.SetLength(numOfVars);
    /* Compute vertices and edges. */
    rationalVector* LP_vertex;
    
    params->vertices_time.start();
    if ((dualApproach[0]=='n') && (Vrepresentation[0] == 'n')) {
      listCone *tmpcones;
      if(LRS[0] == 'n')
	tmpcones=computeVertexCones(filename.c_str(),matrix,numOfVars);
      else
	tmpcones=computeVertexConesViaLrs(filename.c_str(),matrix,numOfVars);
      if(maximum[0] == 'y'){ 
	LP_vertex = LP(matrix, cost, numOfVars);
	vec_RR Rat_cost;  Rat_cost.SetLength(numOfVars);
	for (i = 0; i < numOfVars; i++){
	  conv(tmp_num[i], LP_vertex->numerators()[i]);
	  conv(tmp_den[i], LP_vertex->denominators()[i]);
	  Rat_solution[i] = tmp_num[i]/tmp_den[i];
	  conv(Rat_cost[i], cost[i]);
	}
	if(Singlecone[0] == 'y')
	  Poly->cones = CopyListCones(tmpcones, numOfVars, LP_vertex);
	else Poly->cones = tmpcones;
	if(lengthListCone(Poly->cones) == 1) 
	  cout <<"\nWe found a single vertex cone for IP.\n" << endl;
	cout <<"A vertex which we found via LP is: " << ProjectingUpRR(ProjU_RR, Rat_solution, numOfVars) << endl;
	//printRationalVector(LP_vertex, numOfVars);
	RR LP_OPT;
	vec_RR holdcost_RR;
	holdcost_RR.SetLength(cost.length());
	for(i = 0; i < cost.length(); i++) conv(holdcost_RR[i], cost[i]);
	if(minimize[0] == 'y') holdcost_RR = - holdcost_RR;
	LP_OPT = Rat_cost*Rat_solution; //cout << cost << endl;
	cout << "The LP optimal value is: " << holdcost_RR*ProjectingUpRR(ProjU_RR, Rat_solution, numOfVars) << endl;
      }
      else {
	Poly->cones = tmpcones;
	cout << "The polytope has " << lengthListCone(Poly->cones) << " vertices." << endl;
	//system_with_error_check("rm -f numOfLatticePoints");
      }
    } 

    params->vertices_time.stop();
    cout << params->vertices_time;

    /* Compute triangulation or decomposition of each vertex cone. */

    if (dualApproach[0]=='y') {
      listVector *rays = NULL, *endRays, *tmpRays;
      Poly->cones=createListCone();
      Poly->cones->vertex = new Vertex(createRationalVector(numOfVars));
      rays=createListVector(createVector(numOfVars));
      endRays=rays;
      tmpRays=matrix;
      while (tmpRays) {
	vec_ZZ v=createVector(numOfVars);
	for (i=0; i<numOfVars; i++) v[i]=-(tmpRays->first)[i+1];
	endRays->rest=createListVector(v);
	endRays=endRays->rest;
	tmpRays=tmpRays->rest;
      }
      Poly->cones->rays = rays->rest;
      delete rays; // deletes dummy head
      Poly->dualized = true;

      //     cout << "Homogenization: " << endl;
      //     printListCone(Poly->cones, numOfVars);
    }

    Poly->numOfVars = numOfVars;
  } /* Not VREP */

  return Poly;
}
