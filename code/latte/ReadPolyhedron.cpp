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
#include "ReadSubcones.h"
#ifdef HAVE_FORTYTWO_LIB
#  include "VertexConesWith4ti2.h"
#endif
#include "print.h"

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
  strcpy(Memory_Save, "yes");

  vertexcones = VertexConesWithCdd;
  redundancycheck = RedundancyCheckWithCdd;
  
  expect_dilation_factor = false;
  dilation_const = 1;
  expect_filename = true;
  degree = 1;

  input_homog_cone = false;
  input_dualized = false;
  have_subcones = false;
  input_listcone_format = false;

  matrix = NULL;
  templistVec = NULL;
}

void ReadPolyhedronData::show_options(ostream &stream)
{
  stream << "Standard input specifications:" << endl
         << "  FILENAME                                 Inequalities in LattE format" << endl
	 << "    dil DILATION-FACTOR                    - Dilate by DILATION-FACTOR" << endl
	 << "    +                                      - Add non-negativity constraints" << endl
         << "    int                                    - Handle the interior of the polyhedron" << endl
	 << "  vrep FILENAME                            Vertices in LattE format" << endl
	 << "  cdd FILENAME.{ext,ine}                   Inequalities or vertices in CDD format" << endl
	 << "Intermediate input specifications:" << endl
	 << "  --input-primal-homog-cone=CONE.ext       The homogenized polyhedron given by a " << endl
	 << "                                           full-dimensional cone in CDD format" << endl
	 << "  --input-dual-homog-cone=CONE.ext         The dual of the homogenized polyhedron given by a " << endl
         << "                                           full-dimensional cone in CDD format" << endl
	 << "  --subcones=FILENAME                      Use a subdivision of the above specified" << endl
         << "                                           cone (up to lower-dimensional cones), given by " << endl
         << "                                           ray indicator vectors" << endl
	 << "  --input-primal-homog-cones=CONES         The homogenized polyhedron given by a " << endl
	 << "                                           union of cones (up to lower-dimensional cones) " << endl
	 << "                                           in LattE's internal format" << endl
	 << "  --input-dual-homog-cones=CONES           The dual of the homogenized polyhedron given by a " << endl
         << "                                           union of cones (up to lower-dimensional cones) " << endl
	 << "                                           in LattE's internal format" << endl
         << "Input handling options:" << endl
	 << "  --compute-vertex-cones={cdd,lrs,4ti2}    Use this method for computing vertex cones" << endl
	 << "  --redundancy-check={none,cdd}            Use this method for computing an irredundant " << endl
	 << "                                           representation" << endl
    ;
}

bool ReadPolyhedronData::parse_option(const char *arg)
{
  /* Parse traditional LattE options. */
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
  else if (strncmp(arg, "lrs", 3)==0)
    vertexcones = ReadPolyhedronData::VertexConesWithLrs;
  /* Parse new options. */
  else if (strncmp(arg, "--input-primal-homog-cone=", 26)==0) {
    filename = arg + 26;
    expect_filename = false;
    input_homog_cone = true;
    input_dualized = false;
    strcpy(dualApproach,"yes");
  }
  else if (strncmp(arg, "--input-dual-homog-cone=", 24)==0) {
    filename = arg + 24;
    expect_filename = false;
    input_homog_cone = true;
    input_dualized = true;
    strcpy(dualApproach,"yes");
  }
  else if (strncmp(arg, "--subcones=", 11) == 0) {
    subcones_filename = string(arg + 11);
    have_subcones = true;
  }
  else if (strncmp(arg, "--input-primal-homog-cones=", 27)==0) {
    filename = arg + 27;
    expect_filename = false;
    input_homog_cone = true;
    input_dualized = false;
    input_listcone_format = true;
    strcpy(dualApproach,"yes");
  }
  else if (strncmp(arg, "--input-dual-homog-cones=", 25)==0) {
    filename = arg + 25;
    expect_filename = false;
    input_homog_cone = true;
    input_dualized = true;
    input_listcone_format = true;
    strcpy(dualApproach,"yes");
  }
  else if (strncmp(arg, "--compute-vertex-cones=", 23) == 0) {
    if (strcmp(arg + 23, "cdd") == 0) 
      vertexcones = VertexConesWithCdd;
    else if (strcmp(arg + 23, "lrs") == 0)
      vertexcones = VertexConesWithLrs;
    else if (strcmp(arg + 23, "4ti2") == 0)
      vertexcones = VertexConesWith4ti2;
    else {
      cerr << "Unknown vertex cone method: " << arg + 23 << endl;
      exit(1);
    }
  }
  else if (strncmp(arg, "--redundancy-check=", 19) == 0) {
    if (strcmp(arg + 19, "none") == 0)
      redundancycheck = NoRedundancyCheck;
    else if (strcmp(arg + 19, "cdd") == 0)
      redundancycheck = RedundancyCheckWithCdd;
    else {
      cerr << "Unknown redundancy check method: " << arg + 19 << endl;
      exit(1);
    }
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

static listCone *
read_cone_cdd_format(const string &filename)
{
  FILE *in = fopen(filename.c_str(), "r");
  if (in == NULL) {
    cerr << "Unable to open CDD-style input file " << filename << endl;
    exit(1);
  }
  dd_MatrixPtr M;
  dd_ErrorType err=dd_NoError;
  M = dd_PolyFile2Matrix(in, &err);
  if (err!=dd_NoError) {
    cerr << "Parse error in CDD-style input file " << filename << endl;
    exit(1);
  }
  listCone *cone = cddlib_matrix_to_cone(M);
  dd_FreeMatrix(M);
  return cone;
}

Polyhedron *
ReadPolyhedronData::read_polyhedron(BarvinokParameters *params)
{
  if (expect_filename) {
    cerr << "The input file name is missing." << endl;
    exit(2);
  }

  if (input_homog_cone)
    return read_polyhedron_from_homog_cone_input(params);
  else
    return read_polyhedron_hairy(params);
}

Polyhedron *
ReadPolyhedronData::read_polyhedron_from_homog_cone_input(BarvinokParameters *params)
{
  /* We are already given a full-dimensional, homogenized cone
     or a list of those. */
  ConeProducer *producer = NULL;
  if (input_listcone_format) {
    if (have_subcones) {
      listCone *cones = readListConeFromFile(filename.c_str());
      if (lengthListCone(cones) != 1) {
	cerr << "A subcones file can only be given for a single-cone file." << endl;
	exit(1);
      }
      producer = new SubconeReadingConeProducer(cones, subcones_filename);
    }
    else {
      producer = new ListConeReadingConeProducer(filename);
    }
  }
  else {
    listCone *cone = read_cone_cdd_format(filename);
    if (have_subcones) {
      // Also a subcones file given.
      producer = new SubconeReadingConeProducer(cone, subcones_filename);
    }
    else {
      producer = new SingletonConeProducer(copyCone(cone));
    }
  }
  /* Use the producer to create the polyhedron. */
  CollectingConeConsumer ccc;
  producer->Produce(ccc);
  delete producer;
  Polyhedron *Poly = new Polyhedron;
  Poly->cones = ccc.Collected_Cones;
  int numOfVars;
  if (Poly->cones == NULL || Poly->cones->rays == NULL)
    numOfVars = 0;
  else
    numOfVars = Poly->cones->rays->first.length();
  Poly->numOfVars = numOfVars;
  Poly->homogenized = true;
  Poly->dualized = input_dualized;
  return Poly;
}

Polyhedron *
ReadPolyhedronData::read_polyhedron_hairy(BarvinokParameters *params)
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
  
    if((cddstyle[0] == 'n') && (Vrepresentation[0] == 'n')) {
      switch (redundancycheck) {
      case NoRedundancyCheck:
	break;
      case RedundancyCheckWithCdd:
	/* this changes FILENAME: */
	CheckRed(filename, equationsPresent, maximum, nonneg, interior, dilation, dilation_const);
	break;
      case RedundancyCheckWithCddlib:
	cerr << "RedundancyCheckWithCddlib unimplemented." << endl;
	exit(1);
	break;
      default:
	cerr << "Unknown redundancy check" << endl;
	abort();
      }
    }

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

      switch (vertexcones) {
      case VertexConesWithCdd:
	tmpcones=computeVertexCones(filename.c_str(),matrix,numOfVars);
	break;
      case VertexConesWithLrs:
	tmpcones=computeVertexConesViaLrs(filename.c_str(),matrix,numOfVars);
	break;
      case VertexConesWith4ti2:
#ifdef HAVE_FORTYTWO_LIB
	tmpcones=computeVertexConesWith4ti2(matrix, numOfVars);
#else
	cerr << "VertexConesWith4ti2 not compiled in, sorry" << endl;
	exit(1);
#endif
	break;
      default:
	cerr << "Bad VertexConesType" << endl;
	abort();
      };
	  
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

  if (dualApproach[0] == 'y')
    Poly->homogenized = true;
  else
    Poly->homogenized = false;
  
  return Poly;
}
