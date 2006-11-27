/* count.cpp -- Master program

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

#include <string.h>
#include <stdio.h>
#include <cassert>

#include "config.h"
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
#include "Residue.h"
#include "Grobner.h"
//  #include "jesus.h"
#include "preprocess.h"
#include "print.h"
#include "ramon.h"
#include "rational.h"
#include "timing.h"
#include "flags.h"
//#include "testing.h"
#include "IntegralHull.h"
#include "ReadingFile.h"
#include "binarySearchIP.h"
#include "CheckEmpty.h"
#include "ExponentialSubst.h"
#include "latte_random.h"
#include "Irrational.h"
#include "ExponentialEhrhart.h"
#ifdef HAVE_EXPERIMENTS
#include "ExponentialApprox.h"
#include "TrivialSubst.h"
#endif

#include "banner.h"
#include "convert.h"
#include "latte_system.h"
#include "Polyhedron.h"
#include "ReadLatteStyle.h"

static void usage(const char *progname)
{
  cerr << "usage: " << progname << " [OPTIONS...] " << "INPUTFILE" << endl;
}

/* ----------------------------------------------------------------- */
int main(int argc, char *argv[]) {
  float z;
  int i, degree = 1;
  unsigned int flags = 0, print_flag = 0, output_cone = 0;
  vec_ZZ dim, v, w;
  int oldnumofvars;
  vec_ZZ *generators;
  char fileName[127], invocation[127], equationsPresent[10],
    assumeUnimodularCones[127], dualApproach[127], taylor[127], printfile[127],
    rationalCone[127], nonneg[127], Memory_Save[127], Save_Tri[127],
    Load_Tri[127], Print[127], inthull[127], cddstyle[127], grobner[127],
    removeFiles[127], command[127], maximum[127],  Singlecone[127], LRS[127],
    Vrepresentation[127], dilation[127], minimize[127], binary[127], interior[127];
  bool approx;
  bool ehrhart_polynomial, ehrhart_series, ehrhart_taylor;
  double sampling_factor = 1.0;
  long int num_samples = -1;
  
  listVector *matrix = NULL, *equations = NULL, *inequalities = NULL, *rays = NULL, *endRays, *tmpRays, *matrixTmp;
  vec_ZZ cost;
  listVector *templistVec;
  listCone *tmpcones;
  mat_ZZ ProjU, ProjU2, AA;
  vec_ZZ bb;
  mat_ZZ AAA;
  RR LP_OPT;
  vec_ZZ holdCost;
  struct BarvinokParameters *params = new BarvinokParameters;

  latte_banner(cout);

  if (argc < 2) {
    usage(argv[0]);
    exit(1);
  }
  
  z=0;
  //setbuf(stdout,0);

  strcpy(invocation,"Invocation: ");
  strcat(invocation,argv[0]);
  strcat(invocation," ");
/*    strcat(invocation,argv[argc-1]); */
/*    strcat(invocation,"\n\n"); */
/*    printf(invocation); */

  strcpy(Vrepresentation,"no");
  strcpy(interior,"no");
  strcpy(dilation,"no");
  strcpy(binary,"no");
  strcpy(Singlecone,"no");
  strcpy(removeFiles,"yes");
  strcpy(grobner,"no");
  strcpy(maximum,"no");
  strcpy(minimize,"no");
  strcpy(dualApproach,"no");
  strcpy(equationsPresent,"no");
  strcpy(assumeUnimodularCones,"no");
  strcpy(printfile,"no");
  strcpy(taylor,"no");
  strcpy(rationalCone,"no");
  strcpy(nonneg, "no");
  strcpy(Memory_Save, "yes");
  strcpy(Save_Tri, "no");
  strcpy(Load_Tri, "no");
  strcpy(Print, "no");
  strcpy(inthull, "no");
  strcpy(cddstyle, "no");
  strcpy(LRS, "no");
  approx = false;
  ehrhart_polynomial = false;
  params->substitution = BarvinokParameters::PolynomialSubstitution;
  params->decomposition = BarvinokParameters::DualDecomposition;
  params->triangulation = BarvinokParameters::RegularTriangulationWithCdd;
  params->max_determinant = 1;

  int last_command_index = argc - 2;
  for (i=1; i<=last_command_index; i++) {
    strcat(invocation,argv[i]);
    strcat(invocation," ");
    if (strncmp(argv[i],"vrep",3)==0) strcpy(Vrepresentation,"yes"); 
    // else if(strncmp(argv[i],"bbs",3)==0) strcpy(binary,"yes");
    else if(strncmp(argv[i],"int",3)==0) strcpy(interior,"yes");
    //else if(strncmp(argv[i],"min",3)==0) strcpy(minimize,"yes");
    //else if(strncmp(argv[i],"gro",3)==0) strcpy(grobner,"yes");
    else if(strncmp(argv[i],"nodecom",3)==0
	    || strncmp(argv[i], "--nodecomposition", 5) == 0
	    || strncmp(argv[i], "--no-decomposition", 7) == 0)
      params->max_determinant = 0;
    else if(strncmp(argv[i],"homog",3)==0) {strcpy(dualApproach,"yes"); flags |= DUAL_APPROACH;}
    else if(strncmp(argv[i],"equ",3)==0) strcpy(equationsPresent,"yes");
    else if(strncmp(argv[i],"uni",3)==0) strcpy(assumeUnimodularCones,"yes");
    //else if(strncmp(argv[i],"simp",4)==0) {strcpy(printfile,"yes"); flags |= PRINT;}
    else if(strncmp(argv[i],"file",4)==0) strcpy(Memory_Save, "no");
    //else if(strncmp(argv[i],"single",6)==0) strcpy(Singlecone,"yes");
    //else if(strncmp(argv[i],"ehrhartsimp",3)==0) strcpy(rationalCone,"yes");
    else if(strncmp(argv[i],"+", 1) ==0) strcpy(nonneg,"yes");
    else if(strncmp(argv[i],"memsave",7)==0) strcpy (Memory_Save, "yes");
    else if(strncmp(argv[i],"printcones",3)==0) strcpy (Print, "yes");
    else if(strncmp(argv[i],"cdd",3)==0) strcpy (cddstyle, "yes");
    //else if(strncmp(argv[i],"hull",3)==0) strcpy (inthull, "yes");
    // else if(strncmp(argv[i],"max",3)==0) strcpy (maximum, "yes");
    else if(strncmp(argv[i],"lrs",3)==0) strcpy (LRS, "yes");
    else if(strncmp(argv[i],"dil",3)==0) {
      strcpy (dilation, "yes");
      last_command_index--;
    }
    else if(strncmp(argv[i],"rem",3)==0) strcpy (removeFiles, "no");
    else if(strncmp(argv[i],"rem",3)==0) strcpy (Memory_Save, "no");
    else if(strncmp(argv[i],"trisave",7)==0) {strcpy (Save_Tri, "yes"); flags |= SAVE;}
    else if(strncmp(argv[i],"triload",7)==0) {strcpy (Load_Tri, "yes"); flags |= LOAD;}
    else if (strncmp(argv[i], "--exponential", 5) == 0)
      params->substitution = BarvinokParameters::ExponentialSubstitution;
    else if (strncmp(argv[i], "--polynomial", 6) == 0)
      params->substitution = BarvinokParameters::PolynomialSubstitution;
    else if (strncmp(argv[i], "--maxdet=", 9) == 0)
      params->max_determinant = atoi(argv[i] + 9);
    else if (strncmp(argv[i], "--irrational-all-primal", 14) == 0
	     || strncmp(argv[i], "--all-primal", 5) == 0)
      params->decomposition = BarvinokParameters::IrrationalAllPrimalDecomposition;
    else if (strncmp(argv[i], "--irrational-primal", 5) == 0)
      params->decomposition = BarvinokParameters::IrrationalPrimalDecomposition;
    else if (strncmp(argv[i], "--dual", 6) == 0)
      params->decomposition = BarvinokParameters::DualDecomposition;
    else if (strncmp(argv[i], "--ehrhart-polynomial", 11) == 0)
      ehrhart_polynomial = true;
    else if (strncmp(argv[i], "--ehrhart-series", 11) == 0) {
      ehrhart_series = true;
      strcpy(dualApproach,"yes");
      flags |= DUAL_APPROACH;
      strcpy(printfile,"yes");
      flags |= PRINT;
    }
    else if (strncmp(argv[i], "--simplified-ehrhart-series", 14) == 0) {
      ehrhart_series = true;
      strcpy(dualApproach,"yes");
      flags |= DUAL_APPROACH;
      strcpy(rationalCone, "yes");
    }
    else if (strncmp(argv[i], "--ehrhart-taylor=", 17) == 0) {
      strcpy(taylor, "yes");
      degree = atoi(argv[i] + 17);
      strcpy(dualApproach,"yes");
      flags |= DUAL_APPROACH;
    }
    else if (strncmp(argv[i], "--avoid-singularities", 7) == 0) {
      params->shortvector = BarvinokParameters::SubspaceAvoidingLLL;
      params->triangulation = BarvinokParameters::SubspaceAvoidingRecursiveTriangulation;
    }
    else if (strncmp(argv[i], "--approximate", 7) == 0)
      approx = true;
    else if (strncmp(argv[i], "--sampling-factor=", 18) == 0)
      sampling_factor = atof(argv[i] + 18);
    else if (strncmp(argv[i], "--num-samples=", 14) == 0)
      num_samples = atol(argv[i] + 14);
    else if (strncmp(argv[i], "--random-seed=", 14) == 0) {
      unsigned int seed = atoi(argv[i] + 14);
      seed_random_generator(seed);
    }
    else {
      cerr << "Unknown command/option " << argv[i] << endl;
      exit(1);
    }
  }
  
  if(minimize[0] == 'y') strcpy(maximum, "yes");
  if(grobner[0] == 'y') strcpy(equationsPresent,"yes");
  if(binary[0] == 'y') {strcpy(maximum,"yes"); strcpy(Memory_Save, "no");}
  if(maximum[0] == 'y') strcpy(Memory_Save, "no");
  if(printfile[0] == 'y') strcpy(Memory_Save, "no");
  if(rationalCone[0] == 'y') strcpy(Memory_Save, "no");
  if(printfile[0] == 'y') print_flag = 1;

  if(rationalCone[0] == 'y'){
    
    //HugInt digit(argv[1]);
    //conv(output_cone, digit.BigInt);
    // User can use only Mode one
    output_cone = 3;
  }
  int dilation_const = 1;

  if(dilation[0] == 'y') dilation_const = atoi(argv[argc-2]);

  strcat(invocation,argv[argc-1]);
  strcat(invocation,"\n\n");
#if 1
  cout << invocation;
#endif

  if(output_cone > 3) output_cone = 0;
  flags |= (output_cone << 1);

  char costFile[127];
  if(maximum[0] == 'y'){
    strcpy(fileName,argv[argc-1]);
    strcpy(costFile, argv[argc - 1]);
    strcat(costFile, ".cost");
  }
  else strcpy(fileName,argv[argc-1]);
  //  cout << fileName << " " << costFile << endl;
  if(maximum[0] == 'y') {
    ifstream ReadTest(fileName);
    if(!ReadTest){
      cerr << "Need a polytope input file." << endl;
      exit(2);
    }
    //    cout << fileName << " " << costFile << endl;
    ifstream INCost(costFile);
    if(!INCost){
      cerr << "Need a cost input file." << endl;
      exit(3);
    }
    int costDim, dummy;
    INCost >> dummy >> costDim; 
    cost.SetLength(costDim);
    for(i = 0; i < costDim; i++)
      INCost >> cost[i]; 
  }
  //strcpy (fileName,"stdin");

  /* INPUT HANDLING. */

  Polyhedron *Poly = NULL;

  if(grobner[0] == 'y'){
    CheckGrobner(fileName, cddstyle);
    SolveGrobner(fileName,  nonneg, dualApproach,
		 grobner, equationsPresent, cddstyle);
    exit(0);
  }
  
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
    Poly = ReadLatteStyleVrep(fileName, dualApproach[0] == 'y');
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
    if((cddstyle[0] == 'n') && (grobner[0] == 'n') && (maximum[0] == 'n')&& (minimize[0] == 'n')){
      CheckInputFile(fileName);
      CheckLength(fileName, equationsPresent);
    }
    if(minimize[0] == 'y')  strcpy (maximum, "yes");   
 
    if((cddstyle[0] == 'n') && (grobner[0] == 'n') && (maximum[0] == 'y')){
      CheckInputFile(fileName);
      CheckLength(fileName,equationsPresent);
    }

    if(cddstyle[0] == 'y')
      { CheckInputFileCDDRep(argv[argc - 1]);
      CheckInputFileCDDRep1(argv[argc - 1]);
      CheckInputFileCDDRep3(argv[argc - 1]);
      CheckInputFileCDDRep4(argv[argc - 1]);
      }
  }else CheckInputFileVrep(fileName);
  if (cddstyle[0] == 'y') {
    cerr << "Warning: Not performing check for empty polytope, "
	 << "because it is unimplemented for the CDD-style input format. "
	 << endl;
  }
  else {
    CheckEmpty(fileName);
  }
  //vec_ZZ cost;
  /* Read problem data. */
  params->read_time.start();
  
  if((cddstyle[0] == 'n') && (Vrepresentation[0] == 'n')) CheckRed(fileName, equationsPresent, maximum, nonneg, interior, dilation, dilation_const); 

  dilation_const = 1;
  if((cddstyle[0] == 'n') && (grobner[0] == 'n'))
    readLatteProblem(fileName,&equations,&inequalities,equationsPresent,
		     &numOfVars, nonneg, dualApproach, grobner, maximum, 
		     cost,Vrepresentation);
//   if((equationsPresent[0] == 'n') && (interior[0] == 'y'))
//     Interior(inequalities);
  // if(minimize[0] == 'y') cost = -cost;
  if(cddstyle[0] == 'y'){
    int tmpoutput, tmpflags;
    CDDstylereadLatteProblem(fileName,&equations,&inequalities,equationsPresent,
			     &numOfVars, nonneg, dualApproach, taylor, degree,
			     rationalCone, tmpoutput, tmpflags, Memory_Save,
			     assumeUnimodularCones, inthull, grobner);
    output_cone = 3;
    flags = tmpflags;
  }
  
  // cout << grobner << endl;
  if(minimize[0] == 'y') cost = - cost;
  holdCost = cost;
  //cout <<"Cost is: " << cost << endl;
  vec_RR holdcost_RR;
  holdcost_RR.SetLength(holdCost.length());
  for(i = 0; i < holdCost.length(); i++) conv(holdcost_RR[i], holdCost[i]);

  if(minimize[0] == 'y') holdcost_RR = - holdcost_RR;

  if((dualApproach[0] == 'y') && (nonneg[0] == 'y')&&(equationsPresent[0] == 'n')){
    cerr<<"You cannot use + and dua at the same time." << endl;
    exit(2);
  }
  
  if((Memory_Save[0] == 'y') && (inthull[0] == 'y')){
    cerr<<"You cannot use int and memsave at the same time." << endl;
    exit(3);
  }
  
  numOfVars--;
  /* Binary seach IP*/

  if(binary[0] == 'y'){
    cout << "The number of optimal solutions: " << binarySearch(equations, inequalities,cost, numOfVars, minimize) << endl;
    cout << "Time: " << GetTime() << endl;
    exit(0);
  }

  ProjU.SetDims(numOfVars, numOfVars);
  ProjU2.SetDims(numOfVars, numOfVars);
  oldnumofvars = numOfVars;
  generators=createArrayVector(numOfVars);
  if (equationsPresent[0]=='y') {
    /*    if(grobner[0] == 'y')
	  {  
     matrixTmp=Grobner(equations,inequalities,&generators,&numOfVars, &templistVec, oldnumofvars);
     
     }*/
    matrixTmp=preprocessProblem(equations,inequalities,&generators,&numOfVars, cost, ProjU, interior, dilation_const);
    ProjU2 = transpose(ProjU);
    bb = ProjU2[0];
    AAA.SetDims(ProjU2.NumRows() - 1, ProjU2.NumCols());
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
   }
  else {
    matrix = matrixTmp;}
/* Now matrix contains the new inequalities. */
  params->read_time.stop();
  cout << params->read_time;
  
    //   cout << "Project down cost function: " << cost << endl;
    vec_RR Rat_solution, tmp_den, tmp_num;
    mat_RR ProjU_RR;
    ProjU_RR.SetDims(ProjU.NumRows(), ProjU.NumCols());
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
    if(LRS[0] == 'n')
    tmpcones=computeVertexCones(fileName,matrix,numOfVars);
    else
      tmpcones=computeVertexConesViaLrs(fileName,matrix,numOfVars);
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
       LP_OPT = Rat_cost*Rat_solution; //cout << cost << endl;
      cout << "The LP optimal value is: " << holdcost_RR*ProjectingUpRR(ProjU_RR, Rat_solution, numOfVars) << endl;
      }else {Poly->cones = tmpcones;
      cout << "The polytope has " << lengthListCone(Poly->cones) << " vertices." << endl;
    //system_with_error_check("rm -f numOfLatticePoints");
    }
  } 

  params->vertices_time.stop();
  cout << params->vertices_time;

  /* Compute triangulation or decomposition of each vertex cone. */

  if (dualApproach[0]=='y') {

    Poly->cones=createListCone();
    Poly->cones->vertex = new Vertex(createRationalVector(numOfVars));
    rays=createListVector(createVector(numOfVars));
    endRays=rays;
    tmpRays=matrix;
    while (tmpRays) {
      v=createVector(numOfVars);
      for (i=0; i<numOfVars; i++) v[i]=-(tmpRays->first)[i+1];
      endRays->rest=createListVector(v);
      endRays=endRays->rest;
      tmpRays=tmpRays->rest;
    }
    Poly->cones->rays = rays->rest;
    Poly->dualized = true;

    cout << "Homogenization: " << endl;
    printListCone(Poly->cones, numOfVars);
  }

  Poly->numOfVars = numOfVars;
  } /* Not VREP */

  
  
  if (ehrhart_polynomial) {
    /* Translate all cones to the origin, saving the original vertex. */
    listCone *cone;
    for (cone = Poly->cones; cone; cone = cone->rest) {
      ZZ scale_factor;
      cone->vertex->ehrhart_vertex
	= scaleRationalVectorToInteger(cone->vertex->vertex,
				       Poly->numOfVars, scale_factor);
      if (scale_factor != 1) {
	cerr << "Computation of Ehrhart polynomials is only implemented "
	     << "for integral polytopes." << endl
	     << "Use `ehrhart' for computing the Ehrhart series "
	     << "of rational polytopes." << endl;
	exit(1);
      }
      delete cone->vertex->vertex;
      cone->vertex->vertex = new rationalVector(Poly->numOfVars);
    }
  }
    
  params->Flags = flags;
  params->File_Name = fileName;
  params->Number_of_Variables = Poly->numOfVars;


  switch (params->decomposition) {
  case BarvinokParameters::DualDecomposition:
  case BarvinokParameters::IrrationalPrimalDecomposition:
    if (not Poly->dualized) {
      params->dualize_time.start();
      Poly->cones = dualizeCones(Poly->cones, Poly->numOfVars);
      params->dualize_time.stop();
      cout << params->dualize_time;
      Poly->dualized = true;
    }
    break;
  case BarvinokParameters::IrrationalAllPrimalDecomposition:
    if (Poly->dualized) {
      cerr << "You cannot use `homog' and `--all-primal' at the same time." << endl;
      exit(3);
    }
#if 1
    cout << "Irrationalizing polyhedral cones... ";
    cout.flush();
    params->dualize_time.start();
    if (Vrepresentation[0] == 'y') {
      cout << "(First computing facets for them... "; cout.flush();
      Poly->cones = dualizeCones(Poly->cones, Poly->numOfVars);
      Poly->cones = dualizeBackCones(Poly->cones, Poly->numOfVars); // just swaps
      cout << "done; sorry for the interruption.) "; cout.flush();
    }      
    else {
      /* Fill in the facets of all cones; we determine them by
	 taking all inequalities tight at the respective vertex. */
      computeTightInequalitiesOfCones(Poly->cones, matrix, Poly->numOfVars);
    }
    params->dualize_time.stop(); cout << params->dualize_time;
    params->irrationalize_time.start();
    {
      listCone *cone;
      for (cone = Poly->cones; cone; cone=cone->rest)
	assert(lengthListVector(cone->facets) >= Poly->numOfVars);
    }
    irrationalizeCones(Poly->cones, Poly->numOfVars);
    params->irrationalize_time.stop();
    cout << params->irrationalize_time;
#endif
    break;
  default:
    cerr << "Unknown BarvinokParameters::decomposition" << endl;
    abort();
  }

  try {
    
  switch (params->substitution) {
  case BarvinokParameters::PolynomialSubstitution:
    if (ehrhart_polynomial) {
      cerr << "Computation of Ehrhart polynomials is only implemented "
	   << "for the exponential substitution." << endl;
      exit(1);
    }
    if (assumeUnimodularCones[0]=='n') {
      if (Memory_Save[0] == 'n') {
	Poly->cones=decomposeCones(Poly->cones, not Poly->dualized,
				   *params);
	/* Compute points in parallelepipeds */
	computePointsInParallelepipeds(Poly->cones, Poly->numOfVars);
      }
      // Iterator through simplicial cones, DFS
      else {
	Standard_Single_Cone_Parameters *standard_params
	  = new Standard_Single_Cone_Parameters(*params);
	delete params; params = standard_params;
	decomposeAndComputeResidue(Poly->cones, degree, false,
				   *standard_params);
      }
    }
    break;
  case BarvinokParameters::ExponentialSubstitution:
    if (dualApproach[0] == 'y') {
      cerr << "Exponential substitution is not yet implemented for the homogenized version."
	   << endl;
      exit(1);
    }
    else {
      params->decompose_time.start();
      if (approx) {
#ifdef HAVE_EXPERIMENTS
	{
	  Write_Exponential_Sample_Formula_Single_Cone_Parameters *write_param
	    = new Write_Exponential_Sample_Formula_Single_Cone_Parameters
	    (*params, "Exponential_Sample_Formula", sampling_factor,
	     num_samples);
	  delete params;
	  params = write_param;
	  decomposeAndWriteExponentialSampleFormula(Poly->cones, *write_param);
	}
#else
	cerr << "Approximation code is not compiled in, sorry." << endl;
	exit(1);
#endif
      }
      else if (ehrhart_polynomial) {
	Exponential_Ehrhart_Parameters *exp_param
	  = new Exponential_Ehrhart_Parameters(*params);
	delete params;
	params = exp_param;
	mpq_vector ehrhart_coefficients
	  = decomposeAndComputeEhrhartPolynomial(Poly->cones, *exp_param);
	cout << endl << "Ehrhart polynomial: ";
	{
	  int i;
	  for (i = 0; i<ehrhart_coefficients.size(); i++) {
	    if (ehrhart_coefficients[i] > 0) 
	      cout << " + " << ehrhart_coefficients[i] << " * t^" << i;
	    else if (ehrhart_coefficients[i] < 0)
	      cout << " - " << abs(ehrhart_coefficients[i]) << " * t^" << i;
	  }
	}
	cout << endl << endl;
      }
      else {
	Exponential_Single_Cone_Parameters *exp_param
	  = new Exponential_Single_Cone_Parameters(*params);
	delete params;
	params = exp_param;
	Integer number_of_lattice_points
	  = decomposeAndComputeExponentialResidue(Poly->cones, *exp_param);
	cout << endl << "****  The number of lattice points is: "
	     << number_of_lattice_points << "  ****" << endl << endl;
	// FIXME: Centralize this output stuff.
	ofstream out("numOfLatticePoints");
	out << number_of_lattice_points << endl;
      }
      params->decompose_time.stop();
    }
    break;
  default:
    cerr << "Unknown BarvinokParameters::substitution" << endl;
    abort();
  }

  if(grobner[0] == 'y'){

 Poly->cones = ProjectUp(Poly->cones, oldnumofvars, Poly->numOfVars, templistVec);
 Poly->numOfVars = oldnumofvars;

  }
 if(Print[0] == 'y')
  printListCone(Poly->cones,Poly->numOfVars);

 if(inthull[0] == 'y')
   ;
 //   printListVector(IntegralHull(Poly->cones,  inequalities, equations, Poly->numOfVars), Poly->numOfVars);
 if(maximum[0] == 'y') {
   listCone * Opt_cones;
   if(Singlecone[0] == 'n'){
   Opt_cones = CopyListCones(Poly->cones, Poly->numOfVars);
   ZZ NumOfLatticePoints; //printListCone(Opt_cones, Poly->numOfVars);
   NumOfLatticePoints = Residue(Opt_cones, Poly->numOfVars);
   cout <<"Finished computing a rational function. " << endl;
   cout <<"Time: " << GetTime() << " sec." << endl;
   if(IsZero(NumOfLatticePoints) == 1){
     cerr<<"Integrally empty polytope.  Check the right hand side."<< endl;
     exit(0);}
      else{
	int singleCone = 0;
	if(Singlecone[0] == 'y') singleCone = 1;
	vec_ZZ Opt_solution; 
	if(minimize[0] == 'y') holdCost = -holdCost;
	Opt_solution = SolveIP(Poly->cones, matrix, cost, Poly->numOfVars, singleCone); 
        if(minimize[0] == 'y') cost = -cost;
	cout << "An optimal solution for " <<  holdCost << " is: " << ProjectingUp(ProjU, Opt_solution, Poly->numOfVars) << "." << endl;
	cout << "The projected down opt value is: " << cost * Opt_solution << endl;
	cout <<"The optimal value is: " << holdCost * ProjectingUp(ProjU, Opt_solution, Poly->numOfVars) << "." << endl;
	ZZ IP_OPT; IP_OPT = cost*Opt_solution;
	RR tmp_RR;

	conv(tmp_RR, cost * Opt_solution);
	// cout << tmp_RR << " " << LP_OPT << endl;
	if(minimize[0] == 'y') LP_OPT = - LP_OPT;
	cout <<"The gap is: "<< abs(tmp_RR - LP_OPT) << endl;
	cout << "Computation done." << endl;
	cout <<"Time: " << GetTime() << " sec." << endl;
        strcpy(command,"rm -f ");
        strcat(command,fileName);
        strcat(command,".ext");
        system_with_error_check(command);

        strcpy(command,"rm -f ");
        strcat(command,fileName);
        strcat(command,".cdd");
        system_with_error_check(command);

        strcpy(command,"rm -f ");
        strcat(command,fileName);
        strcat(command,".ead");
        system_with_error_check(command);

	if(cddstyle[0] == 'n' && Vrepresentation[0] == 'n'){
        strcpy(command,"rm -f ");
        strcat(command,fileName);
        system_with_error_check(command);
	}

	exit(0);
      }
   }
   else{
	int singleCone = 0;
	if(Singlecone[0] == 'y') singleCone = 1;
	vec_ZZ Opt_solution; 
	if(minimize[0] == 'y') holdCost = -holdCost;
	Opt_solution = SolveIP(Poly->cones, matrix, cost, Poly->numOfVars, singleCone); 
	cout << "An optimal solution for " <<  holdCost << " is: " << ProjectingUp(ProjU, Opt_solution, Poly->numOfVars) << "." << endl;
        if(minimize[0] == 'y') cost = -cost;
	cout << "The projected down opt value is: " << cost * Opt_solution << endl;
	cout <<"The optimal value is: " << holdCost * ProjectingUp(ProjU, Opt_solution, Poly->numOfVars) << "." << endl;
	ZZ IP_OPT; IP_OPT = cost*Opt_solution;
	RR tmp_RR;
	conv(tmp_RR, IP_OPT);
	// cout << cost * Opt_solution << endl;
	if(minimize[0] == 'y') LP_OPT = - LP_OPT;
	cout <<"The gap is: "<< abs(tmp_RR - LP_OPT) << endl;
	cout << "Computation done." << endl;
	cout <<"Time: " << GetTime() << " sec." << endl;
        strcpy(command,"rm -f ");
        strcat(command,fileName);
        strcat(command,".ext");
        system_with_error_check(command);

        strcpy(command,"rm -f ");
        strcat(command,fileName);
        strcat(command,".cdd");
        system_with_error_check(command);

        strcpy(command,"rm -f ");
        strcat(command,fileName);
        strcat(command,".ead");
        system_with_error_check(command);

	if(cddstyle[0] == 'n' && Vrepresentation[0] == 'n'){
        strcpy(command,"rm -f ");
        strcat(command,fileName);
        system_with_error_check(command);
	}

	exit(0);
   }
 }else{
if(Memory_Save[0] == 'n')
{

	if(dualApproach[0] == 'n'){
	  cout << "Creating generating function.\n"; 
	  //printListVector(templistVec, oldnumofvars); cout << ProjU << endl;
	if(equationsPresent[0] == 'y'){ Poly->cones = ProjectUp2(Poly->cones, oldnumofvars, Poly->numOfVars, AA, bb);
	Poly->numOfVars = oldnumofvars;}

	  createGeneratingFunctionAsMapleInput(fileName,Poly->cones,Poly->numOfVars);  }
        //printListCone(cones, Poly->numOfVars);

	cout << "Printing decomposed cones to decomposed_cones." << endl;
	printListConeToFile("decomposed_cones", Poly->cones, Poly->numOfVars);

	if(dualApproach[0] == 'n'){
	cout << "Starting final computation.\n";
	cout << endl << "****  The number of lattice points is: " << Residue(Poly->cones,Poly->numOfVars) << "  ****" << endl << endl;}


	if(dualApproach[0] == 'y')
	{
	  cout << "Starting final computation.\n";
	  //cout << "output_cone: " << output_cone;
	  switch (params->decomposition) {
	  case BarvinokParameters::IrrationalPrimalDecomposition:
	  case BarvinokParameters::IrrationalAllPrimalDecomposition: {
#ifdef HAVE_EXPERIMENTS
	    ofstream out("func.rat");
	    out << "HS := ";
	    TrivialMonomialSubstitutionMapleOutput(out, Poly->cones, Poly->numOfVars);
	    out << ";";
#else
	    cerr << "Trivial monomial subtitution not compiled in, sorry." << endl;
#endif
	    break;
	  }
	  case BarvinokParameters::DualDecomposition:
	    ResidueFunction(Poly->cones,Poly->numOfVars, print_flag, degree, output_cone);
	    break;
	  default:
	    assert(0);
	  }
	//  Else we have already computed the residue.
	}

}
 }
  } catch (NotIrrationalException) {
    cerr << "Bug: Irrationalization failed" << endl;
    exit(1);
  };

 if(rationalCone[0] == 'y') {
   cout << endl <<"Rational function written to " << argv[argc - 1] << ".rat" << endl << endl;
   strcpy(command, "mv ");
   strcat(command, "simplify.sum ");
   strcat(command, argv[argc - 1]);
   strcat(command, ".rat");
   system_with_error_check(command);
 }

 if(printfile[0] == 'y'){
   cout << endl <<"Rational function written to " << argv[argc - 1] << ".rat" << endl << endl;
   strcpy(command, "mv ");
   strcat(command, "func.rat ");
   strcat(command, argv[argc - 1]);
   strcat(command, ".rat");
   system_with_error_check(command);
 }
 if((removeFiles[0] == 'y') && (dualApproach[0] == 'n')){
   
  strcpy(command,"rm -f ");
  strcat(command,fileName);
  strcat(command,".ext");
  system_with_error_check(command);
  
  strcpy(command,"rm -f ");
  strcat(command,fileName);
  strcat(command,".cdd");
  system_with_error_check(command); 
  
  if(Memory_Save[0] == 'n'){
    strcpy(command,"rm -f ");
    strcat(command,fileName);
    strcat(command,".maple");
    system_with_error_check(command); 
  }

  strcpy(command,"rm -f ");
  strcat(command,fileName);
  strcat(command,".ead");
  system_with_error_check(command); 
  
  if(cddstyle[0] == 'n' && Vrepresentation[0] == 'n'){
    strcpy(command,"rm -f ");
    strcat(command,fileName);
    system_with_error_check(command);
  }
 }

  if((dualApproach[0] == 'y') && (cddstyle[0] == 'n') && Vrepresentation[0] == 'n'){

    strcpy(command,"rm -f ");
    strcat(command,fileName);
    system_with_error_check(command);

  }
  //cout << "Computation done. " << endl;

  params->total_time.stop();
  cout << params->total_time;
  
  {
    // until we have a more sophisticated test script --mkoeppe
    ofstream totalTime("totalTime");
    totalTime << params->total_time.get_seconds()
	      << " (" << params->read_time.get_seconds() << "r"
	      << ", " << params->vertices_time.get_seconds() << "v"
	      << ", " << params->irrationalize_time.get_seconds() << "i"
	      << ", " << params->dualize_time.get_seconds() << "d"
	      << ", " << params->triangulate_time.get_seconds() << "t"
	      << ", " << params->decompose_time.get_seconds() << "b"
	      << ")" << endl;
    ofstream stats("latte_stats");
    params->print_statistics(stats);
  }
 
 return(0);
}
/* ----------------------------------------------------------------- */



