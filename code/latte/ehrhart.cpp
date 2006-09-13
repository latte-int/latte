/* ----------------------------------------------------------------- */
/*                                                                   */
/* LattE (Lattice Point Enumeration)                                 */
/*                                                                   */
/* Master program                                                    */
/*                                                                   */
/* Author     : Raymond Hemmecke, Ruriko Yoshida                     */
/*                                                                   */
/* Created    : 07-JUN-02                                            */
/* Last Update: 03-Mar-03                                            */
/*                                                                   */
/* ----------------------------------------------------------------- */

#include <string.h>
#include <stdio.h>

#include "myheader.h"
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

#include "banner.h"
#include "convert.h"
#include "latte_system.h"

/* ----------------------------------------------------------------- */
int main(int argc, char *argv[]) {
#ifdef SUN
  struct tms tms_buf;
#endif
  float z;
  int i,numOfVars,numOfAllVars, degree = 1;
  unsigned int flags = 0, print_flag = 0, output_cone = 0;
  vec_ZZ dim, v, w;
  int oldnumofvars;
  vec_ZZ *generators;
  char fileName[127], invocation[127], decompose[10], equationsPresent[10],
    assumeUnimodularCones[127], dualApproach[127], taylor[127], printfile[127],
    rationalCone[127], nonneg[127], Memory_Save[127], Save_Tri[127],
    Load_Tri[127], Print[127], inthull[127], cddstyle[127], grobner[127],
    removeFiles[127], command[127], maximum[127],  Singlecone[127], LRS[127],
    Vrepresentation[127], dilation[127], minimize[127], binary[127], 
    interior[127];
  listVector *matrix, *equations, *inequalities, *rays, *endRays, *tmpRays, *matrixTmp;
  vec_ZZ cost;
  listVector *templistVec;
  listCone *cones, *tmp, *tmpcones;

  latte_banner(cout);

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
  strcpy(removeFiles,"no");
  strcpy(grobner,"no");
  strcpy(maximum,"no");
  strcpy(minimize,"no");
  strcpy(decompose,"yes");
  strcpy(dualApproach,"yes");
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

  flags |= DUAL_APPROACH;
 
  for (i=1; i<argc-1; i++) {
    strcat(invocation,argv[i]);
    strcat(invocation," ");
    //if (strncmp(argv[i],"vrep",3)==0) strcpy(Vrepresentation,"yes"); 
    // if (strncmp(argv[i],"bbs",3)==0) strcpy(binary,"yes");
    //if (strncmp(argv[i],"int",3)==0) strcpy(interior,"yes");
    //if (strncmp(argv[i],"min",3)==0) strcpy(minimize,"yes");
    //if (strncmp(argv[i],"gro",3)==0) strcpy(grobner,"yes");
    //if (strncmp(argv[i],"nodecom",3)==0) strcpy(decompose,"no");
    //if (strncmp(argv[i],"homog",3)==0) {strcpy(dualApproach,"yes"); flags |= DUAL_APPROACH;}
    //if (strncmp(argv[i],"equ",3)==0) strcpy(equationsPresent,"yes");
    //if (strncmp(argv[i],"uni",3)==0) strcpy(assumeUnimodularCones,"yes");
    //if (strncmp(argv[i],"ehrh",4)==0) {strcpy(printfile,"yes"); flags |= PRINT;}
    //if(strncmp(argv[i],"file",4)==0) strcpy(Memory_Save, "no");
    //if (strncmp(argv[i],"single",6)==0) strcpy(Singlecone,"yes");
    // if (strncmp(argv[i],"tay",3)==0) strcpy(taylor,"yes");
    if (strncmp(argv[i],"simp",3)==0) strcpy(rationalCone,"yes");
    //if (strncmp(argv[i],"+", 1) ==0) strcpy(nonneg,"yes");
    //if (strncmp(argv[i],"memsave",7)==0) strcpy (Memory_Save, "yes");
    //if (strncmp(argv[i],"printcones",3)==0) strcpy (Print, "yes");
    if (strncmp(argv[i],"cdd",3)==0) strcpy (cddstyle, "yes");
    //if (strncmp(argv[i],"hull",3)==0) strcpy (inthull, "yes");
    //if (strncmp(argv[i],"max",3)==0) strcpy (maximum, "yes");
    if (strncmp(argv[i],"lrs",3)==0) strcpy (LRS, "yes");
    //if (strncmp(argv[i],"dil",3)==0) strcpy (dilation, "yes");
    //if (strncmp(argv[i],"rem",3)==0) strcpy (removeFiles, "yes");
    if (strncmp(argv[i],"trisave",7)==0) {strcpy (Save_Tri, "yes"); flags |= SAVE;}
    if (strncmp(argv[i],"triload",7)==0) {strcpy (Load_Tri, "yes"); flags |= LOAD;}
    } 

  /* if cdd and lrs are NOT called. */
  if((cddstyle[0] == 'n') && (LRS[0] == 'n')){
    if(argc == 2) {strcpy(printfile,"yes"); flags |= PRINT;}
    else if(argc == 3) {
      if(strncmp(argv[1],"sim",3)==0) strcpy(rationalCone,"yes");
      else strcpy(taylor,"yes");
    }else{
      cerr << "Too many arguments.  Check the manual for command line." << endl;
      exit(1);
    }
  }

  /* if cdd is called but and lrs is NOT called. */
  else if((cddstyle[0] == 'y') && (LRS[0] == 'n')){
    if(argc == 3) {strcpy(printfile,"yes"); flags |= PRINT;}
    else if(argc == 4) {
      for(i = 1; i < 3; i++)
	if(strncmp(argv[i],"sim",3)==0) strcpy(rationalCone,"yes");
      if (rationalCone[0] == 'n') strcpy(taylor,"yes");
    }else{
      cerr << "Too many arguments.  Check the manual for command line." << endl;
      exit(1);
    }
  }
  

  /* if cdd is NOT called and lrs is called. */
  else if((cddstyle[0] == 'n') && (LRS[0] == 'y')){
    if(argc == 3) {strcpy(printfile,"yes"); flags |= PRINT;}
    else if(argc == 4) {
      for(i = 1; i < 3; i++)
       if(strncmp(argv[i],"sim",3)==0) strcpy(rationalCone,"yes");
      if (rationalCone[0] == 'n') strcpy(taylor,"yes");
    }else{
      cerr << "Too many arguments.  Check the manual for command line." << endl;
      exit(1);
    }
  }

  /* if cdd and lrs are called. */
  else if((cddstyle[0] == 'y') && (LRS[0] == 'y')){
    if(argc == 4) {strcpy(printfile,"yes"); flags |= PRINT;}
    else if(argc == 5) {
      for(i = 1; i < 4; i++)
	if(strncmp(argv[i],"sim",3)==0) strcpy(rationalCone,"yes");
      // else if(strncmp(argv[3],"sim",3)==0) strcpy(rationalCone,"yes");
      if (rationalCone[0] == 'n')  strcpy(taylor,"yes");
    }else{
      cerr << "Too many arguments.  Check the manual for command line." << endl;
      exit(1);
    }
  }



  flags |= DUAL_APPROACH;
  if(minimize[0] == 'y') strcpy(maximum, "yes");
  if(grobner[0] == 'y') strcpy(equationsPresent,"yes");
  if(binary[0] == 'y') {strcpy(maximum,"yes"); strcpy(Memory_Save, "no");}
  if(maximum[0] == 'y') strcpy(Memory_Save, "no");
  if(printfile[0] == 'y') strcpy(Memory_Save, "no");
  if(rationalCone[0] == 'y') strcpy(Memory_Save, "no");
  if(printfile[0] == 'y') print_flag = 1;
  if(taylor[0] == 'y'){
    degree = atoi(argv[argc - 2]);
  }
  if(rationalCone[0] == 'y'){
    
    //HugInt digit(argv[1]);
    //conv(output_cone, digit.BigInt);
    // User can use only Mode one
    output_cone = 3;
  }
  int dilation_const = 1;

  if(dilation[0] == 'y') dilation_const = atoi(argv[argc-2]);

  if(output_cone > 3) output_cone = 0;
  flags |= (output_cone << 1);
  
  if((dualApproach[0] == 'y') && (nonneg[0] == 'y')&&(equationsPresent[0] == 'n')){
    cerr<<"You cannot use + and dua at the same time." << endl;
    exit(2);
  }

  if((cddstyle[0] == 'y') && (Vrepresentation[0] == 'y')){
    cerr << "Use not cdd style and v-representation." << endl;
    exit(2);
  }
  
  if((Memory_Save[0] == 'y') && (inthull[0] == 'y')){
    cerr<<"You cannot use int and memsave at the same time." << endl;
    exit(3);
  }

  strcat(invocation,argv[argc-1]);
  strcat(invocation,"\n\n");
  cout << invocation;
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
    CheckEmpty(fileName);
  //vector cost;
  /* Read problem data. */
  if((cddstyle[0] == 'n') && (Vrepresentation[0] == 'n'))
   CheckRed(fileName, equationsPresent, maximum, nonneg, interior, dilation, dilation_const); 

  dilation_const = 1;
  if((cddstyle[0] == 'n') && (grobner[0] == 'n'))
    readLatteProblem(fileName,&equations,&inequalities,equationsPresent,
		     &numOfVars, nonneg, dualApproach, grobner, maximum, 
		     cost,Vrepresentation);
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
  vec_ZZ holdCost;
  if(minimize[0] == 'y') cost = - cost;
  holdCost = cost;
  //cout <<"Cost is: " << cost << endl;
  vec_RR holdcost_RR;
  holdcost_RR.SetLength(holdCost.length());
  for(i = 0; i < holdCost.length(); i++) conv(holdcost_RR[i], holdCost[i]);

  if(minimize[0] == 'y') holdcost_RR = - holdcost_RR;
  if(grobner[0] == 'y'){

    CheckGrobner(fileName, cddstyle);
    SolveGrobner(fileName,  nonneg, dualApproach,
		 grobner, equationsPresent, cddstyle);}
  else{
  if((dualApproach[0] == 'y') && (nonneg[0] == 'y')&&(equationsPresent[0] == 'n')){
    cerr<<"You cannot use + and dua at the same time." << endl;
    exit(2);
  }
  
  if((Memory_Save[0] == 'y') && (inthull[0] == 'y')){
    cerr<<"You cannot use int and memsave at the same time." << endl;
    exit(3);
  }
  
  if((Vrepresentation[0] == 'y') && (equationsPresent[0] == 'y')){
    cerr<<"You cannot use vrep and equ at the same time." << endl;
    exit(4);
  }

  numOfVars--;
  /* Binary seach IP*/

  if(binary[0] == 'y'){
    cout << "The number of optimal solutions: " << binarySearch(equations, inequalities,cost, numOfVars, minimize) << endl;
    cout << "Time: " << GetTime() << endl;
    exit(0);
  }

  numOfAllVars=numOfVars;
  mat_ZZ ProjU;
  ProjU.SetDims(numOfVars, numOfVars);
  oldnumofvars = numOfVars;
  generators=createArrayVector(numOfVars);
  if (equationsPresent[0]=='y') {
    /*    if(grobner[0] == 'y')
	  {  
     matrixTmp=Grobner(equations,inequalities,&generators,&numOfVars, &templistVec, oldnumofvars);
     
     }*/
    matrixTmp=preprocessProblem(equations,inequalities,&generators,&numOfVars, cost, ProjU, interior, dilation_const);
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
  RR LP_OPT;
    cout << "\nTime: " << GetTime() << " sec\n\n";
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
       cones = CopyListCones(tmpcones, numOfVars, LP_vertex);
     else cones = tmpcones;
       if(lengthListCone(cones) == 1) 
	cout <<"\nWe found a single vertex cone for IP.\n" << endl;
       cout <<"A vertex which we found via LP is: " << ProjectingUpRR(ProjU_RR, Rat_solution, numOfVars) << endl;
      //printRationalVector(LP_vertex, numOfVars);
       LP_OPT = Rat_cost*Rat_solution; //cout << cost << endl;
      cout << "The LP optimal value is: " << holdcost_RR*ProjectingUpRR(ProjU_RR, Rat_solution, numOfVars) << endl;
      }else {cones = tmpcones;
    cout << "\nThe polytope has " << lengthListCone(cones) << " vertices.\n";
    //system_with_error_check("rm -f numOfLatticePoints");
    cout << endl;}
  } 

  /* Reading from the vertex representation. */

  if(Vrepresentation[0] == 'y')
     cones=computeVertexConesFromVrep(fileName,numOfVars); 

/* Compute triangulation or decomposition of each vertex cone. */

  if (dualApproach[0]=='y') {
    cones=createListCone();
    cones->vertex=createRationalVector(numOfVars);
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
    cones->rays=rays->rest;


	if (Memory_Save[0] == 'n' )
	{
    		cones=decomposeCones(cones,numOfVars, flags, fileName,
				     1, false,
				     BarvinokParameters::DualDecomposition); 
  	}

	else
	{
		decomposeCones_Single(cones, numOfVars, degree, flags, fileName, 1, false,
				      BarvinokParameters::DualDecomposition);
	}

 } 
  
   else 
   {

    if (assumeUnimodularCones[0]=='n') {
      if (decompose[0]=='y') 
      {
	if(Memory_Save[0] == 'n')      
	  cones=decomposeCones(cones,numOfVars, flags, fileName,
			       1, true,
			       BarvinokParameters::DualDecomposition);
	// Iterator through simplicial cones, DFS
	else
	  decomposeCones_Single(cones,numOfVars, degree, flags, fileName,
				1, true,
				BarvinokParameters::DualDecomposition);	
      }
    }
   }

/* Compute points in parallelepipeds, unless we already did using memsave version!  */


 if(Memory_Save[0] == 'n')
 {
	 cout << "Computing the points in the Parallelepiped of the unimodular Cones." << endl;
  	tmp=cones;
	int	Cones_Processed_Count = 0;
 	while (tmp) 
	{
	        tmp->latticePoints=pointsInParallelepiped(tmp, numOfVars);
    		tmp=tmp->rest;

		Cones_Processed_Count++;

		if ((Cones_Processed_Count % 1000) == 0 )
			cout << Cones_Processed_Count << " cones processed." << endl;
  	}
}
  



  if(grobner[0] == 'y'){

 cones = ProjectUp(cones, oldnumofvars, numOfVars, templistVec);
 numOfVars = oldnumofvars;

  }
 if(Print[0] == 'y')
  printListCone(cones,numOfVars);

 if(inthull[0] == 'y')
   ;
 //    printListVector(IntegralHull(cones, inequalities, equations,  numOfVars),numOfVars);

  if(maximum[0] == 'y') {
   listCone * Opt_cones;
   if(Singlecone[0] == 'n'){
   Opt_cones = CopyListCones(cones, numOfVars);
   ZZ NumOfLatticePoints; //printListCone(Opt_cones, numOfVars);
   NumOfLatticePoints = Residue(Opt_cones, numOfVars);
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
	Opt_solution = SolveIP(cones, matrix, cost, numOfVars, singleCone); 
        if(minimize[0] == 'y') cost = -cost;
	cout << "An optimal solution for " <<  holdCost << " is: " << ProjectingUp(ProjU, Opt_solution, numOfVars) << "." << endl;
	cout << "The projected down opt value is: " << cost * Opt_solution << endl;
	cout <<"The optimal value is: " << holdCost * ProjectingUp(ProjU, Opt_solution, numOfVars) << "." << endl;
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

	if(cddstyle[0] == 'n'){
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
	Opt_solution = SolveIP(cones, matrix, cost, numOfVars, singleCone); 
	cout << "An optimal solution for " <<  holdCost << " is: " << ProjectingUp(ProjU, Opt_solution, numOfVars) << "." << endl;
        if(minimize[0] == 'y') cost = -cost;
	cout << "The projected down opt value is: " << cost * Opt_solution << endl;
	cout <<"The optimal value is: " << holdCost * ProjectingUp(ProjU, Opt_solution, numOfVars) << "." << endl;
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

	if(cddstyle[0] == 'n'){
        strcpy(command,"rm -f ");
        strcat(command,fileName);
        system_with_error_check(command);
	}

	exit(0);
   }
 }else{
   /* if(maximum[0] == 'y') {
   listCone * Opt_cones;
   if(Singlecone[0] == 'n'){
   Opt_cones = CopyListCones(cones, numOfVars);
   ZZ NumOfLatticePoints; //printListCone(Opt_cones, numOfVars);
   NumOfLatticePoints = Residue(Opt_cones, numOfVars);
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
	Opt_solution = SolveIP(cones, inequalities, equations,  cost, numOfVars, singleCone); 
	cout << "An optimal solution for " <<  holdCost << " is: " << ProjectingUp(ProjU, Opt_solution, numOfVars) << "." << endl;
	cout << "The projected down opt value is: " << cost * Opt_solution << endl;
	cout <<"The optimal value is: " << holdCost * ProjectingUp(ProjU, Opt_solution, numOfVars) << "." << endl;
	ZZ IP_OPT; IP_OPT = cost*Opt_solution;
	RR tmp_RR;
	conv(tmp_RR, IP_OPT);
	// cout << cost * Opt_solution << endl;
	cout <<"The gap is: "<< abs(tmp_RR - LP_OPT) << endl;
	cout << "Computation done." << endl;
	cout <<"Time: " << GetTime() << " sec." << endl;
	exit(0);
      }
   }
   else{
	int singleCone = 0;
	if(Singlecone[0] == 'y') singleCone = 1;
	vec_ZZ Opt_solution; 
	if(minimize[0] == 'y') holdCost = -holdCost;
	Opt_solution = SolveIP(cones, inequalities, equations,  cost, numOfVars, singleCone); 
	cout << "An optimal solution for " <<  holdCost << " is: " << ProjectingUp(ProjU, Opt_solution, numOfVars) << "." << endl;
	cout << "The projected down opt value is: " << cost * Opt_solution << endl;
	cout <<"The optimal value is: " << holdCost * ProjectingUp(ProjU, Opt_solution, numOfVars) << "." << endl;
	ZZ IP_OPT; IP_OPT = cost*Opt_solution;
	RR tmp_RR;
	conv(tmp_RR, IP_OPT);
	// cout << cost * Opt_solution << endl;
	cout <<"The gap is: "<< abs(tmp_RR - LP_OPT) << endl;
	cout << "Computation done." << endl;
	cout <<"Time: " << GetTime() << " sec." << endl;
	exit(0);
   }
   }else{*/
if(Memory_Save[0] == 'n')
{

	if(dualApproach[0] == 'n'){
	  cout << "Creating generating function.\n"; 
	  //printListVector(templistVec, oldnumofvars); cout << ProjU << endl;
	if(equationsPresent[0] == 'y'){ cones = ProjectUp(cones, oldnumofvars, numOfVars, templistVec);
	numOfVars = oldnumofvars;}

	  createGeneratingFunctionAsMapleInput(fileName,cones,numOfVars);  }
        //printListCone(cones, numOfVars);
	if(dualApproach[0] == 'n'){
	cout << "Starting final computation.\n";
	cout << endl << "****  The number of lattice points is: " << Residue(cones,numOfVars) << "  ****" << endl << endl;}


	if(dualApproach[0] == 'y')
	{
		cout << "Starting final computation.\n";
		//cout << "output_cone: " << output_cone;
  		ResidueFunction(cones,numOfVars, print_flag, degree, output_cone);
	//  Else we have already computed the residue.
	}

}
 }

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
 if(removeFiles[0] == 'y'){
   
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

  if(cddstyle[0] == 'n'){
  strcpy(command,"rm -f ");
  strcat(command,fileName);
  strcat(command,".ead");
  system_with_error_check(command); 
  }
 }
  }

  if(cddstyle[0] == 'n'){
    strcpy(command,"rm -f ");
    strcat(command,fileName);
    system_with_error_check(command);
  }

 cout << "Computation done. " << endl;
 cout << "Time: " << GetTime() << " sec\n\n";
 
 return(0);
}
/* ----------------------------------------------------------------- */



