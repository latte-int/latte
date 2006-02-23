/* ----------------------------------------------------------------- */
/*                                                                   */
/* LattE (Lattice Point Enumeration)                                 */
/*                                                                   */
/* Master program                                                    */
/*  -- simplified version (mkoeppe, margulies) --                    */
/*                                                                   */
/* Author     : Raymond Hemmecke, Ruriko Yoshida                     */
/*                                                                   */
/* Created    : 07-JUN-02                                            */
/* Last Update: 03-Mar-03                                            */
/*                                                                   */
/* ----------------------------------------------------------------- */

#include <string.h>
#include <stdio.h>
#include <NTL/mat_ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_RR.h>
#include <NTL/vec_RR.h>
#include <NTL/config.h>
#include <NTL/LLL.h>
#include <NTL/HNF.h>
#include <NTL/ZZ.h>

#include "myheader.h"
#include "barvinok/dec.h"
#include "barvinok/barvinok.h"
#include "barvinok/Cone.h"
#include "barvinok/ConeDecom.h"
#include "barvinok/Triangulation.h"
#include "vertices/cdd.h"
#include "genFunction/maple.h"
#include "genFunction/piped.h"
#include "cone.h"
#include "ConeDeterminant.h"
#include "dual.h"
#include "RudyResNTL.h"
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

/* ----------------------------------------------------------------- */
int main(int argc, char *argv[]) {
#ifdef SUN
  struct tms tms_buf;
#endif
  float z;
  int i,numOfVars,numOfAllVars, degree = 1;
  unsigned int flags = 0, print_flag = 0, output_cone = 0;
  vector dim, v, w;
  int oldnumofvars;
  vector *generators;
  char fileName[127], invocation[127], decompose[10], equationsPresent[10],
    assumeUnimodularCones[127], dualApproach[127], taylor[127], printfile[127],
    rationalCone[127], nonneg[127], Memory_Save[127], Save_Tri[127],
    Load_Tri[127], Print[127], inthull[127], cddstyle[127], grobner[127],
    removeFiles[127], command[127], maximum[127],  Singlecone[127], LRS[127],
    Vrepresentation[127], dilation[127], minimize[127], binary[127], interior[127];
  listVector *matrix, *equations, *inequalities, *rays, *endRays, *tmpRays, *matrixTmp;
  vector cost;
  listVector *templistVec;
  listCone *cones, *tmpcones;

  latte_banner(cout);

  z=0;
  setbuf(stdout,0);

  strcpy(invocation,"Invocation: ");
  strcat(invocation,argv[0]);
  strcat(invocation," ");

  strcpy(Vrepresentation,"no");
  strcpy(interior,"no");
  strcpy(dilation,"no");
  strcpy(binary,"no");
  strcpy(Singlecone,"no");
  strcpy(removeFiles,"yes");
  strcpy(grobner,"no");
  strcpy(maximum,"no");
  strcpy(minimize,"no");
  strcpy(decompose,"yes");
  strcpy(dualApproach,"no");
  strcpy(equationsPresent,"no");
  strcpy(assumeUnimodularCones,"no");
  strcpy(printfile,"no");
  strcpy(taylor,"no");
  strcpy(rationalCone,"no");
  strcpy(nonneg, "no");
  strcpy(Memory_Save, "no");
  strcpy(Save_Tri, "no");
  strcpy(Load_Tri, "no");
  strcpy(Print, "no");
  strcpy(inthull, "no");
  strcpy(cddstyle, "no");
  strcpy(LRS, "no");

  for (i=1; i<argc-1; i++) {
    strcat(invocation,argv[i]);
    strcat(invocation," ");
    if (strncmp(argv[i],"vrep",3)==0) strcpy(Vrepresentation,"yes"); 
    else if (strncmp(argv[i],"int",3)==0) strcpy(interior,"yes");
    else if (strncmp(argv[i],"equ",3)==0) strcpy(equationsPresent,"yes");
    else if(strncmp(argv[i],"file",4)==0) strcpy(Memory_Save, "no");
    else if (strncmp(argv[i],"+", 1) ==0) strcpy(nonneg,"yes");
    else if (strncmp(argv[i],"printcones",3)==0) strcpy (Print, "yes");
    else if (strncmp(argv[i],"cdd",3)==0) strcpy (cddstyle, "yes");
    else if (strncmp(argv[i],"lrs",3)==0) strcpy (LRS, "yes");
    else if (strncmp(argv[i],"dil",3)==0) strcpy (dilation, "yes");
    else if (strncmp(argv[i],"rem",3)==0) {
      strcpy (removeFiles, "no");
      strcpy (Memory_Save, "no");
    }
    else if (strncmp(argv[i],"trisave",7)==0) {
      strcpy (Save_Tri, "yes");
      flags |= SAVE;
    }
    else if (strncmp(argv[i],"triload",7)==0) {
      strcpy (Load_Tri, "yes");
      flags |= LOAD;
    }
    else {
      cerr << "Unknown command/option " << argv[i] << endl;
      exit(1);
    }
  }
  if(printfile[0] == 'y') strcpy(Memory_Save, "no");
  if(printfile[0] == 'y') print_flag = 1;
  int dilation_const = 1;

  if(dilation[0] == 'y') dilation_const = atoi(argv[argc-2]);

  if(output_cone > 3) output_cone = 0;
  flags |= (output_cone << 1);
  if((cddstyle[0] == 'y') && (Vrepresentation[0] == 'y')){
    cerr << "Use not cdd style and v-representation." << endl;
    exit(2);
  }
  
  strcat(invocation,argv[argc-1]);
  strcat(invocation,"\n\n");
  cout << invocation;
  char costFile[127];
  strcpy(fileName,argv[argc-1]);
  //  cout << fileName << " " << costFile << endl;
  //strcpy (fileName,"stdin");
  
  /* Check input file. */
  if(Vrepresentation[0] == 'n'){
    if((cddstyle[0] == 'n') && (grobner[0] == 'n') && (maximum[0] == 'n')&& (minimize[0] == 'n')){
      CheckInputFile(fileName);
      CheckLength(fileName, equationsPresent);
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
  if((cddstyle[0] == 'n') && (Vrepresentation[0] == 'n')) CheckRed(fileName, equationsPresent, maximum, nonneg, interior, dilation, dilation_const); 

  dilation_const = 1;
  if((cddstyle[0] == 'n') && (grobner[0] == 'n'))
    readLatteProblem(fileName,&equations,&inequalities,equationsPresent,
		     &numOfVars, nonneg, dualApproach, grobner, maximum, 
		     cost,Vrepresentation);
//   if((equationsPresent[0] == 'n') && (interior[0] == 'y'))
//     Interior(inequalities);
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
  vector holdCost;
  holdCost = cost;
  //cout <<"Cost is: " << cost << endl;
  vec_RR holdcost_RR;
  holdcost_RR.SetLength(holdCost.length());
  for(i = 0; i < holdCost.length(); i++) conv(holdcost_RR[i], holdCost[i]);

  if((Vrepresentation[0] == 'y') && (equationsPresent[0] == 'y')){
    cerr<<"You cannot use vrep and equ at the same time." << endl;
    exit(4);
  }

  numOfVars--;

  numOfAllVars=numOfVars;
  mat_ZZ ProjU, ProjU2, AA;
  vec_ZZ bb;
  mat_ZZ AAA;

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
  matrix = matrixTmp;
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
    if (Vrepresentation[0] == 'n') {
    if(LRS[0] == 'n')
    tmpcones=computeVertexCones(fileName,matrix,numOfVars);
    else
      tmpcones=computeVertexConesViaLrs(fileName,matrix,numOfVars);
    {cones = tmpcones;
    cout << "\nThe polytope has " << lengthListCone(cones) << " vertices.\n";
    //system("rm numOfLatticePoints");
    cout << endl;}
  } 

    /* Compute triangulation or decomposition of each vertex cone. */

    cones=dualizeCones(cones,numOfVars);
    cones=decomposeCones(cones,numOfVars, flags, fileName);
    cones=dualizeBackCones(cones,numOfVars);

    /* Compute points in parallelepipeds, unless we already did using memsave version!  */
    
    cout << "Computing the points in the Parallelepiped of the unimodular Cones." << endl;
    computePointsInParallelepipeds(cones, numOfVars);

 if(Print[0] == 'y')
   printListCone(cones,numOfVars);

 cout << "Creating generating function.\n"; 
 //printListVector(templistVec, oldnumofvars); cout << ProjU << endl;
 if(equationsPresent[0] == 'y') {
   cones = ProjectUp2(cones, oldnumofvars, numOfVars, AA, bb);
   numOfVars = oldnumofvars;
 }
 createGeneratingFunctionAsMapleInput(fileName,cones,numOfVars);
 //printListCone(cones, numOfVars);
 cout << "Starting final computation.\n";
 cout << endl << "****  The number of lattice points is: " << Residue(cones,numOfVars) << "  ****" << endl << endl;


 if(rationalCone[0] == 'y') {
   strcpy(command, "mv ");
   strcat(command, "simplify.sum ");
   strcat(command, argv[argc - 1]);
   strcat(command, ".rat");
   system(command);
 }

 if(printfile[0] == 'y'){
   strcpy(command, "mv ");
   strcat(command, "func.rat ");
   strcat(command, argv[argc - 1]);
   strcat(command, ".rat");
   system(command);
 }
 if((removeFiles[0] == 'y')){
   
  strcpy(command,"rm ");
  strcat(command,fileName);
  strcat(command,".ext");
  system(command);
  
  strcpy(command,"rm ");
  strcat(command,fileName);
  strcat(command,".cdd");
  system(command); 
  
    strcpy(command,"rm ");
    strcat(command,fileName);
    strcat(command,".maple");
    system(command); 

  strcpy(command,"rm ");
  strcat(command,fileName);
  strcat(command,".ead");
  system(command); 
  
  if(cddstyle[0] == 'n'){
    strcpy(command,"rm ");
    strcat(command,fileName);
    system(command);
  }
 }
  

 cout << "Computation done. " << endl;
 cout << "Time: " << GetTime() << " sec\n\n";
 
 return(0);
}
/* ----------------------------------------------------------------- */



