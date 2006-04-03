/*************************************************************************

This is for computing a Gr\"obner basis of a toric ideal.

Author: Ruriko Yoshida
Date:   June 18th, 2003.
Update: July 1st, 2003.

I modified Raymond's code.

**************************************************************************/
#include "myheader.h"
#include "cone.h"
#include "print.h"
#include "ramon.h"
#include "ReadingFile.h"
#include "preprocess.h"
#include "Residue.h"
#include "vertices/cdd.h"
#include "dual.h"
#include "genFunction/piped.h"
#include "genFunction/maple.h"
#include "barvinok/dec.h"
#include "latte_system.h"

/* ---------------------------------------------------------------------- */
void CheckGrobner(char* filename, char * cdd){

  ifstream in(filename);
  ifstream in2(filename);

  int d, n, flag = 0;
  string tmp;
  if(cdd[0] == 'n'){
  getline(in2, tmp);

  if(tmp == "universal") flag = 1;
  if(flag == 1) getline(in, tmp);
  in >> d >> n;

  int counter = 0, int_tmp, correct;

  if(flag == 1) correct = d * n;
  else correct = d * n + n * n;
  while(in >> int_tmp){
    counter++;
  }
  if(counter != correct){
    cerr <<"Wrong input file!!" << endl;
    exit(1);}
  }
  else {
    ;

  }


}
/* ---------------------------------------------------------------------- */
void ReadGrobnerProblem(char* filename, mat_ZZ & A, mat_ZZ & W, ZZ & bignum, int & d, int & n, char * uni, char * cdd){
  if(cdd[0] == 'n'){
  ifstream in2(filename);
  ifstream in(filename);
  strcpy(uni, "no");
  string tmp;
  getline(in2, tmp);
  if(tmp == "universal") strcpy(uni, "yes");
  if(tmp == "universal") getline(in, tmp); 
  in >> d >> n;
  A.SetDims(d, n);
  W.SetDims(n, n);
 
  for(int i = 0; i < d; i++)
    for(int j = 0; j < n; j++)
      in >> A[i][j];
  if(uni[0] == 'n')
   for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++)
      in >> W[i][j];
  }
  else{
    string tmp, universal;
  ifstream in2(filename);
  ifstream in(filename);
  strcpy(uni, "no"); 
  while(!in2.eof()){
    getline(in2, tmp);
    if(tmp == "universal")  strcpy(uni, "yes");
  }
    while(tmp != "begin") getline(in, tmp);

  in >> d >> n >> tmp;
 
  A.SetDims(d, n);
  W.SetDims(n, n);

  for(int i = 0; i < d; i++)
    for(int j = 0; j < n; j++)
      in >> A[i][j];

  while(tmp != "grobner") getline(in, tmp);

  if(uni[0] == 'n')
   for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++)
      {in >> W[i][j];}
  }
  mat_ZZ B;
  B = transpose(A);
  for(int i = 0; i < n; i++)
     if(bignum < B[i] * B[i]) bignum = B[i] * B[i];
  bignum = (n + 1)*(n - d)* power(bignum, d/2);
} 

/* ------------------------------------------------------------------ */

void WriteGrobnerProblem(mat_ZZ A, mat_ZZ W, int index, char* filename, 
                         int d, int n){
 
  ofstream out(filename);
  int tmpindex = index;
   if(IsZero(W) != 0) index = 0;
  out << d + index << " " << n << endl;
  for(int i = 0; i < d; i++)
    {for(int j = 0; j < n; j++)
      { out << A[i][j] << " "; }
    out <<endl;
    }

  if(IsZero(W) == 0)
   for(int i = 0; i < index; i++)
    {for(int j = 0; j < n; j++) 
      out << W[i][j] << " ";
    out << endl;
    }
  index = tmpindex;
   out << d + index - 1 << " ";
   for(int i = 0; i < (d + index - 1); i++)
     out << i + 1 << " ";
   out << endl;

}
/* -------------------------------------------------------------- */


void readLatteGrobProblem(char *fileName, listVector **equations,
		      listVector **inequalities, 
		      char *equationsPresent,
                      int *numOfVars, char *nonneg, char* dual,
		      char* grobner, ZZ bignum) {
  int i,j,eq,ind,numOfVectors,numOfEquations;
  vec_ZZ indexEquations, tmpVector;
  listVector *basis, *endBasis, *tmp, *endEquations, *endInequalities;
  vec_ZZ b;

  /* Reads numOfVars, matrix A, and rhs b. */


  cout << "Reading problem.\n";

  //setbuf(stdout,0);

  ifstream in(fileName);
  if(!in){
    cerr << "Cannot open input file " << fileName << " in readLatteGrobProblem." << endl;
    exit(1);
  }

  if(grobner[0] == 'y') strcpy(equationsPresent, "yes");
  in >> numOfVectors;
  in >> (*numOfVars);

  if((dual[0] == 'y')&&(equationsPresent[0] == 'n')) *numOfVars = *numOfVars+1;

  int number = 0, oldNumOfVars = 0;

  if(grobner[0] == 'y'){
    ifstream in2(fileName);
    int dim, row;
    in2 >> dim >> row;
    int Matrix[dim][row];
    for(i = 0; i < dim; i++)
      for(j = 0; j < row; j++) in2 >> Matrix[i][j];
    in2 >> number;
  }

  if(grobner[0] == 'y'){
    oldNumOfVars = (*numOfVars);
    (*numOfVars) = 2 *(*numOfVars) + 1; 
  }

  mat_ZZ NonNeg;
  mat_ZZ A, B;
  A.SetDims(numOfVectors, *numOfVars);
  
  if((nonneg[0] == 'y')||(grobner[0] == 'y')) {
   NonNeg.SetDims(*numOfVars, *numOfVars);
   for(i = 0; i < *numOfVars - 1; i++)
     NonNeg[i][i+1] = 1;
    }

  if((dual[0] == 'n')||(equationsPresent[0] == 'y')){
  b=createVector(*numOfVars);
  if(grobner[0] == 'y'){
    for (j=1; j<=oldNumOfVars; j++){ in >> b[j]; b[j] = -b[j]; A[0][j] = b[j]; }
    for (j=oldNumOfVars + 1; j<(*numOfVars); j++){b[j] = - b[j - oldNumOfVars]; A[0][j] = b[j];}
    basis = createListVector(b);
    endBasis = basis;
    
    for (i=1; i<number; i++) {
      b=createVector(*numOfVars);
    for (j=1; j<=oldNumOfVars; j++){ in >> b[j]; b[j] = -b[j]; A[i][j] = b[j]; }
    for (j=oldNumOfVars + 1; j<(*numOfVars); j++){b[j] = - b[j - oldNumOfVars]; A[i][j] = b[j];}

      endBasis = updateBasis(createListVector(b), endBasis);
    }

    for (i=number; i<numOfVectors; i++) {
      b=createVector(*numOfVars); 
      b[0] = -1;
    for (j=1; j<=oldNumOfVars; j++){ in >> b[j]; }
    for (j=oldNumOfVars + 1; j<(*numOfVars); j++){b[j] = - b[j - oldNumOfVars];}

      endBasis = updateBasis(createListVector(b), endBasis);
    }


  }
  else{
    for (j=0; j<(*numOfVars); j++){ in >> b[j]; A[0][j] = b[j];}
    basis = createListVector(b);
    endBasis = basis;
    
    for (i=1; i<numOfVectors; i++) {
      b=createVector(*numOfVars);
      
      for (j=0; j<(*numOfVars); j++) {in >> b[j];
      A[i][j] = b[j];
      }
      endBasis = updateBasis(createListVector(b), endBasis);
    }

   
  }

  if(grobner[0] == 'y'){
  mat_ZZ UB;
  UB.SetDims(*numOfVars, *numOfVars);
  for(i = 0; i < *numOfVars; i++){
    UB[i][0] = bignum;
   
  }

  for(i = 0; i < (*numOfVars - 1); i++){
    UB[i][i + 1] = -1;
  }
   for(i = 0; i < *numOfVars - 1; i++){
     b=createVector(*numOfVars);
     for(j = 0; j < *numOfVars; j++) b[j] = UB[i][j];
     endBasis = updateBasis(createListVector(b), endBasis);
   }
   for(i = 0; i < *numOfVars - 1; i++){
     b=createVector(*numOfVars);
     for(j = 0; j < *numOfVars; j++) b[j] = NonNeg[i][j];
     endBasis = updateBasis(createListVector(b), endBasis);
   }
  }
 

  if(nonneg[0] == 'y'){
   for(i = 0; i < *numOfVars - 1; i++){
     b=createVector(*numOfVars);
     for(j = 0; j < *numOfVars; j++) b[j] = NonNeg[i][j];
     endBasis = updateBasis(createListVector(b), endBasis);
   }
  }
  }
  /*
  in >> numOfVectors;
  in >> (*numOfVars);
  if((dual[0] == 'y')&&(equationsPresent[0] == 'n')) *numOfVars = *numOfVars+1;
  mat_ZZ NonNeg;
  if(nonneg[0] == 'y') {
   NonNeg.SetDims(*numOfVars, *numOfVars);
   for(i = 0; i < *numOfVars - 1; i++)
     NonNeg[i][i+1] = 1;
    }
  if((dual[0] == 'n')||(equationsPresent[0] == 'y')){
  b=createVector(*numOfVars);
  for (j=0; j<(*numOfVars); j++) in >> b[j];
  basis = createListVector(b);
  endBasis = basis;

  for (i=1; i<numOfVectors; i++) {
    b=createVector(*numOfVars);
    for (j=0; j<(*numOfVars); j++) in >> b[j];
    endBasis = updateBasis(createListVector(b), endBasis);
  }

  if(nonneg[0] == 'y'){
   for(i = 0; i < *numOfVars - 1; i++){
     b=createVector(*numOfVars);
     for(j = 0; j < *numOfVars; j++) b[j] = NonNeg[i][j];
     endBasis = updateBasis(createListVector(b), endBasis);
   }
  }
  }
  */
  if((dual[0] == 'y')&&(equationsPresent[0] == 'n')){
 
  b=createVector(*numOfVars);
  ZZ hold;
  in >> hold;// cout << hold << endl;
  for (j=1; j<(*numOfVars)-1; j++) in >> b[j];
  b[*numOfVars-1] = hold;// cout << b << endl;
  basis = createListVector(b);
  endBasis = basis;

  for (i=1; i<numOfVectors; i++) {
    b=createVector(*numOfVars);
    in >> hold;//cout << hold << endl;
    for (j=1; j<(*numOfVars)-1; j++) in >> b[j];
  b[*numOfVars-1] = hold;
    endBasis = updateBasis(createListVector(b), endBasis);
  }

  }
  if (equationsPresent[0]=='n') {
    (*inequalities)=basis;
    (*equations)=0;
  } else {

    /* Read indices of equations and split basis into list of
       equations and inequalities. */

    in >> numOfEquations;
    indexEquations=createVector(numOfEquations);

    for (i=0; i<numOfEquations; i++) in >> indexEquations[i];
    cout << "\nEquation indices: ";
    printVector(indexEquations,numOfEquations);

    tmpVector=createVector(*numOfVars);
    createListVector(tmpVector);
    (*equations)=createListVector(tmpVector);

    (*inequalities)=createListVector(createVector(*numOfVars));
    endEquations=(*equations);
    endInequalities=(*inequalities);

    eq=0;
    ind=1;

    tmp=basis;
    while (tmp) {
      if (ind==indexEquations[eq]) {
	endEquations->rest=createListVector(tmp->first);
	endEquations=endEquations->rest;
	eq++;
	tmp=tmp->rest;
	if (eq==numOfEquations) {
	  endInequalities->rest=tmp;
	  tmp=0;
	}	
      } else {
	endInequalities->rest=createListVector(tmp->first);
	endInequalities=endInequalities->rest;
	tmp=tmp->rest;
      }
      ind++;
    }
    (*equations)=(*equations)->rest;
    (*inequalities)=(*inequalities)->rest;
  }

  cout << endl;
  cout << "Ax <= b, given as (b|-A):\n";
  cout << "=========================\n";
  printListVector(*inequalities,*numOfVars);
  
  cout << endl;

  cout << "Ax = b, given as (b|-A):\n";
  cout << "========================\n";
  printListVector(*equations,*numOfVars);

  cout << endl;

  return;
}
/******************************************************************/
void CheckFeasibility(listVector* basis, vec_ZZ infeas, int & flag) {
  flag = 0;
  if (basis==0) cout << "[]\n";
  while(basis) { 
    if(basis->first == infeas) flag = 1;
    basis = basis->rest; 
  }

/*  printf("\n"); */
  return ;
}


/******************************************************************/
listVector* Grobner(listVector *equations, 
		    listVector *inequalities, vec_ZZ **generators,
		    int *numOfVars, listVector **newVec,
		    int & oldnumofvars, int & flag) {
  int i,j,k,ind,ind2,indSol,lenOfMatrix,lenOfBasis,numOfIndependentRows,
    numOfRows,numOfVectors,newNumOfVars;
  ZZ det;
  vec_ZZ a,b,bas,rhs,A,U,H,sol,particularSolution;
  listVector *tmp, *tmp2, *basis, *endBasis, *newInequalities, 
    *endNewInequalities, *endtemplistVec, *tmplist, *templistVec;
  mat_ZZ M,unimodM;
  oldnumofvars = (*numOfVars);
  numOfRows=lengthListVector(equations);

  lenOfMatrix = (*numOfVars) * numOfRows;
  lenOfBasis  = (*numOfVars) * (*numOfVars);
  H=createVector(lenOfMatrix);
  rhs=createVector(numOfRows);

  tmp=equations;
  ind=0;
  ind2=0;
  while (tmp) {
    rhs[ind2]=(tmp->first)[0];
    for (i=0; i<(*numOfVars); i++) {
      H[ind]=(tmp->first)[i+1];
      ind++;
    }
    ind2++;
    tmp=tmp->rest;
  }
  H=negativeVector(H,lenOfMatrix);
  A=copyVector(H,lenOfMatrix);

  H=transpose(H,*numOfVars,numOfRows);

  bas=createVector(lenOfBasis);
  numOfIndependentRows=ihermite(&H,&bas,&rhs,numOfRows,*numOfVars);

  U=bas;
  ind=numOfIndependentRows*(*numOfVars);
  numOfVectors = (*numOfVars)-numOfIndependentRows;

  basis = createListVector(createVector(*numOfVars));
  endBasis = basis;
  
//   templistVec = createListVector(createVector(*numOfVars - 1));
//   endtemplistVec = templistVec;

  for (i=0; i<numOfVectors; i++) {
    b=createVector(*numOfVars);
    for (j=0; j<(*numOfVars); j++) b[j] = bas[ind+j];
    endBasis->rest=createListVector(b);
    endBasis = endBasis->rest;
    ind+=(*numOfVars);
  }

  basis=basis->rest;

  H=transpose(H,numOfRows,*numOfVars);
  U=transpose(U,*numOfVars,*numOfVars);

  /* Now basis contains the generators of the integer lattice.
     A contains the original matrix,
     U contains the unimodular transformation matrix, and
     H contains the Hermite normal form. 
     We have A.U = H. */

  sol=createVector(*numOfVars);

  for (i=0; i<(*numOfVars); i++) sol[i]=0;
//    cout << "sol:\n";
//    printVector(sol,*numOfVars);
//    cout << "numOfRows " << numOfRows << endl;
//    cout << "numOfVars " << *numOfVars << endl;

  indSol=0;
  for (i=0; i<numOfRows; i++) {
    if (H[(*numOfVars)*i+i]!=0) {
//    cout << "numOfRows " << numOfRows << endl;
//    cout << "numOfVars " << *numOfVars << endl;
//    cout << "lenOfMatrix " << lenOfMatrix << endl;
//    cout << "(i,i) " << (*numOfVars)*i+i << endl;
//        cout << i << " " << indSol << " " << rhs[i] << " " 
//  	   << H[(*numOfVars)*i+i] << endl;
      // cout << "Mmm..." << endl;
      sol[i]=rhs[i]/H[(*numOfVars)*i+i];
      indSol++;
      for (j=i+1; j<numOfRows; j++) {
	rhs[j]=rhs[j]-sol[i]*H[j*(*numOfVars)+i];
	//	H[j*(*numOfVars)+i]=0;
      }
    }
  }
  //   cout << "sol:\n";
  //   printVector(sol,*numOfVars);

  particularSolution=createVector(*numOfVars);
  for (i=0; i<(*numOfVars); i++) particularSolution[i]=0;

  for (i=0; i<(*numOfVars); i++) {
    particularSolution[i]=0;
    for (j=0; j<(*numOfVars); j++) {
      particularSolution[i]=particularSolution[i]+U[i*(*numOfVars)+j]*sol[j];
    }
  }
//    cout << "Particular solution:\n";
//    printVector(particularSolution,*numOfVars);
//    cout << "Basis:\n";
//    printListVector(basis,*numOfVars);

  newNumOfVars=lengthListVector(basis)+1;

  (*generators)=createArrayVector(newNumOfVars-1);
  tmp=basis;
  for (i=0;i<newNumOfVars-1;i++) {
    (*generators)[i]=copyVector(tmp->first,*numOfVars);
    tmp=tmp->rest;
  }

  M.SetDims(newNumOfVars,*numOfVars);
  for (i=0; i<newNumOfVars-1; i++) M[i]=(*generators)[i];
  LLL(det, M, unimodM);

  for (i=0; i<newNumOfVars-1; i++) (*generators)[i]=M[i];



  newInequalities=createListVector(createVector(*numOfVars));
  endNewInequalities=newInequalities;
  tmp=inequalities;

  for (i=0; i<lengthListVector(inequalities); i++) {
    a=tmp->first;
    b=createVector(newNumOfVars+1);
    b[0]=a[0];
    for (k=0; k<(*numOfVars); k++) {
      b[0]=b[0]+a[k+1]*particularSolution[k];
    }

    tmp2=basis;
    for (j=1; j<newNumOfVars; j++) {
      b[j]=0;
      for (k=0; k<(*numOfVars); k++) {
	b[j]=b[j]+a[k+1]*(tmp2->first)[k];
      }
      tmp2=tmp2->rest;
    }
    endNewInequalities->rest=createListVector(b);
    endNewInequalities=endNewInequalities->rest;
    tmp=tmp->rest;
  }
  newInequalities=newInequalities->rest;
/*    printf("OriginalInequalities:\n"); */
/*    printListVector(inequalities,(*numOfVars)+1); */
  cout << "New inequalities:\n";
  printListVector(newInequalities,newNumOfVars);
  vec_ZZ infeas;
  infeas.SetLength(newNumOfVars + 1);
  infeas[0] = -1;
  CheckFeasibility(newInequalities, infeas, flag);
    if(flag == 0){
  int check1 = 0;
  int check2 = 0;  
  templistVec = createListVector(createVector(newNumOfVars - 1));
  endtemplistVec = templistVec;
  int count = 0;
  // listVector *tmplist;
  tmplist = newInequalities;
  //  cout << numOfRows << endl;
  vec_ZZ HoldingUp;
  HoldingUp.SetLength(newNumOfVars - 1);

  // count = 1;
   for(count = 0; count < *numOfVars; ){
     //count++;
    // cout << tmplist -> first[0] << endl;
    if(tmplist -> first[0] == 0){
      // endtemplistVec -> rest = createListVector(createVector(newNumOfVars - 1));
      for(i = 0; i < (newNumOfVars - 1); i++){
	HoldingUp[i] = (tmplist -> first[i + 1]);
	//	endtemplistVec -> rest ->first[i] = (tmplist -> first[i + 1]);
	if(HoldingUp[i] == 0) check1++;
        if(HoldingUp[i] == 1) check2++;
      }
      // cout << HoldingUp << check1 << " " << check2 << " " << count<<endl;
      if((check1 != newNumOfVars - 2) || (check2 != 1)){
	//count++;
	endtemplistVec -> rest = createListVector(createVector(newNumOfVars - 1));
	endtemplistVec -> rest ->first = HoldingUp;
	//count++;
	endtemplistVec = endtemplistVec -> rest;
      }
      count++; 
      check1 = 0;
      check2 = 0;
    }
    tmplist = tmplist -> rest;
    
   }
   
  (*numOfVars)=newNumOfVars-1;
  templistVec = templistVec -> rest;
  *newVec = templistVec;
  //    printListVector(templistVec, newNumOfVars-1);
  return (newInequalities);
    }
    return 0;
}



/* -------------------------------------------------------------- */
void SolveGrobner(char * filename, char * nonneg, char * dualApproach,
	char * grobner, char * equationsPresent, char * cdd){
  
 ZZ bignum;
 mat_ZZ A, W;
 char uni[127];
 int numOfVars, oldnumofvars, d;
 ReadGrobnerProblem(filename, A, W, bignum, d, numOfVars, uni, cdd);
 char  outFile[127];
 strcpy(outFile, filename);
 strcat(outFile, ".gro");
 ofstream out(outFile);
 int CheckFeas = 0;
 char* File = "Gro.latte";
 unsigned int flags = 0;
 int nn = 0, n = numOfVars;
 if(IsZero(W) == 0) nn = numOfVars;
 else nn = 1;
 ZZ TotalNumLattice; 
 // for(int i = 0; i < nn; i++){
 //  numOfVars = n;
  vec_ZZ *generators;
  listVector *templistVec;
  listCone *cones, *tmp;
 listVector *matrix, *equations, *inequalities;
 for(int i = 0; i < nn; i++){
 numOfVars = n;
 WriteGrobnerProblem(A, W, i + 1,  File, d, numOfVars);

 readLatteGrobProblem(File,&equations,&inequalities,equationsPresent,
		      &numOfVars, nonneg, dualApproach, grobner, bignum);
 cones = 0;
  numOfVars--;
  //  numOfAllVars=numOfVars;
  CheckFeas = 0;  
  generators=createArrayVector(numOfVars);
  matrix=Grobner(equations,inequalities,&generators,&numOfVars, &templistVec, oldnumofvars, CheckFeas);
  if(CheckFeas == 0){
  cones=computeVertexConesViaLrs(File,matrix,numOfVars);

  cones=decomposeCones(cones,numOfVars, flags, File, 1, true,
		       BarvinokParameters::DualDecomposition);
  
  tmp=cones;

  while(tmp) 
    { 
      tmp->latticePoints=pointsInParallelepipedOfUnimodularCone(tmp->vertex,tmp->rays,numOfVars);    
      tmp=tmp->rest;  
    }
  
  cones = ProjectUp(cones, oldnumofvars, numOfVars, templistVec);

  numOfVars = oldnumofvars;

  createGeneratingFunctionAsMapleInputGrob(cones, numOfVars, out);

  TotalNumLattice += Residue(cones, numOfVars);

  system_with_error_check("rm -f Gro.latte*"); 
  }
 }
 //  system_with_error_check("rm Gro.latte*"); 
  out << "0;" << endl;
  cout << "*******************************************************" << endl;
  cout <<"\nThe total number of lattice points is: " << TotalNumLattice << "." << endl << endl;
  cout << "*******************************************************" << endl;
 
}

