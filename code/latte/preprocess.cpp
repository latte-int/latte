/* preprocess.cpp -- Preprocessing and projecting polytopes

   Copyright 2002 Raymond Hemmecke, Ruriko Yoshida

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

#include "cone.h"
#include "print.h"
#include "ramon.h"

/* ----------------------------------------------------------------- */
vec_ZZ transpose(vec_ZZ mat, int numOfVars, int numOfRows){
  int i,j,k,lenOfMatrix;
  vec_ZZ transposedMat;

/*    printf("Transposing matrix\n"); */

  lenOfMatrix = numOfVars*numOfRows;
  transposedMat=createVector(lenOfMatrix);
  k=0;
  for (j=0; j<numOfVars; j++)
    for (i=0; i<numOfRows; i++)
      transposedMat[k++] = mat[i*numOfVars+j];
  return(transposedMat);
}
/* ----------------------------------------------------------------- */
/* Code originally from TiGERS. This is an adapted version from MLP. */

int ihermite(vec_ZZ *S, vec_ZZ *U, vec_ZZ* rhs, int m, int n){
  int i,j,k,done,sign,c,mc,mn,crk;
  ZZ mv,t;

  mv=0;
  t=0;
  c=1;
  crk=0;

  cerr << "Computing hermitean normal form.\n";

  if (m>n) mn=n; else mn=m;

  /* Initialize U to nxn identity */
  for(i=1; i<=n; i++) {
    for(j=1; j<=n; j++) {
      if (i==j) (*U)[(i-1)*n+j-1]=1;
      else (*U)[(i-1)*n+j-1]=0;
    }     
  }
 
  while (c<=mn) {
    /* find minimum entry in col c */
    mv=(*S)[(c-1)*m+c-1];
    if (mv<0) mv*=-1;
    mc=c;
    for(i=c+1; i<=n; i++) {
      t=(*S)[(i-1)*m+c-1];
      if (t<0) t*=-1; 
      if(mv==0 || (mv>t && t!=0)) {
           mv=t;
           mc=i;
      }
    }

    /* if nescesary pivot to put min in row c and multiply by+-1
       to ensure diagonal entry is positive */
    if (mc!=c || (*S)[(mc-1)*m+c-1]<0) {
      if ((*S)[(mc-1)*m+c-1]<0) sign=-1;
      else sign=+1;
      for(j=c; j<=m; j++) {
        t=(*S)[(mc-1)*m+j-1];
        (*S)[(mc-1)*m+j-1]=(*S)[(c-1)*m+j-1];
        (*S)[(c-1)*m+j-1]=sign*t;
      }
      for(j=1; j<=n; j++) {
        t=(*U)[(mc-1)*n+j-1];
        (*U)[(mc-1)*n+j-1]=(*U)[(c-1)*n+j-1];
        (*U)[(c-1)*n+j-1]=sign*t;
      }
    }

    /* if column is not zero do a reduction step */

    if ((*S)[(c-1)*m+c-1]==0) {
      done=0;
      for(i=c; i<=m; i++) 
	for(j=c; j<=n; j++) {
	  if ((done==0) && ((*S)[(j-1)*m+i-1]!=0)) {
	    done=1;
	    for(k=1; k<=n; k++) {
	      t=(*S)[(k-1)*m+i-1];
	      (*S)[(k-1)*m+i-1]=(*S)[(k-1)*m+c-1];
	      (*S)[(k-1)*m+c-1]=t;
	    }
	    t=(*rhs)[i-1];
	    (*rhs)[i-1]=(*rhs)[c-1];
	    (*rhs)[c-1]=t;
	  }
      }
      if (done==0) {return (crk);}
      done=0;
    } else {
      done=1;
      crk=c;
      for(i=c+1; i<=n; i++) {
        t=(*S)[(i-1)*m+c-1]/(*S)[(c-1)*m+c-1];
        for(j=c; j<=m; j++) {
           (*S)[(i-1)*m+j-1]-=t*(*S)[(c-1)*m+j-1];
        }
        for(j=1; j<=n; j++) {
           (*U)[(i-1)*n+j-1]-=t*(*U)[(c-1)*n+j-1];
        }
        if ((*S)[(i-1)*m+c-1]!=0) done=0;
      }
    }
    /* if all entrees of col c below row c are zero go to next col */
    if (done==1) c++;
  }

  return (crk);
 }

/* ----------------------------------------------------------------- */
void checkListVector(listVector* basis, int numOfVars) {
  if (basis==NULL){
    cerr << "\n\n**** Total number of lattice points is: 0 ****\n" << endl;
    ofstream out("numOfLatticePoints");
    out << 0 << endl;
    exit(0);
  } 

  ZZ counter, RHS;
  while(basis) {
    counter = 0;
    RHS = basis -> first[0];
    for(int i = 1; i < numOfVars; i++) counter += abs(basis -> first[i]);
    basis = basis->rest;
    if((IsZero(counter) == 1) &&(RHS < 0)){
      cerr << "\n\n**** Total number of lattice points is: 0 **** \n" << endl;
      ofstream out("numOfLatticePoints");
      out << 0 << endl;
      exit(0);
    }
    else if((IsZero(counter) == 1) &&(RHS > 0))  removeListVector(basis);
  }
  /*  printf("\n"); */
  return ;
}

void dilateListVector(listVector* basis, int numOfVars, int dil){
  ZZ dil_ZZ;
  conv(dil_ZZ, dil);

  while(basis) {
    basis -> first[0] = dil_ZZ * basis -> first[0];
    basis = basis->rest;
  }
  return ;
}

void Interior(listVector* basis){
  while(basis){
    basis->first[0]--;
    basis = basis->rest;
  }
  
}

/* ----------------------------------------------------------------- */
listVector* preprocessProblem(listVector *equations, 
			      listVector *inequalities, vec_ZZ **generators,
			      int *numOfVars, vec_ZZ & cost, mat_ZZ & ProjU, char* interior, int dil) {
  int i,j,k,ind,ind2,indSol,lenOfMatrix,lenOfBasis,numOfIndependentRows,
    numOfRows,numOfVectors,newNumOfVars;
  ZZ det;
  vec_ZZ a,b,bas,rhs,A,U,H,sol,particularSolution;
  listVector *tmp, *tmp2, *basis, *endBasis, *newInequalities, 
    *endNewInequalities;
  mat_ZZ M,unimodM, Solve;
  //  cerr << *numOfVars << lengthListVector(equations) << endl;
  if(inequalities == 0){ 
    if(lengthListVector(equations) == *numOfVars){
      tmp = equations; i = 0;
      sol.SetLength(*numOfVars);
      Solve.SetDims(*numOfVars, *numOfVars);
      for(i=0; i < *numOfVars; i++){
	sol[i] = tmp -> first[0]; 
	for(j = 1; j < *numOfVars + 1; j++) Solve[i][j - 1] = -tmp -> first[j];
	tmp = tmp -> rest;
      }
      // cerr << Solve << sol << endl;
      vec_ZZ x;
      x.SetLength(*numOfVars);
      ZZ d, sum, sum2;
      mat_ZZ Inv;
      inv(d, Inv, Solve); x = Inv * sol;
      for(i = 0; i < *numOfVars; i++) sum += x[i];
      for(i = 0; i < *numOfVars; i++) sum2 += x[i]/d;
      //cerr << d << " " <<sum << " " <<d*sum2 << x << endl; exit(0);
      if(sum == d*sum2){
	ofstream OUT("numOfLatticePoints");
        cerr << "The number of lattice points is 1." << endl;
	OUT << 1 << endl;
	exit(0);
      }else{
	cerr << "The number of lattice points is 0." << endl;
	ofstream OUT("numOfLatticePoints");
	OUT << 0 << endl;
	exit(0);}
    }
    else{
      cerr << "The polytope is not bounded." << endl;
      exit(1);
    }
  }
  numOfRows=lengthListVector(equations);
  
  lenOfMatrix = (*numOfVars) * numOfRows;
  lenOfBasis  = (*numOfVars) * (*numOfVars);
  H=createVector(lenOfMatrix);
  rhs=createVector(numOfRows);
  
  tmp=equations;
  ind=0;
  ind2=0;
  int flag = 0;
  if(cost.length() != 0) flag = 1;

  while (tmp) {
    rhs[ind2]=(tmp->first)[0];
    for (i=0; i<(*numOfVars); i++) {
      H[ind]=(tmp->first)[i+1];
      ind++;
    }
    ind2++;
    tmp=tmp->rest;
  }
  H=-H;
  A=H;

  H=transpose(H,*numOfVars,numOfRows);

  bas=createVector(lenOfBasis);
  numOfIndependentRows=ihermite(&H,&bas,&rhs,numOfRows,*numOfVars);

  U=bas;
  ind=numOfIndependentRows*(*numOfVars);
  numOfVectors = (*numOfVars)-numOfIndependentRows;

  basis = createListVector(createVector(*numOfVars));
  endBasis = basis;

  for (i=0; i<numOfVectors; i++) {
    b=createVector(*numOfVars);
    for (j=0; j<(*numOfVars); j++) b[j] = bas[ind+j];
    endBasis->rest=createListVector(b);
    endBasis = endBasis->rest;
    ind+=(*numOfVars);
  }

  {
    // Drop the dummy head.
    listVector *b = basis->rest;
    delete basis;
    basis=b;
  }

  H=transpose(H,numOfRows,*numOfVars);
  U=transpose(U,*numOfVars,*numOfVars);

  /* Now basis contains the generators of the integer lattice.
     A contains the original matrix,
     U contains the unimodular transformation matrix, and
     H contains the Hermite normal form. 
     We have A.U = H. */
  mat_ZZ UU, HH, HHH;
  UU.SetDims(*numOfVars, *numOfVars);
  HH.SetDims(numOfRows, *numOfVars);
  HHH.SetDims(numOfRows, numOfRows);

  int counter = 0;
  for(i = 0; i < *numOfVars; i++){
    for(j = 0; j < *numOfVars; j++){ UU[i][j] = U[counter];
    counter++;
    }
  }

  counter = 0;
  for(i = 0; i < numOfRows; i++){
    for(j = 0; j < *numOfVars; j++){ HH[i][j] = H[counter];
    counter++;
    }
  }

  for(i = 0; i < numOfRows; i++){
    for(j = 0; j < numOfRows; j++){ HHH[i][j] = HH[i][j];
    }
  }

  // cerr << HH << UU << endl;
  //  cerr << rhs << endl;
  ZZ DD, zeros;
  mat_ZZ invHH;
  inv(DD, invHH, HHH);
  //  cerr <<  invHH <<endl << rhs << endl << DD << endl;
  vec_ZZ sol2;
  sol=createVector(*numOfVars);
  sol2=createVector(*numOfVars);
  if(DD != 0)
    sol2 = (invHH * rhs);

  for (i=0; i<(*numOfVars); i++) sol[i]=0;
//    cerr << "sol:\n";
//    printVector(sol,*numOfVars);
//    cerr << "numOfRows " << numOfRows << endl;
//    cerr << "numOfVars " << *numOfVars << endl;

  indSol=0;
  for (i=0; i<numOfRows; i++) {
    if (H[(*numOfVars)*i+i]!=0) {
//    cerr << "numOfRows " << numOfRows << endl;
//    cerr << "numOfVars " << *numOfVars << endl;
//    cerr << "lenOfMatrix " << lenOfMatrix << endl;
//    cerr << "(i,i) " << (*numOfVars)*i+i << endl;
//        cerr << i << " " << indSol << " " << rhs[i] << " " 
//  	   << H[(*numOfVars)*i+i] << endl;
      // cerr << "Mmm..." << endl;
      sol[i]=rhs[i]/H[(*numOfVars)*i+i];
      indSol++;
      for (j=i+1; j<numOfRows; j++) {
	rhs[j]=rhs[j]-sol[i]*H[j*(*numOfVars)+i];
	//	H[j*(*numOfVars)+i]=0;
      }
    }
  }
    cerr << "sol:\n";
    //    int flag_sol = 0;
    printVectorToFile(cerr,sol,*numOfVars);
    if(DD != 0){
      for(i = 0; i < *numOfVars; i++){
	zeros = abs(sol2[i] - sol[i]*DD);
	if(zeros != 0) { 
	cerr << "Integrally empty polytope." << endl;
	cerr << "\n\n**** Total number of lattice points: 0 ****" << endl << endl;
	ofstream OutZero("numOfLatticePoints");
	OutZero << 0 << endl;
	exit(0);}
	zeros = 0;
      } 
    }//cerr << sol << endl;
  particularSolution=createVector(*numOfVars);
  for (i=0; i<(*numOfVars); i++) particularSolution[i]=0;

  for (i=0; i<(*numOfVars); i++) {
    particularSolution[i]=0;
    for (j=0; j<(*numOfVars); j++) {
      particularSolution[i]=particularSolution[i]+U[i*(*numOfVars)+j]*sol[j];
    }
  }
//    cerr << "Particular solution:\n";
//    printVector(particularSolution,*numOfVars);
//    cerr << "Basis:\n";
//    printListVector(basis,*numOfVars);

  newNumOfVars=lengthListVector(basis)+1;

  (*generators)=createArrayVector(newNumOfVars-1);
  tmp=basis;
  for (i=0;i<newNumOfVars-1;i++) {
    (*generators)[i]=tmp->first;
    tmp=tmp->rest;
  }

  M.SetDims(newNumOfVars,*numOfVars);
  for (i=0; i<newNumOfVars-1; i++) M[i]=(*generators)[i];
  LLL(det, M, unimodM);

  for (i=0; i<newNumOfVars-1; i++) (*generators)[i]=M[i];



  newInequalities=createListVector(createVector(*numOfVars));
  endNewInequalities=newInequalities;
  tmp=inequalities;
  vec_ZZ tmpcost;
  if((flag == 1)||(flag == 0))
    { 
      ProjU.SetDims((*numOfVars), newNumOfVars);
      mat_ZZ tmpProjU, Proj;
      tmpProjU.SetDims((*numOfVars), (*numOfVars) + 1);
      Proj.SetDims((*numOfVars), newNumOfVars + 1);

      for(i = 0; i < *numOfVars; i++) tmpProjU[i][i + 1] = 1;
      tmpcost = cost;
      cost.kill();
      cost.SetLength(newNumOfVars - 1);

      for(i = 0; i <*numOfVars; i++){
	Proj[i][0]=tmpProjU[i][0];
	for (k=0; k<(*numOfVars); k++) {
	  Proj[i][0]=Proj[i][0]+tmpProjU[i][k + 1]*particularSolution[k];
	}
	
	tmp2=basis;
	for (j=1; j<newNumOfVars; j++) {
	  ProjU[i][j]=0;
	  for (k=0; k<(*numOfVars); k++) {
	    Proj[i][j]=Proj[i][j]+tmpProjU[i][k + 1]*(tmp2->first)[k];
	  }
	  tmp2=tmp2->rest;
	}
	for(int m = 0; m < newNumOfVars; m++) ProjU[i][m] = Proj[i][m];
      } //cerr << ProjU << endl;

      //  cerr << tmpcost << endl;
      /*   cost[0]=tmpcost[0];
	   for (k=0; k<(*numOfVars); k++) {
	   cost[0]=cost[0]+tmpcost[k]*particularSolution[k];
	   }*/
      
      tmp2=basis;
     if(flag == 1){ for (j=0; j<newNumOfVars - 1; j++) {
	cost[j]=0;
	for (k=0; k<(*numOfVars); k++) {
	  cost[j]=cost[j]+tmpcost[k]*(tmp2->first)[k];
	  // cerr << (tmp2->first)[k] << " ";
	}
	tmp2=tmp2->rest;
     }}
    }
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
    // if(interior[0] == 'y') b[0]--;
    if(interior[0] == 'y') ;
    ZZ check_Tri;
    for(int i = 1; i < newNumOfVars; i++)
      check_Tri += abs(b[i]);
    ZZ tmp_dil;
    conv(tmp_dil, dil); 
    b[0] = tmp_dil * b[0];
    if(check_Tri != 0){
      endNewInequalities->rest=createListVector(b);
      endNewInequalities=endNewInequalities->rest;
    } 
    if((check_Tri == 0) && (b[0] < 0))
      {
	cerr << "this polytope is empty!" << endl;
	ofstream OutPut("numOfLatticePoints");
	OutPut << 0 << endl;
	exit(0);
      }
    tmp=tmp->rest;
  }
  {
    // Drop the dummy head.
    listVector *t = newInequalities->rest;
    delete newInequalities;
    newInequalities = t;
  }
  freeListVector(basis);
/*    printf("OriginalInequalities:\n"); */
/*    printListVector(inequalities,(*numOfVars)+1); */
  cerr << "New inequalities:\n";
  printListVectorToFile(cerr, newInequalities,newNumOfVars + 1);
  checkListVector(newInequalities, newNumOfVars);
  if(newInequalities == NULL){

    cerr << "\n\n****  Total number of lattice points is: 0 ****\n" << endl;
    ofstream out("numOfLatticePoints");
    out << 0 << endl;
    exit(0);
  }
  (*numOfVars)=newNumOfVars-1;
  return (newInequalities);
}

/******************************************************************/



listVector* transformZZMatToListVector(mat_ZZ A, int numOfVectors,
						int numOfVars) {
  int i;
  vec_ZZ v;
  listVector *L, *endL;

  v=createVector(numOfVars);
  L=createListVector(v);
  endL=L;

  for (i=0; i<numOfVectors; i++) {
    v=A[i];
    endL->rest = createListVector(v);
    endL = endL->rest;
  }

  listVector *result = L->rest;
  delete L; // deletes dummy head
  return result;
}
 listVector* TransformToDualCone(listVector* matrix, int& numOfVars){
 int numOfCons = lengthListVector(matrix);
 mat_ZZ ConeMatrix, tmpMatrix;
 vec_ZZ tmp;
 { tmp = matrix->first; //cerr << tmp <<" " << numOfVars << endl;
    matrix = matrix -> rest;}
 int number = tmp.length();
 tmpMatrix.SetDims(numOfCons, number);
 ConeMatrix.SetDims(numOfCons, number + 1);
 //listVector* rays;
 //rays = matrix;
 tmpMatrix[0] = tmp;
 for(int i = 1; i < numOfCons; i++)
    { tmpMatrix[i] = matrix->first;// cerr << tmp <<" " << numOfVars << endl;
    matrix = matrix -> rest;}     // exit(3);
 for(int i = 0; i < numOfCons; i++)
     ConeMatrix[i][number-1] = tmpMatrix[i][0];

 for(int i = 0; i < numOfCons; i++)
    for(int j = 1; j < number-1; j++){
       ConeMatrix[i][j] = tmpMatrix[i][j];
    }   
  cerr << endl << "After projecting up for a dual cone:" << endl;
  cerr <<"===================================" << endl;
  for(int i = 0; i < numOfCons; i++){
     cerr <<"[";
     for(int j = 0; j < number - 1; j++)
       cerr << ConeMatrix[i][j] << " ";
     cerr << ConeMatrix[i][number - 1] << "]" << endl;
     }
   cerr << "===================================" << endl;
  numOfVars = numOfVars + 1;
  return(transformZZMatToListVector(ConeMatrix, numOfCons, number + 1));
}
/* ----------------------------------------------------------------- */

vec_ZZ ProjectingUp(mat_ZZ ProjU, vec_ZZ cost, int numOfVars){
  // cerr << ProjU << endl;
  int numOfrows = ProjU.NumRows(), numOfcolms = ProjU.NumCols();
  vec_ZZ answer; 
  if(IsZero(ProjU))
    answer = cost;
  else{
    answer.SetLength(numOfrows);
    for(int i = 0; i < numOfrows; i++){
      answer[i] = ProjU[i][0];
      for(int j = 1; j < numOfcolms; j++){
      answer[i] += cost[j - 1]*ProjU[i][j];
      }
      
    }
  }// cerr << answer << endl;
  return answer;
}

/* ----------------------------------------------------------------- */

vec_RR ProjectingUpRR(mat_RR ProjU, vec_RR cost, int numOfVars){
  // cerr << ProjU << endl;
  int numOfrows = ProjU.NumRows(), numOfcolms = ProjU.NumCols();
  vec_RR answer; 
  if(IsZero(ProjU))
    answer = cost;
  else{
    answer.SetLength(numOfrows);
    for(int i = 0; i < numOfrows; i++){
      answer[i] = ProjU[i][0];
      for(int j = 1; j < numOfcolms; j++){
      answer[i] += cost[j - 1]*ProjU[i][j];

      }
      if(answer[i] < power2_RR(-10)) answer[i] = 0;      
    }
  }// cerr << answer << endl;
  return answer;
}

/* ----------------------------------------------------------------- */
