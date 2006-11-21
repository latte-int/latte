/* Residue.cpp -- Polynomial substitution and residue calculations

   Copyright 2002-2004 Jesus A. De Loera, David Haws, Raymond
      Hemmecke, Peter Huggins, Jeremy Tauzer, Ruriko Yoshida
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

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <time.h>
#include <list>
#include <vector>
#include "config.h"
#include "cone.h"
#include "ramon.h"
using namespace std;

/* From here, Rudy edited ----------------------------------------- */
/* ----------------------------------------------------------------- */

ZZ
Residue(listCone* cones, int numOfVars)
{
  int numOfTerms;
//  char outFileName[127];
  listVector *tmp;
  listCone *C, * cones1;
  int dim, noGsPerC,noCones; //noGsPerC is number of generators per cone
  clock_t t,sc=0,sc2=0;

 // strcpy(outFileName,fileName);
 // strcat(outFileName,".residue");

 /* ofstream out(outFileName);
  if (!out) {
    printf("Error opening output file for writing in printResidueFile!");
    exit(1);
  }
  if (cones==0) out << "No cones in list.\n";     */

  numOfTerms=0;

  C=cones;
  while (C) {
    numOfTerms=numOfTerms+lengthListVector(C->latticePoints);
    C=C->rest;
  }


 /* out << numOfVars << " " << lengthListVector(cones->rays) << " " <<
    numOfTerms << "\n\n";  */

  dim=numOfVars;
  noGsPerC=lengthListVector(cones->rays);
  noCones=numOfTerms;
  int i,j; // index or loop vars
  long int k, m;//n=0,p; // extra vars to use as needed
  long int p, n;
  int E[noCones];  // E is the vector of epsilons, each 1 or -1
  long int totalNoGs=noGsPerC*noCones; //total no. of generators,ie,rowdim of B
  list<Integer> A[noCones];  // A is the numerator vectors
  // long int B[totalNoGs][dim];  // B is the denominator vectors
  //  cout<<"tNG: "<<totalNoGs<<endl;

  class denom {
  public:
    Integer *D;
    denom * next;
    denom(int noGPC){
      D = new Integer [noGPC];
    }
  };

  Integer tmp_A;
  denom * B=new denom(noGsPerC*dim);

  denom * Bitr=B;
  listVector* basis, *listtmp1, *listtmp2;
  listCone *listtmp3;
  cones1 = cones;
  i = 0;
  while (cones1) {
    tmp=cones1->latticePoints;
    /* The code is not prepared for non-unimodular cones (segfault) --
       check this out later! --mkoeppe */
    assert(tmp != NULL && tmp->rest == NULL);
    while (tmp) {
      E[i] = cones1->coefficient;
     // printVectorToFileWithoutBrackets(out,tmp->first,numOfVars);
    for (j=0; j<(numOfVars); j++) {
    tmp_A = tmp->first[j];  A[i].push_back(tmp_A);
    //cout << tmp_A << " ";
    
    }
    //cout << endl;
     // printListVectorToFileWithoutBrackets(out,cones->rays,numOfVars);
     basis = cones1->rays;
     while(basis) {
    //printVectorToFileWithoutBrackets(out,basis->first,numOfVars);
    for (j=0; j<noGsPerC; j++) {
    for(k = 0; k <dim; k++){
    Bitr->D[j*dim+k] = basis->first[k];
      }
    listtmp1 = basis;
    basis = basis->rest;
    delete listtmp1;
    }
    Bitr->next=new denom(noGsPerC*dim);
    Bitr=Bitr->next;
  // i++;
  }
    //  out << endl;
    listtmp2 = tmp;
      tmp=tmp->rest; i++;
      delete listtmp2;
    }
    listtmp3 = cones1;
    cones1 = cones1->rest;
    delete listtmp3;
  }
 // out << endl;
  i = 0;
 /* denom * Bitr=B;
  for(i=0;i<noCones;i++) {
    input>>E[i];
    for(j=0;j<dim;j++) {input>>tmp_A; A[i].push_back(tmp_A);}
    for(j=0;j<noGsPerC;j++)
      {
	for(k=0;k<dim;k++) input>>Bitr->D[j*dim+k];
      }
    Bitr->next=new denom(noGsPerC*dim);
    Bitr=Bitr->next;
  }
  input.close();  */

  
    /* Bitr=B;
     Bitr=Bitr->next;
     cout<<"B[1]: ";
     for(i=0;i<noGsPerC;i++) {cout<<endl;
     for(j=0;j<dim;j++) cout<<Bitr->D[i*dim+j]<<" ";
     }
     return ; */    


  // cout<<".done"<<endl<<"Number of cones: "<<noCones<<endl;cin.get();
//        <<"Number of generators in total: "<<totalNoGs<<endl
//        <<endl<<"Now start calculations and stopwatch."<<endl;
  t=clock();
  //  cout<<"  Clock at start reads "<<t<<"."<<endl;

  //------------------------------------------------------------------------------
  //---FIND LAMBDA AND DENOMINATOR EXPONENTS: We want to make substitution
  //---x_i -> t^lambda[i] to leave everything in terms of just 1 variable (t).  We
  //---can only do this if no term in the denominator becomes 0.  Since a factor
  //---in the denominator of a Brion series is of the form 1-x^(row j of B), where
  //---the exponent is multivariate, we must ensure that lambda dotted with any
  //---row of B does not yield zero.
  //---  To search for a feasible lambda, first I run through lambdas with entries
  //---in {-1,0,1,2} by increasing lexicographical order.  When I reach a row of B
  //---(row tracked by var. 'halt') which yields 0 when dotted with lambda, I need
  //---to change an entry of lambda corresponding to a nonzero entry in the halting
  //---vector, which may let me skip some lambdas.
  //------------------------------------------------------------------------------

  //  cout<<"Getting lambda..."<<endl;
  //cout<<"k is "<<k;cin.get();

  //VARS:
  //Integer k;
  Integer tmp_lambda;
  Integer lambda[dim];
  Integer dlambda[dim]; // dlambda tracks change in 2 successive test-lambdas
  //for(i=0;i<dim;i++) lambda[i]=0;  // lambda starts at 0
  vector<Integer> dotProducts(totalNoGs); // ith entry tracks lambda dot row i of B
  //  dotProducts and dlambda used to try to improve calculational efficiency
  long int halt, haltCone;
  halt = 0;
  haltCone=-1; // will track the index where the dot product is 0
  k=0; //cout<<"k is "<<k;// loop control var
  // Also, n tracks the lowest index where lambda changed (dlambda[n] not zero).
 // cout << totalNoGs << endl;
  //--------LOOP 1: try up to 5000 lambdas with entries in {-1,0,1,2}
  //cout<<"k is "<<k;cin.get();
 while(k<5000) {
   //cout << k << " ";cin.get();
    //--FIRST check if lambda works.
    Bitr=B;
    for(i=0;i<haltCone;i++) {
      for(j=0;j<noGsPerC;j++) {
	m=i*noGsPerC+j;
	for(p=n;p<dim;p++) dotProducts[m]+= dlambda[p] * Bitr->D[j*dim+p];
	if(dotProducts[m]==0) {haltCone=i; halt=j; j=noGsPerC; i=noCones+2;}
      }
      if(i<noCones) Bitr=Bitr->next;
    }
    for(;i<noCones;i++) { // keep checking if not yet halted
      for(j=0;j<noGsPerC;j++) {
	m=i*noGsPerC+j;
	dotProducts[m]=0;
	for(p=0;p<dim;p++) dotProducts[m] += lambda[p] * Bitr->D[j*dim+p];
	if(dotProducts[m]==0) {haltCone=i; halt=j; j=noGsPerC; i=noCones+2;}
      }
      if(i<noCones) Bitr=Bitr->next;
    }
    if(i==noCones) { // found lambda
      k=200000; }
    else {

      //--SECOND find next lambda.
      m= (Bitr->D[halt*dim+dim-1] ==0 ? 0 : 1);
	n=dim-1;
	while(n>=1 && (lambda[n]==2 || m==0)) {
	  n--;
	  if(Bitr->D[halt*dim+n] !=0) m++;
	}
	if(lambda[n]==2 || k==4999) {// failed: no lambda with -1,0,1,2 entries found
	  k=60000; }
	else { // now I set n and up entries of lambda and dlambda
	    lambda[n] ++;
	    dlambda[n]=1;
	    for(j=n+1;j<dim;j++) {
	      dlambda[j]=-1-lambda[j];
	      lambda[j]=-1;
	    }
	  } // end if/else
      } // end if/else
    k++;
  } // end while
 // cout << totalNoGs << endl;

 // for(j = 0; j < totalNoGs; j++) cout << dotProducts[j] << " ";
  // cout << endl;

  //--IN CASE no lambda found yet...

// srand(time(0));  //This algorithm not used... it picks a random lambda.
// while(k<-10300) {
//   k++;
//   n=(k-9960)/10;
//   for(j=0;j<dim;j++) lambda[j]=rand()% n - n/2;
//   for(i=0;i<totalNoGs;i++) {
//     dotProducts[i]=0;
//     for(j=0;j<dim;j++) dotProducts[i]+=lambda[j]*B[i][j];
//     if(dotProducts[i]==0) i=totalNoGs+2;
//   }
//   if(i==totalNoGs) k=200000;
// } cout<<"k: "<<k<<endl<<"clk: "<<clock()<<endl;

 //-------LOOP 2: in case no lambda found yet, try bigger lambda entries...
  int d=0; // d runs through odds... it gets added to a single entry of lambda
       // each time through the loop, to fix the problem generator.
  //if(k<15000) {for(i=0;i<dim;i++) lambda[i]=-1;} if running random algorithm
  //halt=0; if running random algorithm
  while(k<150000) {
    d++;
    //FIRST find new lambda.
   /* n=0;
    while(Bitr->D[halt*dim+n]==0) n++;
    d+=2;
    lambda[n]+=d; */
    for(int s = 0; s < dim; s++){
       lambda[s] = rand() % 1000;
       if((rand() % 2) == 1) lambda[s] = - lambda[s];
       else ;
       }
    //SECOND check if lambda works
    Bitr=B;
 /*   for(i=0;i<haltCone;i++) {
      for(j=0;j<noGsPerC;j++) {
	m=i*noGsPerC+j;
	dotProducts[m]+=d*Bitr->D[j*dim+n];
	if (dotProducts[m]==0) {haltCone=i; halt=j; j=noGsPerC; i=noCones+2;}
      }
      if(i<noCones) Bitr=Bitr->next;
    }   */
    for(i = 0;i<noCones;i++) {
      for(j=0;j<noGsPerC;j++) {
	m=i*noGsPerC+j;
	dotProducts[m]=0;
	for(p=0;p<dim;p++) dotProducts[m]+=lambda[p]*Bitr->D[j*dim+p];
	if(dotProducts[m]==0) {haltCone=i; halt=j; j=noGsPerC;
   i=noCones+2;}
      }
      if(i<noCones) Bitr=Bitr->next;
    }
    if(i==noCones) k=200000;
  }

  // cout<<"lambda: "; cin.get();
  // for(i=0;i<dim;i++) cout<<lambda[i]<<" ";
  // cout<<endl<<"  Clock at checkpoint 2: "<<clock()<<endl;
  //cin.get();
  //  cout<<"denominator: ";
 //   for(i=0;i<noGsPerC*noCones;i++) cout<<dotProducts[i]<<" ";
 //   cout<<i <<endl;
  //  return 0;
  //----------------------------------------------------------------------------
  //---CALCULATE NUMERATOR EXPONENTS of t under substitution x[i]->t^lambda[i].
  //---After getting initial exponents, I translate so as to try to minimize
  //---numerator exponents (and promote efficientcy).
  //----------------------------------------------------------------------------


  /****************************************************************************
    Rudy's comment:
    This is a process to simplify a univariable rational function.
    For example, if I have a rational function f(K(v)) = t/(1-t^-1)(1-t^-2),
    then, the following loop simplifies as f(K(v)) = t^4 /(1-t)(1-t^2).

  *****************************************************************************/
  //  cout<<"Getting numerator...";
  Integer translation, numExps[noCones], kk;
  int tmp_k = k;
  kk = tmp_k;
  // Get dot product and add negative denominator exponents for each cone.
  for(i=0;i<noCones;i++) {
    numExps[i]=0;
    for(j=0;j<dim;j++){ numExps[i]+=lambda[j]*A[i].front(); A[i].pop_front(); }
    for(j=0;j<noGsPerC;j++) {
      kk=dotProducts[noGsPerC*i+j];
      if(kk<0) {dotProducts[noGsPerC*i+j]=-kk; numExps[i]+=-kk; E[i]*= -1;}
    }
  }
  // cout<<"numexps: ";cin.get();
   // for(i=0;i<noCones;i++) cout<<numExps[i]<<" ";
  // cout<<endl;
   // return 0;
  // Now I Calculate into j the least translation possible for numExps...
  /*************************************************************************

   Rudy's comment:
   This step simplifies numerators (smaller numerators).
   In this process, he tries to get zeros in numerators as many as possible.
   For example, if we have the list of exponents for numerators for a list of
   cones,  4 4 3 5 3 6, then translation is 4 and after translating,
   these become 0 0 -1 2 -1 2 which are much smaller.  But, if we have a list
   of numerators -1 4 1, then a translation is 1 and the result becomes
   0 5 2.  If we have an element 1 or -1, then a translation becomes just 1.
   This step factors out t^c for some constant c.  It does not affest
   the result b/c at the end we substitute t = 1 (b/c t^c = 1 for all c).
   For example, if we have -1 4 1 means we have numerators t^(-1) + t^4 + t.
   So, at this step we get t^(-1)(1 + t^5 + t^2).  At the end, we substitute
   t = 1, so t^(-1) becomes 1.

  **************************************************************************/
  j=10000;
  int tmp_j = 0;

  for(i=0;i<dim;i++) {
    if(lambda[i]==1 || lambda[i]==-1) {
      i=dim+2; }
    else {
      if(lambda[i]<j && lambda[i]>0) {
       conv(tmp_j, lambda[i]);	j=tmp_j; }
      else {
	if(-lambda[i]<j && lambda[i]<0)  {
       conv(tmp_j,lambda[i]);	j=-tmp_j; }//j=-lambda[i];
      }
    }
  }
  if(i>dim) j=1;

  // Translate numExps.
  for(i=0;i<noCones;i+=100) translation+=numExps[i];
  translation/=1+noCones/100;
  translation-=translation%j;
  for(i=0;i<noCones;i++) numExps[i]-=translation;
 //    cout<<"translation: "<<translation<<"  numexps: ";
 //  for(i=0;i<noCones;i++) cout<<numExps[i]<<" ";
  //  cout<<endl;
  //   return 0;
  // cout<<".done"<<endl<<"  clock at checkpoint 3: "<<clock()<<endl;

  //-------------------------------------------------------------------------------
  //---CALCULATE AND SUM UP CONSTANT COEFFICIENTS IN SERIES EXPANSIONS ABOUT T=1---
  //---For each cone, which has form t^numExps[i]/prod(1-t^dotProducts[i*noGsPerC+j])
  //---we need to get its contribution to the number of lattice pts--
  //---this contribution is its constant coefficient in its expansion about t=1,
  //---so I make substitution t->s+1 and calculate coefficient noGsPerC of
  //---t^numExps[i]/prod((1-(1+s)^dotProducts[i*noGsPerC+j])/s).
  //---This involves expanding numerator and denominator and dividing.  I need
  //---only track up to coefficient noGsPerC in s in these expansions and division,
  //---so calculation time is of the order noGsPerC^3, bounded above by dim^3,
  //---for each cone.
  //---(Division is by play on the recursion q[k]=(n[k]-d[k]q[0]-...-d[1]q[k-1])/d[0]
  //---where q is quotient, n is numerator, d is denominator.)
  //-------------------------------------------------------------------------------

  // cout<<"Getting contributions in the big loop..." << endl;

  //VARS
  long int tenPow=10;  // To be used in getting extra digits of precision
  for(i=noCones;i>0;i/=10) tenPow*=10; // We will necessarily get enough precision.
  ZZ tempSum,temp,temp2,temp3,tempVec[1+noGsPerC];
  //for(i=0;i<=noGsPerC;i++) mpz_init(tempVec[i]);
  /*mpz_init(temp);
  mpz_init(temp2);
  mpz_init(temp3);
  mpz_init(tempSum);   */

  /***************************************************************************
   Rudy's comment:
   These arrays cause crash for latte.....  Example: 3x3x4_1.equ.residue.
   So, I changed to single arrays instead of double arrays.  We do not
   to have double arrays anyway b/c we conpute the coeffecient of the
   constant term at each time and add it to a variable.  So, it's never
   really used double array anyway.  It was wasting memory....
  ***************************************************************************/
  ZZ N[1+noGsPerC];
  ZZ D[1+noGsPerC];

  ZZ noLatticePts, nn;
  //mpz_init(noLatticePts);
  //mpz_set_str(noLatticePts,"0",10);

  /***************************************************************************
    Rudy's comment:
    Now, expanding each rational function and getting the coefficient of
    the constant term.....
  ****************************************************************************/

  //BIG LOOP
  for(i=0;i<noCones;i++) {
  for(int t = 0;t <= noGsPerC; t++) {N[t]=0; D[t]=0;}

  /***************************************************************************
    Rudy's comment:
    This is the first step.  It seems that this step substitute t -> s + 1.
    Also, in this step, his code expands denominator for each cone.
    So, I think D[i][j] stores the coefficients of a polynomial of
    each denominator.
    For example, if we have a rational function for the ith cone such that
    1/(1 - t)*( 1 - t^2), then after substitution, we have
    1/s^2(s + 2) = 1/s^3 + 2 * s^2.  SO, D[i][0] = 2, D[i][1] = 1,
    and D[i][2] = 0.  It seems that D[i][0] is the coefficient of s^noGsPerC
    and in this example, noGsPerC = 2.  Since polynomial in each denominator
    starts from noGsPerC, the degree of polynomial in the denominator is
    noGsPerC + 1 at most.  Thus, we allocate each denominator for each cone
    noGsPerC * mpz_t.

  ****************************************************************************/
    //MULTIPLY OUT DENOMINATOR
    sc=sc-clock();
    nn=dotProducts[i*noGsPerC];
    D[0]=nn; // get initial values for D
    for(j=1;j<=noGsPerC;j++) {
       temp=D[j-1]*(nn-j);
       D[j]=temp/(j+1);
    }
    for(j=noGsPerC*i+1;j<noGsPerC*(i+1);j++) { // multiply each factor into D
      nn=dotProducts[j];
      tempVec[0]=nn;   // First expand this denom. factor into tempVec
      for(k=1;k<=noGsPerC;k++) {  // k is what exp of t
   temp=tempVec[k-1]*(nn-k);
	tempVec[k]=temp/(k+1);
      }
      for(k=noGsPerC;k>=0;k--) {  // And then multimply by D; put product in D.
	tempSum=0;
	for(m=0;m<=k;m++) tempSum=tempSum+(D[m]*tempVec[k-m]);
	D[k]=tempSum;
      }
    }sc=sc+clock();

  /***************************************************************************
    Rudy's comment:
    The second step starts from here. The second step simplifies that
    if the exponent of numerator is negative, nultiply by (1+s)^(-numExps[i])
    and expand the denominator.
    If positive, expand the numerator.  For example, if we have a rational
    function for the ith cone t^2/s^2.  Then, N[i][0] = 1, N[i][1] = 2, and
    N[i][3] = 1.  However, If we have t^5/(2*s^2 + s^3),
    then, N[i][0] = 1, N[i][1] = 5, N[i][2] = 10.  I think we just don't need
    coefficients for the higher degree than noGsPerC b/c we just need the
    coefficient of the noGsPerC term.  SO, all we need is the coefficients
    of the first three terms.

  ****************************************************************************/

    //MULTIPLY (1+s)^abs(numExps[i]) INTO NUMERATOR OR DENOMINATOR
    sc2=sc2-clock();
    nn=numExps[i];
    N[0]=1;
    if(nn<0) { // Case: numExps[i]<0, so multiply denominator by (1+s)^(-numExps[i])
      nn=-nn;
      tempVec[0]=1; // put (1+s)^n into tempVec
      for(j=1;j<=noGsPerC;j++) {
	temp=tempVec[j-1]*(nn+1-j);
	tempVec[j]=temp/j;
      }
      for(k=noGsPerC;k>=0;k--) { // multiply tempVec into D
	tempSum=0;
	for(m=0;m<=k;m++) tempSum=tempSum+(D[m]*tempVec[k-m]);
	D[k]=tempSum;
      }
    }
    else { // Case numExps[i]>0; multiply out the numerator
      for(j=1;j<=noGsPerC;j++) {
	temp=N[j-1]*(nn+1-j);
	N[j]=temp/j;
      }
    }

  /***************************************************************************
    Rudy's comment:
    The third step starts from here. I think he stores all information into
    tempVec.  I think this is not a good idea.....  I do not understand the
    comment below too....  tempVec[j] stores the coefficient of the noGsPerc
    term for each cone.  So, at the end, we add up.....
  ****************************************************************************/
    //DIVISION: I track D[i][j]*D[i][0]^(j-1), N[i][j]*D[i][0]^j, and
    //--[jth coef. of quotient]*D[i][0]^(j+1), in row i of D, row i of N, and
    //--tempVec respectively, in order to preserve integer arithmetic.
    temp=1;   // track D[i][0]^(j-1) below
    tempVec[0]=N[0];
    for(j=1;j<=noGsPerC;j++) { // j is power of original D[i][0] we are on
      temp2=temp;
      temp=temp*D[0];
      tempVec[j]=N[j]*temp;
      D[j]=D[j]*temp2;
      for(m=1;m<=j;m++) tempVec[j]=tempVec[j]-(tempVec[j-m]*D[m]);
    }
    j--;
    temp=temp*D[0];
    tempVec[j]=tempVec[j]*tenPow*E[i];
    tempVec[j]=tempVec[j]/temp;

    //ADD CONTRIBUTION
    noLatticePts=noLatticePts+tempVec[j]; sc2=sc2+clock();
    //   cout<<" contrib: "<<mpz_get_str(NULL,10,tempVec[j])<<endl;
  } // end this long contribution loop

  //-----------------------------------------------------------------------------
  //--------------------FINISH UP AND DISPLAY RESULTS----------------------------
  //-----------------------------------------------------------------------------

  //REFINE (take abs and round) noLatticePts.
  noLatticePts=abs(noLatticePts); // case noGsPerC is odd (denom. factors)
  // cout<<".done"<<endl<<endl<<"tenPow: "<<tenPow<<endl;
//        <<"noLatticePts before division: "<<mpz_get_str(NULL,10,noLatticePts)<<endl;
  noLatticePts=noLatticePts+tenPow/2;
  noLatticePts=noLatticePts/tenPow;
  ofstream out("numOfLatticePoints");
  out << noLatticePts << endl;

  //OUTPUT TIMES AND RESULT
//    cout<<"denominator subclock: "<<sc<<endl<<"numerator subclock: "<<sc2<<endl
//        <<"clocks per second: "<<CLOCKS_PER_SEC<<endl<<endl;
//    cout<<"TIME FOR CALCULATION (in seconds): "<<(clock()-t)/CLOCKS_PER_SEC<<".";
  i=((clock()-t)%CLOCKS_PER_SEC)*100/CLOCKS_PER_SEC;
//    if(i<10) cout<<"0";
//    cout<<i<<endl;
//  cout<<endl<<"  **** THE NUMBER OF LATTICE POINTS IS: "<< noLatticePts <<" ****"<<endl<<endl;

  return noLatticePts;
}

/* ----------------------------------------------------------------- */


int
Residue_Single_Cone(listCone* cones, int numOfVars,
		    const vec_ZZ &Random_Lambda,
		    ZZ *Total_Lattice_Points, ZZ *Ten_Power) 
{

//  char outFileName[127];
  listCone *C, * cones1;
  int dim, noGsPerC,noCones; //noGsPerC is number of generators per cone
  clock_t t,sc=0,sc2=0;

 // strcpy(outFileName,fileName);
 // strcat(outFileName,".residue");

 /* ofstream out(outFileName);
  if (!out) {
    printf("Error opening output file for writing in printResidueFile!");
    exit(1);
  }
  if (cones==0) out << "No cones in list.\n";     */

  // Here we implicitly decompose each cone of index k into k "cones"
  // with just one lattice point.  --mkoeppe, Fri Mar 31 23:37:30 PST 2006
  noCones = 0;
  C=cones;
  while (C) {
    noCones += lengthListVector(cones->latticePoints);
    C=C->rest;
  }

  dim=numOfVars;
  noGsPerC=lengthListVector(cones->rays);
  int i,j;			// index or loop vars
  long int k, m;	       //n=0,p; // extra vars to use as needed
  vector<int> E(noCones);	  // E is the vector of epsilons, each 1 or -1
  vector<list<Integer> > A(noCones);	// A is the numerator vectors

  Integer tmp_A;
  int result = 1;

  long int totalNoGs=noGsPerC*noCones; //total no. of generators,ie,rowdim of B
  vector<Integer> dotProducts(totalNoGs);
  
  listVector* basis;
  listCone *listtmp3;
  cones1 = cones;
  i = 0;

  while (cones1) {
    listVector *tmp = cones1->latticePoints;
    while (tmp) {
      E[i] = cones1->coefficient;
      for (j=0; j<(numOfVars); j++) {  
	tmp_A = tmp->first[j];  
	A[i].push_back(tmp_A);
      }
      for (j=0, basis = cones1->rays; j<noGsPerC; j++, basis = basis->rest) {
	dotProducts[i*noGsPerC + j] = 0;
	for(k = 0; k < dim; k++) {
	  dotProducts[i*noGsPerC + j] += basis->first[k] * Random_Lambda[k];
	}
	// if the dot product is zero in the denominator, then barf
	if(dotProducts[i*noGsPerC + j] == 0)
	  result = -1;
      }
      assert(basis == NULL);
      tmp=tmp->rest; 
      i++;
    }
    listtmp3 = cones1;
    cones1 = cones1->rest;
    freeCone(listtmp3);
  }
  assert(i == noCones);

  if(result == -1)
    return result;

  i = 0;
 /* denom * Bitr=B;
  for(i=0;i<noCones;i++) {
    input>>E[i];
    for(j=0;j<dim;j++) {input>>tmp_A; A[i].push_back(tmp_A);}
    for(j=0;j<noGsPerC;j++)
      {
	for(k=0;k<dim;k++) input>>Bitr->D[j*dim+k];
      }
    Bitr->next=new denom(noGsPerC*dim);
    Bitr=Bitr->next;
  }
  input.close();  */

  
    /* Bitr=B;
     Bitr=Bitr->next;
     cout<<"B[1]: ";
     for(i=0;i<noGsPerC;i++) {cout<<endl;
     for(j=0;j<dim;j++) cout<<Bitr->D[i*dim+j]<<" ";
     }
     return ; */    


  // cout<<".done"<<endl<<"Number of cones: "<<noCones<<endl;cin.get();
//        <<"Number of generators in total: "<<totalNoGs<<endl
//        <<endl<<"Now start calculations and stopwatch."<<endl;
  t=clock();
  //  cout<<"  Clock at start reads "<<t<<"."<<endl;

  //----------------------------------------------------------------------------
  //---CALCULATE NUMERATOR EXPONENTS of t under substitution x[i]->t^lambda[i].
  //---After getting initial exponents, I translate so as to try to minimize
  //---numerator exponents (and promote efficientcy).
  //----------------------------------------------------------------------------


  /****************************************************************************
    Rudy's comment:
    This is a process to simplify a univariable rational function.
    For example, if I have a rational function f(K(v)) = t/(1-t^-1)(1-t^-2),
    then, the following loop simplifies as f(K(v)) = t^4 /(1-t)(1-t^2).

  *****************************************************************************/
  //  cout<<"Getting numerator...";
  Integer translation, numExps[noCones], kk;
  int tmp_k = k;
  kk = tmp_k;
  // Get dot product and add negative denominator exponents for each cone.
  for(i=0;i<noCones;i++) {
    numExps[i]=0;
    for(j=0;j<dim;j++){ numExps[i]+=Random_Lambda[j]*A[i].front(); A[i].pop_front(); }
    for(j=0;j<noGsPerC;j++) {
      kk=dotProducts[noGsPerC*i+j];
      if(kk<0) {dotProducts[noGsPerC*i+j]=-kk; numExps[i]+=-kk; E[i]*= -1;}
    }
  }
  // cout<<"numexps: ";cin.get();
   // for(i=0;i<noCones;i++) cout<<numExps[i]<<" ";
  // cout<<endl;
   // return 0;
  // Now I Calculate into j the least translation possible for numExps...
  /*************************************************************************

   Rudy's comment:
   This step simplifies numerators (smaller numerators).
   In this process, he tries to get zeros in numerators as many as possible.
   For example, if we have the list of exponents for numerators for a list of
   cones,  4 4 3 5 3 6, then translation is 4 and after translating,
   these become 0 0 -1 2 -1 2 which are much smaller.  But, if we have a list
   of numerators -1 4 1, then a translation is 1 and the result becomes
   0 5 2.  If we have an element 1 or -1, then a translation becomes just 1.
   This step factors out t^c for some constant c.  It does not affest
   the result b/c at the end we substitute t = 1 (b/c t^c = 1 for all c).
   For example, if we have -1 4 1 means we have numerators t^(-1) + t^4 + t.
   So, at this step we get t^(-1)(1 + t^5 + t^2).  At the end, we substitute
   t = 1, so t^(-1) becomes 1.

  **************************************************************************/
  /*
  j=10000;
  int tmp_j = 0;

  for(i=0;i<dim;i++) {
    if(Random_Lambda[i]==1 || Random_Lambda[i]==-1) {
      i=dim+2; }
    else {
      if(Random_Lambda[i]<j && Random_Lambda[i]>0) {
       conv(tmp_j, Random_Lambda[i]);	j=tmp_j; }
      else {
	if(-Random_Lambda[i]<j && Random_Lambda[i]<0)  {
       conv(tmp_j,Random_Lambda[i]);	j=-tmp_j; }//j=-[i];
      }
    }
  }
  if(i>dim) j=1;
*/
  // Translate numExps.
  /*for(i=0;i<noCones;i+=100) translation+=numExps[i];
  translation/=1+noCones/100;
  translation-=translation%j;
  for(i=0;i<noCones;i++) numExps[i]-=translation;
  */
  //    cout<<"translation: "<<translation<<"  numexps: ";
 //  for(i=0;i<noCones;i++) cout<<numExps[i]<<" ";
  //  cout<<endl;
  //   return 0;
  // cout<<".done"<<endl<<"  clock at checkpoint 3: "<<clock()<<endl;

  //-------------------------------------------------------------------------------
  //---CALCULATE AND SUM UP CONSTANT COEFFICIENTS IN SERIES EXPANSIONS ABOUT T=1---
  //---For each cone, which has form t^numExps[i]/prod(1-t^dotProducts[i*noGsPerC+j])
  //---we need to get its contribution to the number of lattice pts--
  //---this contribution is its constant coefficient in its expansion about t=1,
  //---so I make substitution t->s+1 and calculate coefficient noGsPerC of
  //---t^numExps[i]/prod((1-(1+s)^dotProducts[i*noGsPerC+j])/s).
  //---This involves expanding numerator and denominator and dividing.  I need
  //---only track up to coefficient noGsPerC in s in these expansions and division,
  //---so calculation time is of the order noGsPerC^3, bounded above by dim^3,
  //---for each cone.
  //---(Division is by play on the recursion q[k]=(n[k]-d[k]q[0]-...-d[1]q[k-1])/d[0]
  //---where q is quotient, n is numerator, d is denominator.)
  //-------------------------------------------------------------------------------

  // cout<<"Getting contributions in the big loop..." << endl;

  //VARS
  //long int tenPow=10;  // To be used in getting extra digits of precision
  //for(i=noCones;i>0;i/=10) tenPow*=10; // We will necessarily get enough precision.
  ZZ tempSum,temp,temp2,temp3,tempVec[1+noGsPerC];
  //for(i=0;i<=noGsPerC;i++) mpz_init(tempVec[i]);
  /*mpz_init(temp);
  mpz_init(temp2);
  mpz_init(temp3);
  mpz_init(tempSum);   */

  /***************************************************************************
   Rudy's comment:
   These arrays cause crash for latte.....  Example: 3x3x4_1.equ.residue.
   So, I changed to single arrays instead of double arrays.  We do not
   to have double arrays anyway b/c we conpute the coeffecient of the
   constant term at each time and add it to a variable.  So, it's never
   really used double array anyway.  It was wasting memory....
  ***************************************************************************/
  ZZ N[1+noGsPerC];
  ZZ D[1+noGsPerC];

  ZZ noLatticePts, nn;
  
  noLatticePts = 0;
  //mpz_init(noLatticePts);
  //mpz_set_str(noLatticePts,"0",10);

  /***************************************************************************
    Rudy's comment:
    Now, expanding each rational function and getting the coefficient of
    the constant term.....
  ****************************************************************************/

  //BIG LOOP
  for(i=0;i<noCones;i++) {
  for(int t = 0;t <= noGsPerC; t++) {N[t]=0; D[t]=0;}

  /***************************************************************************
    Rudy's comment:
    This is the first step.  It seems that this step substitute t -> s + 1.
    Also, in this step, his code expands denominator for each cone.
    So, I think D[i][j] stores the coefficients of a polynomial of
    each denominator.
    For example, if we have a rational function for the ith cone such that
    1/(1 - t)*( 1 - t^2), then after substitution, we have
    1/s^2(s + 2) = 1/s^3 + 2 * s^2.  SO, D[i][0] = 2, D[i][1] = 1,
    and D[i][2] = 0.  It seems that D[i][0] is the coefficient of s^noGsPerC
    and in this example, noGsPerC = 2.  Since polynomial in each denominator
    starts from noGsPerC, the degree of polynomial in the denominator is
    noGsPerC + 1 at most.  Thus, we allocate each denominator for each cone
    noGsPerC * mpz_t.

  ****************************************************************************/
    //MULTIPLY OUT DENOMINATOR
    sc=sc-clock();
    nn=dotProducts[i*noGsPerC];
    D[0]=nn; // get initial values for D
    for(j=1;j<=noGsPerC;j++) {
       temp=D[j-1]*(nn-j);
       D[j]=temp/(j+1);
    }
    for(j=noGsPerC*i+1;j<noGsPerC*(i+1);j++) { // multiply each factor into D
      nn=dotProducts[j];
      tempVec[0]=nn;   // First expand this denom. factor into tempVec
      for(k=1;k<=noGsPerC;k++) {  // k is what exp of t
   temp=tempVec[k-1]*(nn-k);
	tempVec[k]=temp/(k+1);
      }
      for(k=noGsPerC;k>=0;k--) {  // And then multimply by D; put product in D.
	tempSum=0;
	for(m=0;m<=k;m++) tempSum=tempSum+(D[m]*tempVec[k-m]);
	D[k]=tempSum;
      }
    }sc=sc+clock();

  /***************************************************************************
    Rudy's comment:
    The second step starts from here. The second step simplifies that
    if the exponent of numerator is negative, nultiply by (1+s)^(-numExps[i])
    and expand the denominator.
    If positive, expand the numerator.  For example, if we have a rational
    function for the ith cone t^2/s^2.  Then, N[i][0] = 1, N[i][1] = 2, and
    N[i][3] = 1.  However, If we have t^5/(2*s^2 + s^3),
    then, N[i][0] = 1, N[i][1] = 5, N[i][2] = 10.  I think we just don't need
    coefficients for the higher degree than noGsPerC b/c we just need the
    coefficient of the noGsPerC term.  SO, all we need is the coefficients
    of the first three terms.

  ****************************************************************************/

    //MULTIPLY (1+s)^abs(numExps[i]) INTO NUMERATOR OR DENOMINATOR
    sc2=sc2-clock();
    nn=numExps[i];
    N[0]=1;
    if(nn<0) { // Case: numExps[i]<0, so multiply denominator by (1+s)^(-numExps[i])
      nn=-nn;
      tempVec[0]=1; // put (1+s)^n into tempVec
      for(j=1;j<=noGsPerC;j++) {
	temp=tempVec[j-1]*(nn+1-j);
	tempVec[j]=temp/j;
      }
      for(k=noGsPerC;k>=0;k--) { // multiply tempVec into D
	tempSum=0;
	for(m=0;m<=k;m++) tempSum=tempSum+(D[m]*tempVec[k-m]);
	D[k]=tempSum;
      }
    }
    else { // Case numExps[i]>0; multiply out the numerator
      for(j=1;j<=noGsPerC;j++) {
	temp=N[j-1]*(nn+1-j);
	N[j]=temp/j;
      }
    }

  /***************************************************************************
    Rudy's comment:
    The third step starts from here. I think he stores all information into
    tempVec.  I think this is not a good idea.....  I do not understand the
    comment below too....  tempVec[j] stores the coefficient of the noGsPerc
    term for each cone.  So, at the end, we add up.....
  ****************************************************************************/
    //DIVISION: I track D[i][j]*D[i][0]^(j-1), N[i][j]*D[i][0]^j, and
    //--[jth coef. of quotient]*D[i][0]^(j+1), in row i of D, row i of N, and
    //--tempVec respectively, in order to preserve integer arithmetic.
    temp=1;   // track D[i][0]^(j-1) below
    tempVec[0]=N[0];
    for(j=1;j<=noGsPerC;j++) { // j is power of original D[i][0] we are on
      temp2=temp;
      temp=temp*D[0];
      tempVec[j]=N[j]*temp;
      D[j]=D[j]*temp2;
      for(m=1;m<=j;m++) tempVec[j]=tempVec[j]-(tempVec[j-m]*D[m]);
    }
    j--;
    temp=temp*D[0];
    tempVec[j]=tempVec[j]*(*Ten_Power)*E[i];
    tempVec[j]=tempVec[j]/temp;

    //ADD CONTRIBUTION
    noLatticePts=noLatticePts+tempVec[j]; sc2=sc2+clock();
    //   cout<<" contrib: "<<mpz_get_str(NULL,10,tempVec[j])<<endl;
  } // end this long contribution loop

  //-----------------------------------------------------------------------------
  //--------------------FINISH UP AND DISPLAY RESULTS----------------------------
  //-----------------------------------------------------------------------------

  //REFINE (take abs and round) noLatticePts.
  //noLatticePts=abs(noLatticePts); // case noGsPerC is odd (denom. factors)
  // cout<<".done"<<endl<<endl<<"tenPow: "<<tenPow<<endl;
//        <<"noLatticePts before division: "<<mpz_get_str(NULL,10,noLatticePts)<<endl;
  //noLatticePts=noLatticePts+tenPow/2;
  //noLatticePts=noLatticePts/tenPow;
	
 // 	cout << "Residue_Single_Cone: " << 
  
  *Total_Lattice_Points += noLatticePts;
  
  //OUTPUT TIMES AND RESULT
//    cout<<"denominator subclock: "<<sc<<endl<<"numerator subclock: "<<sc2<<endl
//        <<"clocks per second: "<<CLOCKS_PER_SEC<<endl<<endl;
//    cout<<"TIME FOR CALCULATION (in seconds): "<<(clock()-t)/CLOCKS_PER_SEC<<".";
  i=((clock()-t)%CLOCKS_PER_SEC)*100/CLOCKS_PER_SEC;
//    if(i<10) cout<<"0";
//    cout<<i<<endl;

  return 1;
  
}

/* ----------------------------------------------------------------- */
