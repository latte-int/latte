/* RudyResNTL.cpp -- Polynomial substitution and residue calculations

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

#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include <cstdlib>
#include <list>
#include <vector>
#include <cassert>
#include "PolyTree.h"
#include "cone.h"
#include "ramon.h"
#include "RudyResNTL.h"
#include <malloc.h>
#include "latte_system.h"

using namespace std;

/***********************************************************************************/
void ResidueFunction(listCone* cones, int numOfVars, int print_flag, 
    int degree, int output_cone) {
  int numOfTerms, DEGREE = 1;
//  char outFileName[127];
  listVector *tmp;
  static int never_printed = 1;
  listCone *C, * cones1;
  int OUTPUT_S_T_FUNCTION = 0, COMPUTE_SUM_MAPLE = 0, OUTPUT_CONE = 0, OUTPUT_CONE_MULTI = 0;
  int dim, noGsPerC,noCones; //noGsPerC is number of generators per cone
  clock_t t;//,sc=0,sc2=0;
  if(degree != 0) DEGREE = degree;
 // strcpy(outFileName,fileName);
 // strcat(outFileName,".residue");
 if(output_cone == 1) OUTPUT_CONE = 1;
 else if(output_cone == 2) COMPUTE_SUM_MAPLE = 1;
 else if(output_cone == 3) OUTPUT_S_T_FUNCTION = 1;

 	int	Max_Cones_File = 1000;
	int	Cones_File_Count = 0;

 /* ofstream out(outFileName);
  if (!out) {
    printf("Error opening output file for writing in printResidueFile!");
    exit(1);
  }
  if (cones==0) out << "No cones in list.\n";     */

  numOfTerms=0;

  C=cones;
  while (C) {
    assert(IsZero(C->vertex->numerators()));
    assert(abs(C->determinant) == 1);
    numOfTerms=numOfTerms+lengthListVector(C->latticePoints);
    C=C->rest;
  }


 /* out << numOfVars << " " << lengthListVector(cones->rays) << " " <<
    numOfTerms << "\n\n";  */

  dim=numOfVars;
  noGsPerC=lengthListVector(cones->rays);
  noCones=numOfTerms;

  vector<Cone_Data> Cones_Array(noCones);  // Create structure to hold all the information of the cones.
  
  for (int q = 0;q < noCones; q++)  // Do some initialization
  {
	  Cones_Array[q].order = 0;
	  Cones_Array[q].Generators_of_Cone = new Generator [noGsPerC];  //Each cone has noGsPerC generators.
  }
  // Added by P/D
  int i,j; // index or loop vars
  long int k, m;//n=0,p; // extra vars to use as needed
  long int p, n;
  vector<int> E(noCones);  // E is the vector of epsilons, each 1 or -1
  long int totalNoGs=noGsPerC*noCones; //total no. of generators,ie,rowdim of B
  vector<list<Integer> > A(noCones);  // A is the numerator vectors
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
    while (tmp) {
      E[i] = cones1->coefficient;
      Cones_Array[i].sign = E[i];   //Added by P/D:  Load what sign each cone is
     // printVectorToFileWithoutBrackets(out,tmp->first,numOfVars);
    for (j=0; j<(numOfVars); j++) {
    tmp_A = tmp->first[j];  A[i].push_back(tmp_A);
      }
      Cones_Array[i].Numerator_Generator.T_Exponent=tmp_A;
     // printListVectorToFileWithoutBrackets(out,cones->rays,numOfVars);
     basis = cones1->rays;
     while(basis) {
    //printVectorToFileWithoutBrackets(out,basis->first,numOfVars);
    for (j=0; j<noGsPerC; j++) {
    for(k = 0; k <dim; k++){
    Bitr->D[j*dim+k] = basis->first[k];
      }
      // Check when k=dim-1 if input is 0, increment order of cone
	Cones_Array[i].Generators_of_Cone[j].T_Exponent=Bitr->D[j*dim + dim -1];  //exponent of t in the denominator	
	Cones_Array[i].Generators_of_Cone[j].Form_Type = ONE_SUB_RT;
	if (Bitr->D[j*dim + dim - 1] == 0)  //if the exponent of t is zero, increment the order for this cone
		Cones_Array[i].order++;
      
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
  vector<Integer> lambda(dim);
  //lambda[0] = 327;
  vector<Integer> dlambda(dim); // dlambda tracks change in 2 successive test-lambdas
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
	for(p=n;p<dim-1;p++) dotProducts[m]+= dlambda[p] * Bitr->D[j*dim+p];
	if(dotProducts[m]==0) {haltCone=i; halt=j; j=noGsPerC; i=noCones+2;}
      }
      if(i<noCones) Bitr=Bitr->next;
    }
    for(;i<noCones;i++) { // keep checking if not yet halted
      for(j=0;j<noGsPerC;j++) {
	m=i*noGsPerC+j;
	dotProducts[m]=0;
	for(p=0;p<dim-1;p++) dotProducts[m] += lambda[p] * Bitr->D[j*dim+p];
	if(dotProducts[m]==0) {haltCone=i; halt=j; j=noGsPerC; i=noCones+2;}
      }
      if(i<noCones) Bitr=Bitr->next;
    }
    if(i==noCones) { // found lambda
      k=200000; }
    else {

      //--SECOND find next lambda.
      m= (Bitr->D[halt*dim+dim-2] ==0 ? 0 : 1);
	n=dim-2;
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
 // cout << totalNoGs << end;
 // 	
 //for (int q = 0;q<dim-1;q++)
//	 cout << lambda[q] << " ";
 //cout << " Lambda " << endl;

  /*for(j = 0; j < noCones; j++) 
	  if(dotProducts[j] == 0) 
		  cout <<  "oops";

  exit(1);*/
  
  // cout << endl;*/

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
 //int d=1; // d runs through odds... it gets added to a single entry of lambda
       // each time through the loop, to fix the problem generator.
  //if(k<15000) {for(i=0;i<dim;i++) lambda[i]=-1;} if running random algorithm
  //halt=0; if running random algorithm
  int ss = 0;
  while(k < 150000) {
//      ss = ( k - 9960) / 10;
//      for(int s = 0;s < dim; s++) lambda[s] = rand() % ss - ss/2;
    //FIRST find new lambda.
    ss++;
    n = 0;
    //  ss = ( k - 9960) / 10;
    for(int s = 0;s < dim; s++)
      {lambda[s] = rand() % 1000;
      if(rand() % 2 == 1)  lambda[s] = -lambda[s];
      else ;
      }
    /*    while(Bitr->D[halt*dim+n]==0) n++;
    //d+=2;
    d = rand() % 1000;
    lambda[n]+=d; */

    //SECOND check if lambda works
    Bitr=B;
    /*    for(i=0;i<haltCone;i++) {
      for(j=0;j<noGsPerC;j++) {
	m=i*noGsPerC+j;
	dotProducts[m]+=d*Bitr->D[j*dim+n];
	if (dotProducts[m]==0) {haltCone=i; halt=j; j=noGsPerC; i=noCones+2;}
      }
      if(i<noCones) Bitr=Bitr->next;
      }*/
    for(i = 0;i<noCones;i++) {
      for(j=0;j<noGsPerC;j++) {
	m=i*noGsPerC+j;
	dotProducts[m]=0;
	for(p=0;p<dim;p++) dotProducts[m]+=lambda[p]*Bitr->D[j*dim+p];
	if(dotProducts[m]==0) {haltCone=i; halt=j; j=noGsPerC; i=noCones+2;}
      }
      if(i<noCones) Bitr=Bitr->next;
    }
    if(i==noCones) k=200000;
  }

  // cout<<"lambda: "; cin.get();
//  for(i=0;i<dim;i++) cout<<lambda[i]<<" ";
 // cout << endl;
 // cout << ss << endl;
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

  
 

  
 	//cout << "Copying dot product into Cones_Array" << endl; 
  //*****************************************************************
  //  PETER/DAVE CODE REALLY BEGINS HERE
  //  Our data structure Cones_Array which holds all of our information
  //  is currently loaded with the sign of each cone, the
  //  exponent of t in the numerators, and the exponent of t
  //  for all the generators.

  //  COPY DOTPRODUCT INTO CONES_ARRAY STRUCTURE.
  //  CALCLUTE DOT PRODUCT OF THE NUMERATOR AND STORE IN CONES_ARRAY
  	for (int q = 0;q < noCones;q++)  // For each cone
	{
		Cones_Array[q].Numerator_Generator.R_Exponent = 0;
	
		//cout << "Calculationg dot product of numerator" << endl;	
		// Calculate dot product of numerator
		for (int f = 0;f < dim-1 ;f++)
		{
			Cones_Array[q].Numerator_Generator.R_Exponent += lambda[f]*A[q].front();
			A[q].pop_front ();
		}
		
		A[q].pop_front ();  //remove last dimension's exponent
		
		//cout << "Store dot product of generators into Cones_Array" << endl;
		// Store dot product of generators into Cones_Array
		for (int t = 0;t < noGsPerC;t++)
		{
			Cones_Array[q].Generators_of_Cone[t].R_Exponent = dotProducts[q*noGsPerC + t];
		}
	}
	
	
	//cout << "Simplifying generators to have nonnegative exponents" << endl;
  //**************************************************************************
  //  Simplify all the generators of each cone such that all the
  //  exponents are nonnegative, changing the sign and Form_Type accordingly.
  //  Also, we look through all the Exponents on R of each cone and 
  //  record the minimum into Numerator_R_Exponent_Minumum.
  //  We use the minumum afterwards to make the exponents of the numerators
  //  nonnegative for every cone
	
	ZZ	Numerator_R_Exponent_Minimum; //Used to store the minumum exponent
	
	for (int q = 0;q < noCones; q++)
	{
		for (int t = 0;t < noGsPerC; t++)
		{
			if (Cones_Array[q].Generators_of_Cone[t].R_Exponent < 0)
			{       	
				// R_Exponent < 0 and T_Exponent < 0
				if (Cones_Array[q].Generators_of_Cone[t].T_Exponent <= 0)
				{
					Cones_Array[q].sign *= -1;
				
					Cones_Array[q].Generators_of_Cone[t].R_Exponent *= -1;
					Cones_Array[q].Numerator_Generator.R_Exponent += 
						Cones_Array[q].Generators_of_Cone[t].R_Exponent;

					Cones_Array[q].Generators_of_Cone[t].T_Exponent *= -1;
					Cones_Array[q].Numerator_Generator.T_Exponent +=
						Cones_Array[q].Generators_of_Cone[t].T_Exponent;
				
				}
				else // R_Exponent < 0 and T_Exponent > 0 
				{
					Cones_Array[q].Generators_of_Cone[t].Form_Type = R_SUB_T; //(r-t)
					
					Cones_Array[q].Generators_of_Cone[t].R_Exponent *= -1;
					Cones_Array[q].Numerator_Generator.R_Exponent += 
						Cones_Array[q].Generators_of_Cone[t].R_Exponent;
				}
			}
			else // R_Exponent > 0   Check T_Exponent
			{
				// R_Exponent > 0 and T_Exponent < 0
				if (Cones_Array[q].Generators_of_Cone[t].T_Exponent < 0)
				{
					Cones_Array[q].sign *= -1;

					Cones_Array[q].Generators_of_Cone[t].Form_Type = R_SUB_T; // (r-t)

					Cones_Array[q].Generators_of_Cone[t].T_Exponent *= -1;
					Cones_Array[q].Numerator_Generator.T_Exponent +=
						Cones_Array[q].Generators_of_Cone[t].T_Exponent;


				}
				// Otherwise R_Exponent and T_Exponent > 0 so we change nothing

			}
		//  Initialize the minimum to the first R_Exponent of the cones
		//  Runs only the first pass of the for loop
		if (q == 0)
			Numerator_R_Exponent_Minimum = Cones_Array[q].Numerator_Generator.R_Exponent;

		//  If we find an element smaller than our current minimum
		//  record it as our new minumum
		if (Cones_Array[q].Numerator_Generator.R_Exponent < Numerator_R_Exponent_Minimum)
				Numerator_R_Exponent_Minimum = Cones_Array[q].Numerator_Generator.R_Exponent;
			
			
		}

	}

	//cout << "Minimum exponent is " << Numerator_R_Exponent_Minimum << endl;
	//cout << "Factoring out minumum exponent of r" << endl;
  
  //************************************************************************
  //  Peter/Dave: Here we make sure all the numerators of the Cones are
  //  nonnegative.  We use the Minumum of the numerators calculated previously.
  //  We do this assuming that we can factor out any number of R's, as long
  //  as we factor out the same amount from each cone.  R = 1 so this should
  //  be ok.   	
	for (int q = 0; q < noCones; q++)
	{
		Cones_Array[q].Numerator_Generator.R_Exponent -= Numerator_R_Exponent_Minimum;
	}


	// Create all the variables we are going to use in our big loop
	
	vector<PolyTree_Node *> Numerator_Vector(noGsPerC + 1); // [i] - coefficient of s^i
	vector<PolyTree_Node *> Denominator_Result(noGsPerC + 1);   // Used to store the running total of generators
	vector<PolyTree_Node *> Denominator_Current_Generator(noGsPerC + 1);  //For each generatore as we iterate through them
	vector<PolyTree_Node *>	Quotient_Coefficient(noGsPerC + 1);  // holds the coefficients for calculting residue
							      // This structure does explicitly hold the value of the
							      // b_0 denominator.  It is implied by its index b_0^i+1	
	
		
	PolyTree_Node	*Coefficient_Addition_Root;		//Used to create new addition roots
	PolyTree_Node	*Coefficient_Multiplication_Root;	//Used to create new mul roots
	PolyTree_Node	*Quotient_Root;				//Used to create new quotient roots
	PolyTree_Node	*Coefficient_Exponent_Root;			//Used to create new exponent roots
	PolyTree_Node	*Final_Cone_Expression;			//We are going to store our final expression
								//here as one addition root
		
	T_Node  *New_T_Node;					//Used to create new T_Node roots
	T_Node  *T_Node_Negative_One;				//A Node that is equal to -1

	ZZ	Numerical_Coefficient; //Used for calculating coefficient of binomial expansion of (s+1)^p
	ZZ	Temp;  // Used for calculation coefficient
	ZZ	Exponent_of_T;  //Form_Type dictates what this will be

	int	Exponent_Reduction_Offset;  //Used when we factor out 1/s^(order) of each cone
					// It is equal to 0 or 1 for each generator depending
					// if it is a pole or not 
	
	// Initialize our final expression


	cout << endl << endl;
	
	Node_Controller	Controller (noGsPerC, DEGREE);

	vector<ZZ>	Final_Taylor_Result(DEGREE + 1);
	
	for (int i = 0; i <= DEGREE; i++)
	{
		Final_Taylor_Result[i] = 0;
	}
	
	Taylor_Parameters *Cone_Taylor_Parameters = new Taylor_Parameters;
	ZZ	*Ten_Power = new ZZ;
	
	*Ten_Power = 10;
	
	// Ten_Power = 10^( ceil[log(noCones)] + 1)
	for (int i = noCones;i > 0; i /= 10)
		*Ten_Power *= 10;
	
	//Calculate the taylor expasion for the first DEGREE + 1 many terms
	Cone_Taylor_Parameters->Result = new ZZ [DEGREE + 1];
	Cone_Taylor_Parameters->Ten_Power = Ten_Power;
	Cone_Taylor_Parameters->Degree_of_Expansion = DEGREE;
	ofstream Rational_Function_Output_File;
		
		
	if(print_flag == 1)
		 {
		   //system_with_error_check("rm func.rat");
	 		 cout << "Outputing rational functions to file" << endl;
			Rational_Function_Output_File.open ("func.rat");
		 }

	cout << "Formulating rational functions and performing taylor expansion on cones." << endl;

	ofstream Simplify_Sum, Simplify_Term;
	
	if (COMPUTE_SUM_MAPLE == 1)
	{
		cout << "Compute maple called" << endl;
		//system_with_error_check("rm simplify.sum");
		
		//Create initial sum file simplify.sum
		Simplify_Sum.open ("simplify.sum");

		Simplify_Sum << "s := 0:";
		
		Simplify_Sum.close ();
	}

	ofstream Rational_Function_S_T;

	if( OUTPUT_S_T_FUNCTION == 1)
	{
		Simplify_Sum.open ("simplify.sum");
		
		Simplify_Sum << "HS := 0:";

		Simplify_Sum.close ();	
		

	}
	
	if ( OUTPUT_CONE == 1)
	{
		Simplify_Sum.open ("simplify.sum");

		Simplify_Sum << "HS := 0:";

		Simplify_Sum.close ();
	}

	if ( OUTPUT_CONE_MULTI == 1)
	{
		Simplify_Sum.open ("simplify.sum");

		Simplify_Sum << "HS := 0:";

		Simplify_Sum.close ();
	}	

	// For each Cone
	for (int i = 0; i < noCones; i++)
	{
		//Calculate the numerator vector
		New_T_Node = Controller.Get_T_Node ();
		
		New_T_Node->Node_Type = POLYTREE_T_NODE;
		New_T_Node->Coefficient = 1;
		New_T_Node->Exponent = Cones_Array[i].Numerator_Generator.T_Exponent;
		Numerator_Vector[0] = New_T_Node;
		
		Temp = Cones_Array[i].Numerator_Generator.R_Exponent;
		Numerical_Coefficient = 1;
		
		
		// The "Choose" function for the coefficients of the numerator
		for (int k = 1; k <= Cones_Array[i].order; k++)
		{
			Numerical_Coefficient *= Temp;
			Numerical_Coefficient /= k;

			New_T_Node = Controller.Get_T_Node ();
			New_T_Node->Node_Type = POLYTREE_T_NODE;
			New_T_Node->Coefficient = Numerical_Coefficient;
			New_T_Node->Exponent = Cones_Array[i].Numerator_Generator.T_Exponent;
			Numerator_Vector[k] = New_T_Node;
			Temp--;	
			
		}

		if ( OUTPUT_CONE == 1 )
		{
			Simplify_Term.open ("simplify.term");
			
			Simplify_Term << " d := " << Cones_Array[i].order << ":" << endl;

			if (Cones_Array[i].sign == 1)
				Simplify_Term << "x := (";
			else
				Simplify_Term << "x := (-1)*(";

			for (int g = 0; g <= Cones_Array[i].order; g++)
			{
				Numerator_Vector[g]->Print_Rational_Functions_to_File (Simplify_Term);
				if ( g != 0)
					Simplify_Term << "*s^" << g;
				if ( g != Cones_Array[i].order)
					Simplify_Term << "+";
			}
				
			Simplify_Term << ")/(";


		}	
		if ( OUTPUT_CONE_MULTI == 1 )
		{
			// Open to appendto simplify.term
			if (Cones_File_Count == 0)
				system_with_error_check ("rm -f simplify.term");

			Simplify_Term.open ("simplify.term");
			
			if (Cones_File_Count == 0)
			{
					
				Simplify_Term << " d := " << Cones_Array[i].order << ":" << endl;
			}
				
			if (Cones_Array[i].sign == 1)
			{
				Simplify_Term << "x" << Cones_File_Count << " := (";
			}
			else
			{
				Simplify_Term << "x" << Cones_File_Count << " := (-1)*(";
			}
			for (int g = 0; g <= Cones_Array[i].order; g++)
			{
				Numerator_Vector[g]->Print_Rational_Functions_to_File (Simplify_Term);
				if ( g != 0)
					Simplify_Term << "*s^" << g;
				if ( g != Cones_Array[i].order)
					Simplify_Term << "+";
			}
				
			Simplify_Term << ")/(";
		
			//Simplify_Term.close ();	
			Cones_File_Count++;

		}
		// Initialize the Denominator result with the information of the first generator, 
		// (only if the T_Exponent in non-zero, i.e. this generator is not a pole)

		
		
		if(Cones_Array[i].Generators_of_Cone[0].T_Exponent != 0)
		{
			Denominator_Result[0] = Controller.Get_PolyTree_Node ();
			Denominator_Result[0]->Node_Type = POLYTREE_ADD; // +
			Denominator_Result[0]->Number_of_Children = 2;
		
			New_T_Node = Controller.Get_T_Node ();
			New_T_Node->Node_Type = POLYTREE_T_NODE;
			New_T_Node->Coefficient = 1;
			New_T_Node->Exponent = 0;
			
			Denominator_Result[0]->Children[0] = New_T_Node;

			New_T_Node = Controller.Get_T_Node ();
			New_T_Node->Node_Type = POLYTREE_T_NODE;
			New_T_Node->Coefficient = -1;
			New_T_Node->Exponent = Cones_Array[i].Generators_of_Cone[0].T_Exponent;
			Denominator_Result[0]->Children[1] = New_T_Node;
			
			// This generator is not pole, thus do not factor out an 's' from it
			Exponent_Reduction_Offset = 0;	
		}
		else
		{
			// This simulates that we are reducing the power of 's' by one since this generator is a pole
			Exponent_Reduction_Offset = 1;
		
		}

		// Expand this cone according to its type, (1-rt) or (r-t)
		if (Cones_Array[i].Generators_of_Cone[0].Form_Type == ONE_SUB_RT) // means (1-rt)
		{
			Exponent_of_T = Cones_Array[i].Generators_of_Cone[0].T_Exponent;
			Numerical_Coefficient = -1;		
		}
		else // Form_Type is (r-t)
		{
			Exponent_of_T = 0;
			Numerical_Coefficient = 1;
		}
		
		Temp = Cones_Array[i].Generators_of_Cone[0].R_Exponent;

		// The "choose" function to calculate the coefficients
		for (int k = 1; k <= Cones_Array[i].order + Exponent_Reduction_Offset; k++)
		{
			Numerical_Coefficient *= Temp;
			Numerical_Coefficient /= k;

			New_T_Node = Controller.Get_T_Node ();
			New_T_Node->Node_Type = POLYTREE_T_NODE;
			New_T_Node->Coefficient = Numerical_Coefficient;
			New_T_Node->Exponent = Exponent_of_T;
			Denominator_Result[k - Exponent_Reduction_Offset] = New_T_Node;
			
			Temp--;
		}

	       
		if ( OUTPUT_CONE == 1 )
		{
		  Simplify_Term << "(";

			for (int g = 0; g <= Cones_Array[i].order; g++)
			{
				Denominator_Result[g]->Print_Rational_Functions_to_File (Simplify_Term);
				if ( g != 0)
					Simplify_Term << "*s^" << g;
				if ( g != Cones_Array[i].order)
					Simplify_Term << "+";
			}
				
			Simplify_Term << ")*";


		}	

		
		
		// Now Denominator_Result holds the information for the first generator.
		// Now we can fold the rest of the generators into Denominator_Result.

		//For each generator starting at 1;		
		for (int q = 1; q < noGsPerC; q++)
		{
			
		        // Take the qth Generator and convert to our Deninator_Current_Generator
			
			
			if (Cones_Array[i].Generators_of_Cone[q].T_Exponent != 0) //Not a pole
			{	
				Denominator_Current_Generator[0] = Controller.Get_PolyTree_Node ();
				Denominator_Current_Generator[0]->Node_Type = POLYTREE_ADD; // +
				Denominator_Current_Generator[0]->Number_of_Children = 2;
		
				New_T_Node = Controller.Get_T_Node ();
				New_T_Node->Node_Type = POLYTREE_T_NODE;
				New_T_Node->Coefficient = 1;
				New_T_Node->Exponent = 0;
				Denominator_Current_Generator[0]->Children[0] = New_T_Node;

				New_T_Node = Controller.Get_T_Node ();
				New_T_Node->Node_Type = POLYTREE_T_NODE;
				New_T_Node->Coefficient = -1;
				New_T_Node->Exponent = Cones_Array[i].Generators_of_Cone[q].T_Exponent;
				Denominator_Current_Generator[0]->Children[1] = New_T_Node;
			
				// This generator is NOT a pole, thus do not factor out an 's'
				Exponent_Reduction_Offset = 0;
			}
			else //Pole, factor out an s 
			{
				
				Exponent_Reduction_Offset = 1;
			}

			
			// Expand the generator according to its type (1-rt) or (r-t)
			if (Cones_Array[i].Generators_of_Cone[q].Form_Type == ONE_SUB_RT) // means (1-rt)
			{
				Exponent_of_T = Cones_Array[i].Generators_of_Cone[q].T_Exponent;
				Numerical_Coefficient = -1;		
			}
			else // Form_Type is (r-t)
			{
				Exponent_of_T = 0;
				Numerical_Coefficient = 1;
			}
		
			Temp = Cones_Array[i].Generators_of_Cone[q].R_Exponent;

			// The "choose" function to calculate the coefficients
			for (int k = 1; k <= Cones_Array[i].order + Exponent_Reduction_Offset; k++)
			{
				Numerical_Coefficient *= Temp;
				Numerical_Coefficient /= k;

				New_T_Node = Controller.Get_T_Node ();
				New_T_Node->Node_Type = POLYTREE_T_NODE;
				New_T_Node->Coefficient = Numerical_Coefficient;
				New_T_Node->Exponent = Exponent_of_T;
				Denominator_Current_Generator[k- Exponent_Reduction_Offset] = New_T_Node;
			
				Temp--;
			}
			
			if ( OUTPUT_CONE == 1 )
			 {
			  	Simplify_Term << "(";

			    	for (int g = 0; g <= Cones_Array[i].order; g++)
			      	{
					Denominator_Current_Generator[g]->Print_Rational_Functions_to_File (Simplify_Term);
					if ( g != 0)
						Simplify_Term << "*s^" << g;
					if ( g != Cones_Array[i].order)
						Simplify_Term << "+";
			      	}
				if ( q+1 < noGsPerC)				
					Simplify_Term << ")*";
				else
					Simplify_Term << ")";
			 }
			if ( OUTPUT_CONE_MULTI == 1 )
			 {
			  	Simplify_Term << "(";

			    	for (int g = 0; g <= Cones_Array[i].order; g++)
			      	{
					Denominator_Current_Generator[g]->Print_Rational_Functions_to_File (Simplify_Term);
					if ( g != 0)
						Simplify_Term << "*s^" << g;
					if ( g != Cones_Array[i].order)
						Simplify_Term << "+";
			      	}
				if ( q+1 < noGsPerC)				
					Simplify_Term << ")*";
				else
					Simplify_Term << ")";
			 }
		

			
			if ( OUTPUT_CONE == 0 )
			{
			//   SHMUSHING BEGINS
			//   Shushing is the process of multiplying our polynomials for each generator but only retaining
			//   the coefficients for powers of s up to s^m, where m is the order of the pole for the cone
			
			for (int r = Cones_Array[i].order; r > 0; r--)  //Calculate the coeffiecient of s^order on down
			{
				Coefficient_Addition_Root = Controller.Get_PolyTree_Node ();
				Coefficient_Addition_Root->Number_of_Children = r + 1;
				Coefficient_Addition_Root->Node_Type = POLYTREE_ADD; // +
				
				for (int f = 0;  f < r + 1; f++)
				{
					Coefficient_Multiplication_Root = Controller.Get_PolyTree_Node ();
					Coefficient_Multiplication_Root->Node_Type = POLYTREE_MUL; // *
					Coefficient_Multiplication_Root->Number_of_Children = 2;
					
					Coefficient_Multiplication_Root->Children[0] = Denominator_Current_Generator[f];
					Coefficient_Multiplication_Root->Children[1] = Denominator_Result[r-f];		
				
					Coefficient_Addition_Root->Children[f] = Coefficient_Multiplication_Root;
				}

				Denominator_Result[r] = Coefficient_Addition_Root;
				
			}
			
			// Handle the s^0 coeffiecient
			
			Coefficient_Multiplication_Root = Controller.Get_PolyTree_Node ();
			Coefficient_Multiplication_Root->Node_Type = POLYTREE_MUL;
			Coefficient_Multiplication_Root->Number_of_Children = 2;
			
			Coefficient_Multiplication_Root->Children[0] = Denominator_Current_Generator[0];
			Coefficient_Multiplication_Root->Children[1] = Denominator_Result[0];

			Denominator_Result[0] = Coefficient_Multiplication_Root;
			
			} // end of if statment OUTPUT_CONE == 0
	
			// DEBUG
			/*
			cout << "Sign " << Cones_Array[i].sign << endl;
			for (int k = 0; k <= Cones_Array[i].order; k++)
			{
				Numerator_Vector[k]->Print ();
				if (k != 0)
					cout << "s^" << k ;
				
				if (k <= Cones_Array[i].order - 1)
					cout << " + " << endl;
			}	
			cout << endl << "  Divided by " << endl;


			for (int k = 0; k <= Cones_Array[i].order; k++)
			{
				Denominator_Result[k]->Print ();
				
				if (k != 0)
					cout << "s^" << k;
				
				if (k <= Cones_Array[i].order - 1)
					cout << " + " << endl;
			}
			
			cout << endl << endl;
			*/
	
		} //End of loop iterating through generators.	
			
			
		
		if (OUTPUT_CONE == 0)
		{
		
		// Coefficient recursion formula. Calculate the numerators N_k of the coefficients using
		// the recursion formula. 
		// First we rewrote the Coefficient formula to be:
		// C_k = (( 1/b_0 )^(k+1)) * N_k     where
		// N_k = ( (b_0^k)*(a_k) - (b_1)*N_(k-1) - (b_0)*(b_2)*N_(k-2) - (b_0^2)*(b_3)*N_(k-3) -...-
		//
		// This frees us from having to do polynomial division until the end when we represent
		// C_(order).  We get C_(order) by finding N_(order) and divide it by b_0^(order + 1) where 
		// N_(order) is simply a polynomial
		//
				
		//cout << "Residue coefficient" << endl;
			
		Quotient_Coefficient[0] = Numerator_Vector[0]; // C_0

		//Used for - sign in coefficient, and later with our Final_Result
		T_Node_Negative_One = Controller.Get_T_Node ();
		T_Node_Negative_One->Node_Type = POLYTREE_T_NODE;
		T_Node_Negative_One->Coefficient = -1;
		T_Node_Negative_One->Exponent = 0;

		for (int k = 1; k < Cones_Array[i].order + 1; k++)
		{
			// This will hold N_k
			Coefficient_Addition_Root = Controller.Get_PolyTree_Node ();
			Coefficient_Addition_Root->Node_Type = POLYTREE_ADD;
			Coefficient_Addition_Root->Number_of_Children = k + 1;
			
			// This holds the first term, (b_0^k)*(a_k)
			Coefficient_Multiplication_Root = Controller.Get_PolyTree_Node ();
			Coefficient_Multiplication_Root->Node_Type = POLYTREE_MUL;
			Coefficient_Multiplication_Root->Number_of_Children = 2;

			// This is (b_0^k)
			Coefficient_Exponent_Root = Controller.Get_PolyTree_Node ();
			Coefficient_Exponent_Root->Node_Type = POLYTREE_EXP;
			Coefficient_Exponent_Root->Number_of_Children = k;
			Coefficient_Exponent_Root->Children[0] = Denominator_Result[0];
				
			Coefficient_Multiplication_Root->Children[0] = Coefficient_Exponent_Root; // (b_0^k)
			Coefficient_Multiplication_Root->Children[1] = Numerator_Vector[k]; //a_k
				
			// Add (b_0^k)*(a_k) to our N_k
			Coefficient_Addition_Root->Children[0] = Coefficient_Multiplication_Root;
					
			//Do second term, -(b_1)(N_(k-1))
			Coefficient_Multiplication_Root = Controller.Get_PolyTree_Node ();
			Coefficient_Multiplication_Root->Node_Type = POLYTREE_MUL;
			Coefficient_Multiplication_Root->Number_of_Children = 3;

			Coefficient_Multiplication_Root->Children[0] = T_Node_Negative_One; // -1
			Coefficient_Multiplication_Root->Children[1] = Denominator_Result[1]; // b_1
			Coefficient_Multiplication_Root->Children[2] = Quotient_Coefficient[k-1]; // N_(k-1)
				
			//Add the term to our N_k
			Coefficient_Addition_Root->Children[1] = Coefficient_Multiplication_Root;

			//Do third term if we have a third term.  -(b_0)*(b_2)*(N_(k-2))
			if ( k > 1)
			{
				Coefficient_Multiplication_Root = Controller.Get_PolyTree_Node ();
				Coefficient_Multiplication_Root->Node_Type = POLYTREE_MUL;
				Coefficient_Multiplication_Root->Number_of_Children = 4;
				Coefficient_Multiplication_Root->Children[0] = T_Node_Negative_One; // -1
				Coefficient_Multiplication_Root->Children[1] = Denominator_Result[0];//b_0
				Coefficient_Multiplication_Root->Children[2] = Denominator_Result[2];//b_1
				Coefficient_Multiplication_Root->Children[3] = Quotient_Coefficient[k-2]; // N_(k-2)
				
				//Add the term to our N_k
				Coefficient_Addition_Root->Children[2] = Coefficient_Multiplication_Root;

			}
			
			//Do remaining terms according to the recursion relation
			for (int j = 3; j < k + 1; j++)
			{
				// Each term is a multiplication of 4 things, one of them being negative one
				Coefficient_Multiplication_Root = Controller.Get_PolyTree_Node ();
				Coefficient_Multiplication_Root->Node_Type = POLYTREE_MUL;
				Coefficient_Multiplication_Root->Number_of_Children = 4;
				Coefficient_Multiplication_Root->Children[0] = T_Node_Negative_One;
					
				// This is (b_0)^(j-1)
				Coefficient_Exponent_Root = Controller.Get_PolyTree_Node ();
				Coefficient_Exponent_Root->Node_Type = POLYTREE_EXP;
				Coefficient_Exponent_Root->Number_of_Children = j-1;
				Coefficient_Exponent_Root->Children[0] = Denominator_Result[0]; // b_0

				Coefficient_Multiplication_Root->Children[1] = Coefficient_Exponent_Root;//b_0^(j-1)
				Coefficient_Multiplication_Root->Children[2] = Denominator_Result[j];//b_j
				Coefficient_Multiplication_Root->Children[3] = Quotient_Coefficient[k-j];//N_(k-j)

				//Add the term to our N_k
				Coefficient_Addition_Root->Children[j] = Coefficient_Multiplication_Root;
					
			}		

			// Now Quotient_Coefficient[k] holds N_k
			Quotient_Coefficient[k] = Coefficient_Addition_Root;
		
				
				
		} //End of for loop going from 0 to "order"			
	
		} //End of if statement for OUTPUT_CONE == 0
		
		
		if ( OUTPUT_CONE == 1)
		{
			Simplify_Term << "):";

			Simplify_Term.close ();

			system_with_error_check ("maple <simplify3.add >out.simplify");
		
			cout << "%";	
		}
		
		if ( OUTPUT_CONE_MULTI == 1)
		{
			Simplify_Term << "):" << endl;

			if (Cones_File_Count == Max_Cones_File)
			{
				Cones_File_Count = 0;

				system_with_error_check ("maple <simplify4.add >out.simplify");

			}	
			

		}
		
		if (OUTPUT_S_T_FUNCTION == 1)
			{
				Rational_Function_S_T.open ("simplify.term");
			
				Rational_Function_S_T << "d :=  " << Cones_Array[i].order << ":" <<endl;
				
				if (Cones_Array[i].sign == 1)
					Rational_Function_S_T << "x := (";
				else
					Rational_Function_S_T << "x := (-1)*(";
				
				for (int g=0; g <= Cones_Array[i].order; g++)
				{
					Numerator_Vector[g]->Print_Rational_Functions_to_File ( Rational_Function_S_T );
					if (g != 0)
						Rational_Function_S_T << "*s^" << g;
					if (g != Cones_Array[i].order)
						Rational_Function_S_T << "+";
				}
					
				Rational_Function_S_T << ")/(";
					
				for (int g=0; g <= Cones_Array[i].order; g++)
				{
					Denominator_Result[g]->Print_Rational_Functions_to_File ( Rational_Function_S_T );
					if (g != 0)
						Rational_Function_S_T << "*s^" << g;
					if (g != Cones_Array[i].order)
						Rational_Function_S_T << "+";
				}

				Rational_Function_S_T << "):" << endl;
				
				Rational_Function_S_T.close ();	

				system_with_error_check ("maple < simplify2.add >out.simplify");
			}

		//Now we have all the N_k's up to N_(order)
		
		//Now we create a division root to hold C_(order) = ( (1/b_0)^(order+1) )*N_(order)
		Quotient_Root = Controller.Get_PolyTree_Node ();
		Quotient_Root->Node_Type = POLYTREE_DIV;
		Quotient_Root->Number_of_Children = 2;
		Quotient_Root->Children[0] = Quotient_Coefficient[Cones_Array[i].order]; //N_order
			
		// This will hold (b_0)^(order + 1)
		Coefficient_Exponent_Root = Controller.Get_PolyTree_Node ();
		Coefficient_Exponent_Root->Node_Type = POLYTREE_EXP;
		Coefficient_Exponent_Root->Number_of_Children = Cones_Array[i].order + 1;
		Coefficient_Exponent_Root->Children[0] = Denominator_Result[0]; //b_0

		Quotient_Root->Children[1] = Coefficient_Exponent_Root; // (b_0)^(order + 1)

			
		if (Cones_Array[i].sign == 1)
			Final_Cone_Expression = Quotient_Root;
		else
		{
			Coefficient_Multiplication_Root = Controller.Get_PolyTree_Node ();
			Coefficient_Multiplication_Root->Node_Type = POLYTREE_MUL;
			Coefficient_Multiplication_Root->Number_of_Children = 2;

			Coefficient_Multiplication_Root->Children[0] = T_Node_Negative_One;
			Coefficient_Multiplication_Root->Children[1] = Quotient_Root;		
			Final_Cone_Expression = Coefficient_Multiplication_Root;
		}	
		
		if (COMPUTE_SUM_MAPLE == 1)
		{
			//system_with_error_check("rm simplify.term");
			cout << "%";
			
			Simplify_Term.open ("simplify.term");

			Simplify_Term << "x :=";
			Final_Cone_Expression->Print_Rational_Functions_to_File ( Simplify_Term );
			Simplify_Term << ":" << endl;
			
			Simplify_Term.close ();

			system_with_error_check ("maple <simplify.add >out.simplify");
			
		}
		
		if(print_flag == 1)
		{
		  if(never_printed)
		    {
		      never_printed = 0;
		      Rational_Function_Output_File << "x := ";
		    }
                  else
                      	Rational_Function_Output_File << " + ";
			Final_Cone_Expression->Print_Rational_Functions_to_File( Rational_Function_Output_File );
			//	Rational_Function_Output_File << ":" << endl;
		}
		
		if (OUTPUT_CONE == 0 )
		{
		
		Final_Cone_Expression->Taylor_Expansion(Cone_Taylor_Parameters);	

		//cout << "Residue: Taylor Expansion of " << i << " cone. ";
		
		//for (int k = 0; k <= DEGREE; k++)
		//{
		//	cout << Cone_Taylor_Parameters->Result[k] << "t^" << k << " + ";
		//}
		//cout << endl;	
		
		for (int k = 0; k <= DEGREE; k++)
		{
			Final_Taylor_Result[k] += Cone_Taylor_Parameters->Result[k];
		}
		
		}
		
		// Reset all the PolyTree_Node and T_Nodes to be reused to save memory :)
		Controller.Reset ();
		

		delete [] Cones_Array[i].Generators_of_Cone;
		
		if (i%50 == 0)
			cout << i << " / " << noCones << " Done " << endl;
		
	} //End of for loop iterating through all the cones

	
	//Output the rational functions to file

	if(print_flag == 1)
	{
		Rational_Function_Output_File << ":" << endl;
		Rational_Function_Output_File.close ();	 
	}
	
	if(DEGREE > 1){
	cout << endl << "This is the taylor expansion of our signed sum of rational functions up to degree ";
	cout << DEGREE << endl << endl;
	
	// The quotient function internaly will multiply its results by Ten_Power
	// So to print out the correct number, divide by Ten_Power
	
	for (int i = 0; i <= DEGREE; i++)
	{
		//Coefficient_Total += (Final_Parameters->Result[i] + *Ten_Power/2) / *Ten_Power;
		cout << (Final_Taylor_Result[i] + *Ten_Power/2) /  *Ten_Power; 
		if ( i != 0)
			cout << "t^" << i;
		if ( i != DEGREE)
			cout << " + ";
	}

	cout << endl << endl; }
	else if(DEGREE == 1){

	  // cout << "\n****  Total number of lattice points is: " << (Final_Taylor_Result[1] + *Ten_Power/2) /  *Ten_Power << "  ****" << endl << endl;
	  ofstream out("numOfLatticePoints");
	  out << (Final_Taylor_Result[1] + *Ten_Power/2) /  *Ten_Power << endl;
	}
	
	delete[] Cone_Taylor_Parameters->Result;
	delete	Cone_Taylor_Parameters;	
  	delete	Ten_Power;
	
	return ;
}

