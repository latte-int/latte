#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include <stdlib.h>
#include <list>
//#include <vector>
#include "PolyTree.h"
#include "myheader.h"
#include "cone.h"
#include "ramon.h"
#include "RudyResNTL.h"
#include <malloc.h>
#include <stdlib.h>
#include "flags.h"

int ResidueFunction_Single_Cone(listCone *cones,
				Standard_Single_Cone_Parameters *Residue_Parameters)
{
  Node_Controller *Controller = Residue_Parameters->Controller;
  	listCone *C, * cones1;
  	int OUTPUT_S_T_FUNCTION = 0, COMPUTE_SUM_MAPLE = 0, OUTPUT_CONE = 0;
  	int DEGREE;
	int dim, noGsPerC,noCones; //noGsPerC is number of generators per cone
  	if(Residue_Parameters->Degree_of_Taylor_Expansion > 1) 
		DEGREE = Residue_Parameters->Degree_of_Taylor_Expansion;
	else
		DEGREE = 1;
 	if( ((Residue_Parameters->Flags & OUTPUT) >> 1) == 1 ) 
		OUTPUT_CONE = 1;
 	else if( ((Residue_Parameters->Flags & OUTPUT) >> 1) == 2) 
		COMPUTE_SUM_MAPLE = 1;
 	else if( ((Residue_Parameters->Flags & OUTPUT) >> 1) == 3) 
		OUTPUT_S_T_FUNCTION = 1;

  	int numOfTerms=0;

  	C=cones;
  	while (C) 
	{
    		numOfTerms++;
    		C=C->rest;
  	}



  	dim=Residue_Parameters->Number_of_Variables;
  	noGsPerC=lengthListVector(cones->rays);
  	noCones=numOfTerms;

  	Cone_Data	Cones_Array[noCones];  // Create structure to hold all the information of the cones.
  
  	for (int q = 0;q < noCones; q++)  // Do some initialization
  	{
		Cones_Array[q].order = 0;
		Cones_Array[q].Generators_of_Cone = new Generator [noGsPerC];  //Each cone has noGsPerC generators.
  	}
  
	// Added by P/D
  	int i,j; // index or loop vars
  	long int k;//n=0,p; // extra vars to use as needed
  	//long int totalNoGs=noGsPerC*noCones; //total no. of generators,ie,rowdim of B
  	// long int B[totalNoGs][dim];  // B is the denominator vectors
  	//  cout<<"tNG: "<<totalNoGs<<endl;
	int result = 1;
	
  	Integer tmp_A;

  	listVector *basis;
  	listCone *listtmp3;
  	cones1 = cones;
  	i = 0;
  	while (cones1) 
	{
    		//tmp=cones1->latticePoints;
    		//while (tmp) 
		//{
      			Cones_Array[i].sign = cones1->coefficient;   //Added by P/D:  Load what sign each cone is
      		
     			basis = cones1->rays;
     			while(basis) 
			{
				//cout << "R_Exponent";
    				for (j=0; j<noGsPerC; j++) 
				{
					Cones_Array[i].Generators_of_Cone[j].R_Exponent = 0;
					
					Cones_Array[i].Generators_of_Cone[j].T_Exponent=basis->first[dim -1];  //exponent of t in the denominator	
    					//cout << "ResidueFunction_Single: Random Lambda = ";
					for(k = 0; k < dim - 1; k++)
					{
    						Cones_Array[i].Generators_of_Cone[j].R_Exponent += basis->first[k] * Residue_Parameters->generic_vector[k];
      						//cout << Residue_Parameters->generic_vector[k] << " ";
					}
					//cout << endl;

					// Test to see if dot product is zero...if so, barf
					if(Cones_Array[i].Generators_of_Cone[j].R_Exponent == 0 && Cones_Array[i].Generators_of_Cone[j].T_Exponent == 0)
					{
						cout << "ResidueFunction_Single: zero dotproduct. ";
						for (int p = 0; p < dim-1;p++)
							cout << basis->first[p] << " ";
						cout << endl;
						result = -1;	
      					
					}
					//else	
					//{
					//	for (int p = 0; p < dim;p++)
					//		cout << basis->first[p] << " ";
					//	cout << endl;
					//}
					//cout << Cones_Array[i].Generators_of_Cone[j].R_Exponent << " ";	
					// Check when k=dim-1 if input is 0, increment order of cone
					//Cones_Array[i].Generators_of_Cone[j].T_Exponent=basis->first[dim -1];  //exponent of t in the denominator	
					Cones_Array[i].Generators_of_Cone[j].Form_Type = ONE_SUB_RT;
					if (basis->first[dim - 1] == 0)  //if the exponent of t is zero, increment the order for this cone
					Cones_Array[i].order++;
      
    					basis = basis->rest;
    				}
  			}
    			//listtmp2 = tmp;
      			//tmp=tmp->rest; i++;
      			//delete listtmp2;
    		//}
    		listtmp3 = cones1;
    		cones1 = cones1->rest;
    		freeCone(listtmp3);
  	}
	//cout << "ResidueFunction_Single: Done reading input." << endl;  
	
	if (result == -1)
		return result;
	
	i = 0;
 	
	//cout << "Copying dot product into Cones_Array" << endl; 
  	//*****************************************************************
  	//  PETER/DAVE CODE REALLY BEGINS HERE
  	//  Our data structure Cones_Array which holds all of our information
  	//  is currently loaded with the sign of each cone, the
  	//  exponent of t in the numerators, and the exponent of t
  	//  for all the generators.

  	
	//  CALCLUTE DOT PRODUCT OF THE NUMERATOR (= 0) AND STORE IN CONES_ARRAY
  	for (int q = 0;q < noCones;q++)  // For each cone
	{
		Cones_Array[q].Numerator_Generator.R_Exponent = 0;
		Cones_Array[q].Numerator_Generator.T_Exponent = 0;	
		
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
	//for (int q = 0; q < noCones; q++)
	//{
	//	Cones_Array[q].Numerator_Generator.R_Exponent -= Numerator_R_Exponent_Minimum;
	//}


	// Create all the variables we are going to use in our big loop
	
	PolyTree_Node	*Numerator_Vector[noGsPerC + 1]; // [i] - coefficient of s^i
	PolyTree_Node 	*Denominator_Result[noGsPerC + 1];   // Used to store the running total of generators
	PolyTree_Node 	*Denominator_Current_Generator[noGsPerC + 1];  //For each generatore as we iterate through them
	PolyTree_Node	*Quotient_Coefficient[noGsPerC + 1];  // holds the coefficients for calculting residue
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

	T_Node_Negative_One = Controller->Get_T_Node ();
	T_Node_Negative_One->Node_Type = POLYTREE_T_NODE;
	T_Node_Negative_One->Coefficient = -1;
	T_Node_Negative_One->Exponent = 0;

	
	ZZ	Numerical_Coefficient; //Used for calculating coefficient of binomial expansion of (s+1)^p
	ZZ	Temp;  // Used for calculation coefficient
	ZZ	Exponent_of_T;  //Form_Type dictates what this will be

	int	Exponent_Reduction_Offset;  //Used when we factor out 1/s^(order) of each cone
					// It is equal to 0 or 1 for each generator depending
					// if it is a pole or not 
	
	// Initialize our final expression

	

	
	
	Taylor_Parameters *Cone_Taylor_Parameters = new Taylor_Parameters;
	
	//Calculate the taylor expasion for the first DEGREE + 1 many terms
	Cone_Taylor_Parameters->Result = new ZZ [DEGREE + 1];
	Cone_Taylor_Parameters->Ten_Power = &Residue_Parameters->Ten_Power;
	Cone_Taylor_Parameters->Degree_of_Expansion = DEGREE;
	ofstream Rational_Function_Output_File;
	
	for (int i = 0;i <= DEGREE; i++)
		Cone_Taylor_Parameters->Result[i] = 0;
		
	if(Residue_Parameters->Flags & PRINT	== 1)
		 {
		   // system ("rm func.rat");
	 		 //cout << "Outputing rational functions to file" << endl;
		   // char File[200];
		   // strcpy(File, fileName);
		   // strcat(File, ".rat");
		    Rational_Function_Output_File.open ("func.rat");
		 }

	//cout << "Formulating rational functions and performing taylor expansion on cones." << endl;

	ofstream Simplify_Sum, Simplify_Term;
	
	if (COMPUTE_SUM_MAPLE == 1)
	{
		//system ("rm simplify.sum");
		
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
	
	// For each Cone
	for (int i = 0; i < noCones; i++)
	{
		//Calculate the numerator vector
		New_T_Node = Controller->Get_T_Node ();
		
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

			New_T_Node = Controller->Get_T_Node ();
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
		
		// Initialize the Denominator result with the information of the first generator, 
		// (only if the T_Exponent in non-zero, i.e. this generator is not a pole)

		
		
		if(Cones_Array[i].Generators_of_Cone[0].T_Exponent != 0)
		{
			Denominator_Result[0] = Controller->Get_PolyTree_Node ();
			Denominator_Result[0]->Node_Type = POLYTREE_ADD; // +
			Denominator_Result[0]->Number_of_Children = 2;
		
			New_T_Node = Controller->Get_T_Node ();
			New_T_Node->Node_Type = POLYTREE_T_NODE;
			New_T_Node->Coefficient = 1;
			New_T_Node->Exponent = 0;
			
			Denominator_Result[0]->Children[0] = New_T_Node;

			New_T_Node = Controller->Get_T_Node ();
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

			New_T_Node = Controller->Get_T_Node ();
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
				Denominator_Current_Generator[0] = Controller->Get_PolyTree_Node ();
				Denominator_Current_Generator[0]->Node_Type = POLYTREE_ADD; // +
				Denominator_Current_Generator[0]->Number_of_Children = 2;
		
				New_T_Node = Controller->Get_T_Node ();
				New_T_Node->Node_Type = POLYTREE_T_NODE;
				New_T_Node->Coefficient = 1;
				New_T_Node->Exponent = 0;
				Denominator_Current_Generator[0]->Children[0] = New_T_Node;

				New_T_Node = Controller->Get_T_Node ();
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

				New_T_Node = Controller->Get_T_Node ();
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

		

			
			if ( OUTPUT_CONE == 0 )
			{
			//   SHMUSHING BEGINS
			//   Shushing is the process of multiplying our polynomials for each generator but only retaining
			//   the coefficients for powers of s up to s^m, where m is the order of the pole for the cone
			
			for (int r = Cones_Array[i].order; r > 0; r--)  //Calculate the coeffiecient of s^order on down
			{
				Coefficient_Addition_Root = Controller->Get_PolyTree_Node ();
				Coefficient_Addition_Root->Number_of_Children = r + 1;
				Coefficient_Addition_Root->Node_Type = POLYTREE_ADD; // +
				
				for (int f = 0;  f < r + 1; f++)
				{
					Coefficient_Multiplication_Root = Controller->Get_PolyTree_Node ();
					Coefficient_Multiplication_Root->Node_Type = POLYTREE_MUL; // *
					Coefficient_Multiplication_Root->Number_of_Children = 2;
					
					Coefficient_Multiplication_Root->Children[0] = Denominator_Current_Generator[f];
					Coefficient_Multiplication_Root->Children[1] = Denominator_Result[r-f];		
				
					Coefficient_Addition_Root->Children[f] = Coefficient_Multiplication_Root;
				}

				Denominator_Result[r] = Coefficient_Addition_Root;
				
			}
			
			// Handle the s^0 coeffiecient
			
			Coefficient_Multiplication_Root = Controller->Get_PolyTree_Node ();
			Coefficient_Multiplication_Root->Node_Type = POLYTREE_MUL;
			Coefficient_Multiplication_Root->Number_of_Children = 2;
			
			Coefficient_Multiplication_Root->Children[0] = Denominator_Current_Generator[0];
			Coefficient_Multiplication_Root->Children[1] = Denominator_Result[0];

			Denominator_Result[0] = Coefficient_Multiplication_Root;
			
			} // end of if statment OUTPUT_CONE == 0
	
	
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

		for (int k = 1; k < Cones_Array[i].order + 1; k++)
		{
			// This will hold N_k
			Coefficient_Addition_Root = Controller->Get_PolyTree_Node ();
			Coefficient_Addition_Root->Node_Type = POLYTREE_ADD;
			Coefficient_Addition_Root->Number_of_Children = k + 1;
			
			// This holds the first term, (b_0^k)*(a_k)
			Coefficient_Multiplication_Root = Controller->Get_PolyTree_Node ();
			Coefficient_Multiplication_Root->Node_Type = POLYTREE_MUL;
			Coefficient_Multiplication_Root->Number_of_Children = 2;

			// This is (b_0^k)
			Coefficient_Exponent_Root = Controller->Get_PolyTree_Node ();
			Coefficient_Exponent_Root->Node_Type = POLYTREE_EXP;
			Coefficient_Exponent_Root->Number_of_Children = k;
			Coefficient_Exponent_Root->Children[0] = Denominator_Result[0];
				
			Coefficient_Multiplication_Root->Children[0] = Coefficient_Exponent_Root; // (b_0^k)
			Coefficient_Multiplication_Root->Children[1] = Numerator_Vector[k]; //a_k
				
			// Add (b_0^k)*(a_k) to our N_k
			Coefficient_Addition_Root->Children[0] = Coefficient_Multiplication_Root;
					
			//Do second term, -(b_1)(N_(k-1))
			Coefficient_Multiplication_Root = Controller->Get_PolyTree_Node ();
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
				Coefficient_Multiplication_Root = Controller->Get_PolyTree_Node ();
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
				Coefficient_Multiplication_Root = Controller->Get_PolyTree_Node ();
				Coefficient_Multiplication_Root->Node_Type = POLYTREE_MUL;
				Coefficient_Multiplication_Root->Number_of_Children = 4;
				Coefficient_Multiplication_Root->Children[0] = T_Node_Negative_One;
					
				// This is (b_0)^(j-1)
				Coefficient_Exponent_Root = Controller->Get_PolyTree_Node ();
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

			system ("maple <simplify3.add >out.simplify");
		
			cout << "%";	
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

				system ("maple < simplify2.add >out.simplify");
			}

		//Now we have all the N_k's up to N_(order)
		
		//Now we create a division root to hold C_(order) = ( (1/b_0)^(order+1) )*N_(order)
		Quotient_Root = Controller->Get_PolyTree_Node ();
		Quotient_Root->Node_Type = POLYTREE_DIV;
		Quotient_Root->Number_of_Children = 2;
		Quotient_Root->Children[0] = Quotient_Coefficient[Cones_Array[i].order]; //N_order
			
		// This will hold (b_0)^(order + 1)
		Coefficient_Exponent_Root = Controller->Get_PolyTree_Node ();
		Coefficient_Exponent_Root->Node_Type = POLYTREE_EXP;
		Coefficient_Exponent_Root->Number_of_Children = Cones_Array[i].order + 1;
		Coefficient_Exponent_Root->Children[0] = Denominator_Result[0]; //b_0

		Quotient_Root->Children[1] = Coefficient_Exponent_Root; // (b_0)^(order + 1)

			
		if (Cones_Array[i].sign == 1)
			Final_Cone_Expression = Quotient_Root;
		else
		{
			Coefficient_Multiplication_Root = Controller->Get_PolyTree_Node ();
			Coefficient_Multiplication_Root->Node_Type = POLYTREE_MUL;
			Coefficient_Multiplication_Root->Number_of_Children = 2;

			Coefficient_Multiplication_Root->Children[0] = T_Node_Negative_One;
			Coefficient_Multiplication_Root->Children[1] = Quotient_Root;		
			Final_Cone_Expression = Coefficient_Multiplication_Root;
		}	
		
		if (COMPUTE_SUM_MAPLE == 1)
		{
			//system ("rm simplify.term");
			cout << "%";
			
			Simplify_Term.open ("simplify.term");

			Simplify_Term << "x :=";
			Final_Cone_Expression->Print_Rational_Functions_to_File ( Simplify_Term );
			Simplify_Term << ":" << endl;
			
			Simplify_Term.close ();

			system ("maple <simplify.add >out.simplify");
			
		}
		
		if(Residue_Parameters->Flags & PRINT == 1)
		{
			Rational_Function_Output_File << "x := ";
			Final_Cone_Expression->Print_Rational_Functions_to_File( Rational_Function_Output_File );
			Rational_Function_Output_File << ":" << endl;
		}
		
		if (OUTPUT_CONE == 0 )
		{
		
		//cout << "ResidueFunction: Computer taylor expansion...";	
		Final_Cone_Expression->Taylor_Expansion(Cone_Taylor_Parameters);	
		//cout << "done." << endl;
	
		//cout << "ResidueFunction: Copying result of taylor expansion into Parameters result...";
		for (int k = 0; k <= DEGREE; k++)
		{
			Residue_Parameters->Taylor_Expansion_Result[k] += Cone_Taylor_Parameters->Result[k];
		}
		//cout << "done." << endl;
		
		}
		
		// Reset all the PolyTree_Node and T_Nodes to be reused to save memory :)
		Controller->Reset ();
		

		delete [] Cones_Array[i].Generators_of_Cone;
		
		
	} //End of for loop iterating through all the cones

	
	//Output the rational functions to file

	if(Residue_Parameters->Flags & PRINT == 1)
	{
			 
		Rational_Function_Output_File.close ();	 
	}
	
	for (int k = 0; k <= DEGREE;k++)
		Cone_Taylor_Parameters->Result[k].kill ();
	
	delete[] Cone_Taylor_Parameters->Result;
	delete	Cone_Taylor_Parameters;	
	
	return 1;
}

