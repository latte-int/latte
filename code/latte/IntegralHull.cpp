#include <string.h>
#include <stdio.h>
#include <string>
#include <math.h>
#include <fstream>
#include <ctype.h>

#include "config.h"
#include "myheader.h"
#include "ramon.h"
#include "print.h"
#include "IntegralHull.h"
#include "cone.h"
#include "rational.h"
#include "ConeInfo.h"

using namespace std;
int IntegralHull_Flag = 1;
int Verify_IP_Flag = 1;

ZZ	Polytope_Max_Width;


#define BOUND 150

listCone* FindRationalFunction(listCone* cones, vec_ZZ a, vec_ZZ cost, int numOfVars)
{
  	listVector * endRays, *rays, *tmpRays;
  	vec_ZZ Vertex, OptVertex, numerator, v;
  	numerator.SetLength(numOfVars);
  	ZZ Opt, tmpDotProd;
  	Opt = -10000000;
  	//  int flag = 0;
  	listCone * Tmpcones, *tmp, *tmpcone, *endcone;
  	listVector* tmpVector;
  	tmpcone = cones;
  	int number = 0, Coeff;
  	list< vec_ZZ > solutions;
  	Tmpcones =createListCone();  
  	endcone = Tmpcones;

  	while(tmpcone)
	{
   		tmpDotProd = 0;
    		number = 0;
    		tmpVector = tmpcone -> rays;

    		numerator = tmpcone ->latticePoints -> first;

    		Coeff = tmpcone -> coefficient;

    		while(tmpVector)
		{
      
		/*      if(((tmpVector -> first)*cost) == 0){
		cout <<"Zero dot product." << endl; exit (0);
		}*/
      		if(((tmpVector -> first)*cost) > 0)
		{ 
			number++;
      			numerator -= tmpVector -> first;

      		}
      		else ;
      		tmpVector = tmpVector -> rest;
   		}

    		if(numerator == a)
		{
  			tmp =createListCone();

  			tmp->vertex=createRationalVector(numOfVars);
  			rays=createListVector(createVector(numOfVars));
  			tmp ->latticePoints=createListVector(createVector(numOfVars));
  			endRays=rays;
  			//  tmp -> rays = tmpcone -> rays;
  
   			tmp ->latticePoints -> first = a;
  
   			tmp -> coefficient = Coeff;
   			tmp -> vertex = createRationalVector(numOfVars);
   			tmp -> vertex = tmpcone -> vertex;
   			tmpRays=tmpcone -> rays;
     			while (tmpRays) 
			{
      				v=createVector(numOfVars);
      				for (int i=0; i<numOfVars; i++) 
					v[i]=(tmpRays->first)[i];
      				endRays->rest=createListVector(v);
      				endRays=endRays->rest;
      				tmpRays=tmpRays->rest;
     			} 
   			tmp->rays=rays->rest;

    			endcone -> rest = tmp;
    			endcone = endcone -> rest;
    		}
    	tmpcone = tmpcone -> rest;
  	}
//  printf("=======================\n"); 
//   printListCone(Tmpcones,numOfVars);
  

  	return Tmpcones -> rest;
} // FindRationalFunction

vec_ZZ SolveIP(listCone* cones,listVector* matrix, vec_ZZ cost, int numOfVars, int SINGLE_CONE)
{
  //int SINGLE_CONE = 1;
	//cout << "SolveIP Called. Cost = " << cost << endl;
  	int	Digging_Count = 0;
	vec_ZZ OptVertex;
  	int  flag = 0;
  	listCone * tmpcone;
  	listVector  *tmpmatrix;
  	tmpcone = cones;
  	list< vec_ZZ > solutions;
  	vec_ZZ		Temp_Vector;
  
	vec_ZZ possible;
  	ZZ RHS;
  
	//  cout << cost << endl;
	//

	ConeInfo *new_info;

	/**************
	// CREATE THE CONE HEAP!!!!
	***************/
	
	ConeInfo_Heap cone_heap;
	
	vec_ZZ Our_Cost;

	Our_Cost = cost;
	
	int	Pertubation_Count = 0;	

	int	Cone_Heap_Count = 0;	
	//cout << "Solve_IP: Creating Cone_Heap...";
	while(tmpcone)
	{

		new_info = new ConeInfo (&Our_Cost, tmpcone, numOfVars);
			
		/*******************
		// THROW THE CONE_INFO IN THE HEAP!
		********************/
		cone_heap.Add_Heap(new_info);
		
		tmpcone = tmpcone->rest;

		if (new_info->S_Values_Zero_Flag == 1)
		{
			Pertubation_Count++;
			if(SINGLE_CONE == 1) {

			  cerr << "Zero Dot product.  Please start IP without a single cone method." << endl;
			  exit(1);
			}
			cout << "S_Value zero for some cone. Pertubating. " << endl;
			
			//if ((Pertubation_Count % 1) == 0)
			//       cout << "%";	
			cone_heap.Clear_Tree ();

			tmpcone = cones;
			//ZZ 	Normalize_Length;
			//Normalize_Length = 10000;			
			//Normalize_Length = Calculate_Polytope_Width (cones, matrix, numOfVars);
			
			Our_Cost = Calculate_Pertubation (cones, &cost, 10, numOfVars);			
			Cone_Heap_Count = 0;
		}	
			
		//if ((Cone_Heap_Count % 100) == 0)
		//	cout << Cone_Heap_Count << " added to Cone_Heap. " << endl;
		Cone_Heap_Count++;
  	}
	//cout << "Done" << endl;
	//cout << endl;
	//cout << "SolveIP: cone_heap ready. Looping until coefficient_Sum is nonzero." << endl;
	
	/***************************
	// LOOP UNTIL OPTIMAL VALUE IS FOUND
	******************************/

	ConeInfo_List	*ConeInfo_List_Highest_Terms, *Head_ConeInfo_List, *Temp_ConeInfo_List;
	int	Coefficient_Sum;

	// This loop only breaks when the coefficient_Sum is nonzero
	while(1)
	{
		//cout << "$% ";
		// pop items off coneinfo heap as long as they have the same
		// current highest term as the first item removed
		
		ConeInfo_List_Highest_Terms = new ConeInfo_List;
		Head_ConeInfo_List = ConeInfo_List_Highest_Terms;
		
		
		ConeInfo_List_Highest_Terms->Next = NULL;

		// pop coneinfo off the top of the heap
		ConeInfo_List_Highest_Terms->ConeInfo_Pointer = cone_heap.Pop_Top_Heap ();

		
		Coefficient_Sum = ConeInfo_List_Highest_Terms->ConeInfo_Pointer->Get_Coefficient ();
		//cout << "First cone coefficient: " << Coefficient_Sum << " " << ConeInfo_List_Highest_Terms->ConeInfo_Pointer->Heap << endl;
		
		//cout << "SolveIP: Highest Term removed, checking for more on heap." << endl;	
		
		while (cone_heap.Check_Top_Heap ( ConeInfo_List_Highest_Terms->ConeInfo_Pointer ) == 1 )
		{
			ConeInfo_List_Highest_Terms->Next = new ConeInfo_List;
			ConeInfo_List_Highest_Terms->Next->ConeInfo_Pointer = cone_heap.Pop_Top_Heap ();
			ConeInfo_List_Highest_Terms->Next->Next = NULL;
			
			Coefficient_Sum += ConeInfo_List_Highest_Terms->Next->ConeInfo_Pointer->Get_Coefficient ();

			//cout << ConeInfo_List_Highest_Terms->Next->ConeInfo_Pointer->Get_Coefficient () << " ";
			//cout << ConeInfo_List_Highest_Terms->Next->ConeInfo_Pointer->Heap << " " << endl;
			ConeInfo_List_Highest_Terms = ConeInfo_List_Highest_Terms->Next;
			
		}
		//cout << endl;	

		//cout << "SolveIP: Checking coefficient_sum. " << endl;	
			
		// if Coefficnt_Sum is 0 we got a problem !!!

		if (Coefficient_Sum == 0)
		{
			if (Digging_Count % 100 == 0)
				cout << "Digging [" << Digging_Count + 1 << "] ";

			Digging_Count++;	
			ConeInfo_List_Highest_Terms = Head_ConeInfo_List;

			while (ConeInfo_List_Highest_Terms != NULL)
			{
				//cout << "Digging: calculating next term: ";
				ConeInfo_List_Highest_Terms->ConeInfo_Pointer->Calculate_Next_Highest_Term ();
				//cout << "Done calculation next highest term." << endl;
			
				cone_heap.Add_Heap (ConeInfo_List_Highest_Terms->ConeInfo_Pointer);
				
				Temp_ConeInfo_List = ConeInfo_List_Highest_Terms;
					
				ConeInfo_List_Highest_Terms = ConeInfo_List_Highest_Terms->Next;
			
				delete Temp_ConeInfo_List;
			}
			//cout << "Done Digging " << Digging_Count << endl;
		}
		else	//if coefficient_sum != 0 then break
		{	
			// Check if feasible, if not dig!
			// Only do if flag is set.
			flag = 0;
			//if (SINGLE_CONE == 1)
			//{
				list< vec_ZZ > solutions2;

				Temp_ConeInfo_List = Head_ConeInfo_List;
				Temp_Vector.SetLength(numOfVars);

				while (Temp_ConeInfo_List != NULL)
				{
					while (Temp_ConeInfo_List->ConeInfo_Pointer->Calculate_Integral_Point (Temp_Vector) == 1 )
					{
						solutions2.push_front(Temp_Vector);	
						//solutions.push_front(Temp_Vector);
					}
		
					Temp_ConeInfo_List = Temp_ConeInfo_List->Next;
				}	
				
				OptVertex.SetLength(numOfVars);
				
				while(!solutions2.empty())
				{
    					flag = 0;
    					possible = solutions2.front(); 
    					solutions2.pop_front(); //cout << possible << endl;
    					//cout << possible << " ";
					tmpmatrix = matrix;      
    					while(tmpmatrix)
					{ 
						RHS = 0; //cout <<"after while: " << possible << endl;
      						for(int i = 0; i < numOfVars; i++)
							RHS += (tmpmatrix -> first)[i + 1] * possible[i];
      						if((tmpmatrix -> first[0] + RHS) < 0)
						{
							// flag = 1 means infeasible
							flag = 1;
          						break;
						}
      						else ;
      						tmpmatrix = tmpmatrix -> rest;
    					}
    					if(flag == 0)
					{	
						OptVertex = possible; 
						//solutions.clear ();
						//possible = solutions2.front ();
						//cout << "Solution feasible, breaking." << endl;
						break;
					}
					if(solutions2.empty())
					{
						if (SINGLE_CONE != 1)
						{
							cout << "Not in SINGLE_CONE mode, coefficient is non zero but no feasibible solutions" << endl;
							exit (1);
						}
						// solution set is empty, so make them dig!
						cout << "Solutions is empty.  Point not feasible, make em dig!";
						ConeInfo_List_Highest_Terms = Head_ConeInfo_List;

						while (ConeInfo_List_Highest_Terms != NULL)
						{
							//cout << "Digging: calculating next term: ";
							ConeInfo_List_Highest_Terms->ConeInfo_Pointer->Calculate_Next_Highest_Term ();
							//cout << "Done calculation next highest term." << endl;
			
							cone_heap.Add_Heap (ConeInfo_List_Highest_Terms->ConeInfo_Pointer);
				
							Temp_ConeInfo_List = ConeInfo_List_Highest_Terms;
					
							ConeInfo_List_Highest_Terms = ConeInfo_List_Highest_Terms->Next;
			
							delete Temp_ConeInfo_List;
						}
						flag = 1;	
					}
				}
		
			//}
			
			// break if coefficient is not zero and if in non-SINGLE_CONE mode solution feasible	
			
			if (flag == 0) 
				break;
		}
		
	}
	if(Pertubation_Count == 0){
	  if(Coefficient_Sum == 1)
	    cout << endl << "There is one optimal solution. \t\t" << endl;
	  else
	    cout << endl << "There are " << Coefficient_Sum << " optimal solutions.\t\t" << endl;
	  //cout << "SolveIP: coefficient_sum is nonzero, exit while loop." << endl;
	}
	
	//if (SINGLE_CONE != 1)
	//{
	
	// feasiblity testing is now inside the loop!
	/*
	Temp_ConeInfo_List = Head_ConeInfo_List;
	Temp_Vector.SetLength(numOfVars);
	
	while (Temp_ConeInfo_List != NULL)
	{
		while (Temp_ConeInfo_List->ConeInfo_Pointer->Calculate_Integral_Point (Temp_Vector) == 1 )
		{
			solutions.push_front(Temp_Vector);	
		}
		
		Temp_ConeInfo_List = Temp_ConeInfo_List->Next;
	}
	*/
	
	//}
	//solutions.push_front( Temp_Vector );
	
  	// cout <<"indicator " << Indicator << endl;
  	/*if(Indicator == 0) 
	{
       		cout << "Error on IP." << endl;
    		//   printf("=======================\n"); 
    		

		//  printListCone(FindRationalFunction(cones, a, cost, numOfVars), numOfVars);
    		return Error;
    		//exit (1);
  	}
  	//   printListVector(matrix, numOfVars+1);*/

  	//else 
	// testing feasibility to fine integral point
	/*
  	OptVertex.SetLength(numOfVars);
  
	//cout << "Possible solutions: ";
	while(!solutions.empty())
	{
    		flag = 0;
    		possible = solutions.front();
	        cout << "Possible solution: " << possible << endl;	
    		solutions.pop_front(); //cout << possible << endl;
    		//cout << possible << " ";
		tmpmatrix = matrix;      
    		while(tmpmatrix)
		{ 
			RHS = 0; //cout <<"after while: " << possible << endl;
      			for(int i = 0; i < numOfVars; i++)
				RHS += (tmpmatrix -> first)[i + 1] * possible[i];
      			if((tmpmatrix -> first[0] + RHS) < 0)
			{
				flag = 1;
          			break;
			}
      			else ;
      			tmpmatrix = tmpmatrix -> rest;
    		}
    		if(flag == 0)
		{	
			OptVertex = possible; 
			break;
		}

		if(solutions.empty())
		{
			cout << "\n OH NO!!!!###########################################################";
  			cout << "\nCost:  " << cost << "\t\tVertex:  " << OptVertex << "\n";
			exit (0);
		}
	}
	*/
	//cout << endl;
	//cout << "\nCost:  " << cost << "\t\tVertex:  " << OptVertex << "\n";
	
	// Delete the ConeInfos that resulted in highest term;
	ConeInfo_List *Delete_ConeInfo_List;	
	Temp_ConeInfo_List = Head_ConeInfo_List;
	
	while (Temp_ConeInfo_List != NULL)
	{
		Delete_ConeInfo_List = Temp_ConeInfo_List;
		Temp_ConeInfo_List = Temp_ConeInfo_List->Next;
		
		delete Delete_ConeInfo_List->ConeInfo_Pointer;
		delete Delete_ConeInfo_List;
	}
	
	if ( Verify_IP_Flag == 1)
	{	
		ZZ	Max_IP;
		Max_IP = OptVertex * cost;
		ZZ	For_Limit, For_Counter;;

		For_Limit = power_ZZ(3,numOfVars);	
	
		for (int i = 0; i < numOfVars; i++)
		{
		 	possible[i] = -1;
		}
	
		for (For_Counter = 0; For_Counter < For_Limit; For_Counter++)
		{
		// Check current possible value agains max_ip
		
			flag = 0;
			tmpmatrix = matrix;      
		    	while(tmpmatrix)
			{ 
				RHS = 0; 
      				for(int l = 0; l < numOfVars; l++)
					RHS += (tmpmatrix -> first)[l + 1] * ( possible[l] + OptVertex[l]);
      				if((tmpmatrix -> first[0] + RHS) < 0)
				{
					flag = 1;
          				break;
				}
      				tmpmatrix = tmpmatrix -> rest;
    			}
			
			if (flag == 0)
				if ( ( (possible + OptVertex )*cost) > Max_IP)
					cout << "Oracle Wrong!" << endl;

		

		//Do ternary addition.
	
			possible[0] += 1;

			for (int j = 0; j < numOfVars - 1; j++)
			{
				if (possible[j] == 2)
				{
					possible[j] = -1;
					possible[j+1] += 1;
				}
			}	
	
		}	
	}
	/*
	for (int i = 0;i < 75; i++)
	{	
		for (int j = 0;j < 105; j++)
			for (int k = -25; k < 35; k++)
			{

			possible[0] = i;
			possible[1] = j;
			possible[2] = k;
			flag = 0;
			tmpmatrix = matrix;      
	    		while(tmpmatrix)
			{ 
				RHS = 0; //cout <<"after while: " << possible << endl;
      				for(int l = 0; l < numOfVars; l++)
					RHS += (tmpmatrix -> first)[l + 1] * possible[l];
      				if((tmpmatrix -> first[0] + RHS) < 0)
				{
					flag = 1;
          				break;
				}
      				else ;
      				tmpmatrix = tmpmatrix -> rest;
    			}
			
			if (flag == 0)
				if ( (possible*cost) > Max_IP)
						Max_IP = possible*cost;

			}
	}
	
	if ( (OptVertex * cost) != Max_IP)
		cout << "Oracle wrong!" << endl;
	*/
	
	return OptVertex;
	
	
}


listVector* GetHRepresentation(listVector* vertices, int numOfVars){

  int i, numOfHyperplane = 0, tmpdim, dummy;
  listVector* tmp, *basis, *endBasis;
  string tmpString;

  ofstream OUT;
  OUT.open("IH.ext");
  OUT << "V-representation" << endl;
  OUT << "begin " << endl;
  OUT << lengthListVector(vertices) << " " << numOfVars + 1 << " integer" << endl;
  tmp=vertices;
  while (tmp) {
    OUT << 1 << " ";
    for (i = 0; i < (numOfVars); i++) OUT << (tmp -> first)[i] << " ";
    OUT << endl;
    tmp=tmp->rest;
  }
  OUT << "end" << endl;
  OUT << "hull" << endl;
  OUT.close();

  system(LRS_PATH " IH.ext > IH.ine");

 ifstream in("IH.ine");
  if(!in){
    cerr << "Cannot open input file in IH.ine file." << endl;
    exit(1);
  }

  while (tmpString!="begin") getline(in,tmpString);
  while (tmpString!="end"){ getline(in,tmpString); numOfHyperplane++;}

  numOfHyperplane = numOfHyperplane - 2;
  // cout << numOfHyperplane << endl;
 ifstream in2("IH.ine");
  if(!in2){
    cerr << "Cannot open input file in IH.ine file." << endl;
    exit(2);
  }

  while (tmpString != "begin") getline(in2,tmpString);
  in2 >> tmpString >> tmpdim >> tmpString;
  vec_ZZ Hyperplane;

 basis = createListVector(createVector(numOfVars));
  endBasis = basis;

  for (i = 0; i < numOfHyperplane; i++) {
    Hyperplane = createVector(numOfVars); in2 >> dummy; 
    for (int j = 0; j < (numOfVars); j++){ in2 >> Hyperplane[j];
                                        Hyperplane[j] = - Hyperplane[j];} 
    if(!IsZero(Hyperplane)){
       endBasis->rest = createListVector(Hyperplane);
       endBasis = endBasis->rest;
    }
  }
  //system("rm IH.*");

//     printf("List of vectors:\n");
//   printf("================\n");
//   printListVector(basis->rest,numOfVars); 
   
  return(basis->rest);

}

listVector* GetVertices(listCone* cones, listVector* matrix, listVector* hyperplane, int numOfVars, int flag)
{
  	listVector* vertices, *endVertex, *tmpVertices;
  	vec_ZZ cost, vertex;
  	
	cost.SetLength(numOfVars);
  	
	if(flag == 0)
	{
    		if(hyperplane == 0)
		{
      			vertices = createListVector(createVector(numOfVars));
      			endVertex = vertices;
      			for(int j = 0; j < numOfVars; j++)
			{
				for(int i = 0; i < numOfVars; i++)
	  			{
					cost[i] = rand() % 100; 
					if(rand() % 2 == 0) 
						cost[i] = - cost[i];
				
				}
				vertex = SolveIP(cones, matrix, cost, numOfVars, 0); // cout << vertex << endl;
				endVertex -> rest = createListVector(vertex);
				endVertex = endVertex -> rest;
      			}
    		}
    
    		else
		{
      			tmpVertices = hyperplane;
      			vertices = createListVector(createVector(numOfVars));
      			endVertex = vertices;
      			//   printListVector(hyperplane, numOfVars); 
      			while(tmpVertices)
			{
				vertex = SolveIP(cones, matrix, tmpVertices -> first, numOfVars, 0);  
				endVertex -> rest = createListVector(vertex);
				endVertex = endVertex -> rest;   
				tmpVertices = tmpVertices -> rest;
      			}
    		}
  	}

  	else if(flag == 1)
	{
    		vertices = createListVector(createVector(numOfVars));
    		endVertex = vertices;
    		cout << "Enter a cost function." << endl;
    		vec_ZZ cost;
    		cost.SetLength(numOfVars);
    		for(int i = 0; i < numOfVars; i++)  
			cin >> cost[i];
    		vertex = SolveIP(cones, matrix, cost, numOfVars, 0);
    		endVertex -> rest = createListVector(vertex);
    		endVertex = endVertex -> rest; 
  	}

	//   printf("List of vectors:\n");
	//   printf("================\n");
	//   printListVector(vertices -> rest,numOfVars); 
  	return	vertices -> rest;
}

int CheckVertices(listVector* vertices, listVector* newVertices)
{
  	int flag = 0, len1 = 0, len2 = 0, counter = 0;
  	listVector * tmpvertices, *tmpnewVertices;
  	vec_ZZ vertex, newvertex;

  	tmpvertices = vertices;

  	len1 = lengthListVector(vertices);
  	len2 = lengthListVector(newVertices);
  	// cout << len1 << " " << len2 << endl;
  	for(int i = 0; i < len1; i++)
	{
    		vertex = tmpvertices -> first; 
    		tmpnewVertices = newVertices;
    		for(int j = 0; j < len2; j++)
		{
      			newvertex = tmpnewVertices -> first; 
      			if(vertex == newvertex) 
				counter++;
      			tmpnewVertices = tmpnewVertices -> rest;
    		}
    		tmpvertices = tmpvertices -> rest;
  	}
  	
	if(counter < len2)  
		flag = 1;// cout << counter << " " << len2 << endl;
  	
	return flag; 
}

listVector* Push_Vector(listVector* head, listVector* tail, int numOfVars){

  listVector* List, *endList;
  int len1, len2;
  vec_ZZ tmp;
  len1 = lengthListVector(head);
  len2 = lengthListVector(tail);
  List = createListVector(createVector(numOfVars));
  endList = List;
  vec_ZZ ArrayVec[len1];
  int flag = 0;

  for(int i = 0; i < (len1); i++) ArrayVec[i].SetLength(numOfVars); 
  for(int i = 0; i < len1; i++){
    tmp = head -> first;
    endList -> rest = createListVector(tmp);
    ArrayVec[i] = tmp;
    endList= endList -> rest;
    head = head -> rest;
  }
  for(int i = 0; i < len2; i++){
    flag = 0;
    tmp = tail -> first;
    for(int j = 0; j < len1; j++){
      if(tmp == ArrayVec[i]) flag = 1;
    }
    //    if(flag = 0){
      endList -> rest = createListVector(tmp);
      endList = endList -> rest;
      // }
    tail = tail -> rest;
  }

  return List -> rest;

}

ZZ Calculate_Polytope_Width (listCone *cones,listVector *matrix,int numOfVars)	
{
	//cout << "SolveIP Called. Cost = " << cost << endl;
	vec_ZZ OptVertex;
  	listCone * tmpcone;
  	listVector *tmpVector;
  	vec_ZZ	numerator, Temp_Vector, Max_Direction;
	ZZ 	S_Min, S_Max;

	S_Min = 0;
	S_Max = 0;
	
  
	Max_Direction.SetLength(numOfVars);
	
	for (int i = 0; i < numOfVars; i++)
		Max_Direction[i] = 0;
	
	for (int i = 0; i < numOfVars; i++)
	{	
	
		tmpcone = cones;
		while (tmpcone)
		{
			tmpVector = tmpcone -> rays;
	    		numerator = tmpcone ->latticePoints -> first;
	
			S_Max = numerator[i];
			S_Min = -numerator[i];

			while (tmpVector)
			{
				if (tmpVector->first[i] > 0)
					S_Max -= tmpVector->first[i];
				else if (tmpVector->first[i] < 0)
					S_Min += tmpVector->first[i];
				
				tmpVector = tmpVector->rest;
			}
			tmpcone = tmpcone->rest;

			if (S_Max < 0)
				S_Max *= -1;
			
			if (S_Min < 0)
				S_Min *= -1;
		
			if (S_Max > S_Min)
			{
				if (S_Max > Max_Direction[i])
					Max_Direction[i] = S_Max;
			}
			else
			{
				if (S_Min > Max_Direction[i])
					Max_Direction[i] = S_Min;
			}
		}
		
	}

	
	ZZ	Root,N;
	
	Root = 1;

	N = Max_Direction*Max_Direction;
	//cout << "N = " << N " ";
	for (int i = 0; i < 1000; i++)
	{
		Root = ((Root + N/Root)/2);
		if(Root == 0)
			Root = 1;
	}
		
	Root += 1;

	Polytope_Max_Width = Root;
	

	//cout << "Calculate_Polytope_Width: Max width is = " << Root << endl;
	return Root;	
}

listVector* IntegralHull(listCone* cones,listVector* matrix, int numOfVars)
{
	listVector* vertices, *hyperplanes, *newVertices;
	
	if (IntegralHull_Flag == 1)
	{
		cout << "Computing Integer Hull: " ;
		Calculate_Polytope_Width (cones, matrix, numOfVars);	
		int counter = 1, len = 0;
	
		vertices = GetVertices(cones, matrix, 0, numOfVars, 0);
	
		for(int i = 0; i < numOfVars; i++)
		{
		
    			vertices = Push_Vector(vertices,GetVertices(cones, matrix, 0, numOfVars, 0), numOfVars);
    		}

  		len = lengthListVector(vertices);
  
		int  Hull_Counter = 0;
   		while(counter != 0)
		{
		       	if ((Hull_Counter % 100) == 0)
				cout << Hull_Counter << " Done. " << endl;
			
    			hyperplanes = GetHRepresentation(vertices, numOfVars);
   	 		newVertices = GetVertices(cones, matrix, hyperplanes, numOfVars, 0);
    			counter = CheckVertices(vertices, newVertices);// cout << counter << endl;
    			vertices = Push_Vector(vertices, newVertices, numOfVars);
       			
			Hull_Counter++;
		}
  	}
	else if (IntegralHull_Flag == 0)
	{
  		// Read in from file cost.fun
		
		ifstream Cost_File("cost.fun");

		if (Cost_File.fail())
			exit (1);

		vec_ZZ Cost_Vector;

		Cost_Vector.SetLength (numOfVars);
		int	Int_Read;
		int	Solve_Count = 0;
		
		cout << "Reading in file." << endl;	
		while (!Cost_File.eof())
		{
			for (int j = 0; j < numOfVars; j++)
			{
				if (Cost_File.eof())
					break;
				Cost_File >> Int_Read;

				Cost_Vector[j] = Int_Read;		
			}
			//cout << "ConeInfo Object_Count: " << ConeInfo::Get_Object_Count () << endl;	
			SolveIP (cones, matrix, Cost_Vector, numOfVars, 0);

			Solve_Count++;

			if ((Solve_Count % 500) == 0)
				cout << "Solve_Count[" << Solve_Count << "]" << endl;	
		}
  	
	}

	return	vertices;
}


