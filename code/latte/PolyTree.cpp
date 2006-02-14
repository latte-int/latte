#include "PolyTree.h"
#include <stdlib.h>
#include <fstream.h>

int PolyTree_Node::Print ()
{
	char	Operator;
	
	if (Node_Type == POLYTREE_EXP)
	{
		cout << "(";
		Children[0]->Print ();
		cout << "^" << Number_of_Children << ")";
		return 1;
	}
	
	if (Node_Type == POLYTREE_ADD) // +
		Operator = '+';
	if (Node_Type == POLYTREE_MUL) // *
	{
		Operator = '*';
		//Check for zero
		
		/*for (int i = 0; i < Number_of_Children; i++)
		{
			
			if (Children[i]->Node_Type == POLYTREE_T_NODE);
			{
				if (Children[i]->Check_For_Zero() == 1)
					return 0;
			}
		}*/
	}

	if (Node_Type == POLYTREE_DIV)
		Operator = '/';
	
	cout << "(";

	for (unsigned int i = 0;i < Number_of_Children; i++)
	{
		if (Children[i]->Print () == 1)
		{
			if (i+1 != Number_of_Children)
				//if (Children[i+1]->Check_For_Zero () == 0 || Node_Type != POLYTREE_ADD)
					cout << " " << Operator << " ";
		}
	}
		
	cout << ")";

	return 1;	
}

int PolyTree_Node::Check_For_Zero ()
{
	if (Node_Type == POLYTREE_MUL)
	{
		for (unsigned int i = 0; i < Number_of_Children; i++)
		{
			if (Children[i]->Check_For_Zero () == 1)
				return 1;
		}

		return 0;
	}

	return 0;	
}

//Returns 1 if zero, else 0
int T_Node::Check_For_Zero ()
{
	if (Exponent == 0 && Coefficient == 0)
		return 1;
	else
		return 0;

}

/*int S_Node::Print ()
{
	cout << Coefficient << "s^" << Exponent;
	
	for (int i = 0;i < Number_of_Children; i++)
	{
		Children[i]->Print ();
	}

	return 1;
}*/

int T_Node::Print ()
{
	if (Node_Type == POLYTREE_T_NODE)
	{
		if (Exponent != 0 )
		{
			if (Coefficient != 1)
				cout << Coefficient << "*t^" << Exponent;
			else
				cout << "t^" << Exponent;
		}
		else 
			cout << Coefficient;
	}
	
	
	return 1;
}


void PolyTree_Node::Taylor_Expansion (Taylor_Parameters *Parameters)
{
	Taylor_Parameters	*Temp_Parameters;

	//cout << "Taylor_Expansion: Node_Type: " << Node_Type << endl;
	//Check if this NODE has already calculated its taylor expansion result
		
	if ( Taylor_Expansion_Result_Dirty == 0)
	{
		//cout << "clean!" << endl;
		for (int i = 0; i <= Parameters->Degree_of_Expansion; i++)
		{
			Parameters->Result[i] = Taylor_Expansion_Result[i];
		}

		return;
	}
	//cout << "dirty :(" << endl;
	
	if (Node_Type == POLYTREE_ADD)
	{
		//cout << "Expansion: ADD" << endl;
		
		Temp_Parameters = new Taylor_Parameters;
		Temp_Parameters->Degree_of_Expansion = Parameters->Degree_of_Expansion;
		Temp_Parameters->Ten_Power = Parameters->Ten_Power;
		
		Temp_Parameters->Result = new ZZ [Parameters->Degree_of_Expansion + 1];
		
		//Assume Result is DIRTY, so zero out all terms before continuing		
		for (int i = 0; i <= Parameters->Degree_of_Expansion;i++)
			Parameters->Result[i] = 0;

		
		for (unsigned int i = 0; i < Number_of_Children; i++)
		{
			Children[i]->Taylor_Expansion(Temp_Parameters);

			for (int j = 0; j <= Parameters->Degree_of_Expansion; j++)
			{
				Parameters->Result[j] += Temp_Parameters->Result[j];
			}
			
		}

		delete Temp_Parameters->Result;
		delete Temp_Parameters;

		

		/*for (int g = 0; g <= Parameters->Degree_of_Expansion; g++)
		  {
		    cout << Parameters->Result[g] << "t^" << g << " + " ;
		  }
		cout << endl;
		cout << "Expansion: ADD END" << endl;
		*/
		
	}
	else if (Node_Type == POLYTREE_MUL)
	{
		//cout << "Expansion: MUL" << endl;

		/*for (int i = 0;i < Number_of_Children; i++)
			if ( Children[i]->Check_For_Zero () == 1)
			{
				for (int k = 0; k <= Parameters->Degree_of_Expansion; k++)
				{
					Parameters->Result[0] = 0;
				}
					
				return;
			}
			*/
		
		Temp_Parameters = new Taylor_Parameters;
		Temp_Parameters->Degree_of_Expansion = Parameters->Degree_of_Expansion;
		Temp_Parameters->Ten_Power = Parameters->Ten_Power;
		
		Temp_Parameters->Result = new ZZ[Parameters->Degree_of_Expansion + 1];
		//cout << "done in MUL" << endl;
	
		// Copy first childs result into Parameters for shmushing to commence.	
		Children[0]->Taylor_Expansion( Parameters );
	 	//cout << "Done with call on first child" << endl;	
		
		// SMUSHING
		for (unsigned int i = 1; i < Number_of_Children; i++)
		{
			Children[i]->Taylor_Expansion( Temp_Parameters );

			for (int j = Parameters->Degree_of_Expansion; j >= 0; j--)
			{
				Parameters->Result[j] *= Temp_Parameters->Result[0];
				for (int k = 1; k <= j; k++)
				{	
					Parameters->Result[j] += Temp_Parameters->Result[k]*Parameters->Result[j-k];
				}
			}
			
		}
		//cout << "Done smushing" << endl;

		delete Temp_Parameters->Result;
		delete Temp_Parameters;

		
		/*for (int g = 0; g <= Parameters->Degree_of_Expansion; g++)
		  {
		    cout << Parameters->Result[g] << "t^" << g << " + " ;
		  }
		cout << endl;
		cout << "Expansion: MUL END " << endl;
		*/
		
	}
	else if (Node_Type == POLYTREE_DIV)
	{
		//cout << "Expansion: DIV" << endl;
		
		Taylor_Parameters	*Temp_Numerator_Parameters = new Taylor_Parameters;
		Taylor_Parameters	*Temp_Denominator_Parameters = new Taylor_Parameters;
		ZZ			B_Zero;
		
		//Get results numerators
		Temp_Numerator_Parameters->Degree_of_Expansion = Parameters->Degree_of_Expansion;
		Temp_Numerator_Parameters->Ten_Power = Parameters->Ten_Power;
		
		Temp_Numerator_Parameters->Result = new ZZ[Parameters->Degree_of_Expansion + 1];
		

		Children[0]->Taylor_Expansion ( Temp_Numerator_Parameters );

		Temp_Denominator_Parameters->Degree_of_Expansion = Parameters->Degree_of_Expansion;
		Temp_Denominator_Parameters->Ten_Power = Parameters->Ten_Power;
		Temp_Denominator_Parameters->Result = new ZZ[Parameters->Degree_of_Expansion + 1];

		Children[1]->Taylor_Expansion ( Temp_Denominator_Parameters );

		Parameters->Result[0] = Temp_Numerator_Parameters->Result[0]; //C_0 = a_0
 
		for (int k = 1; k <= Parameters->Degree_of_Expansion; k++ )
		{
			B_Zero = 1;
			
			Parameters->Result[k] = 0;
			
			for (int j = 1; j <= k; j++)
			{	
				Parameters->Result[k] -= Temp_Denominator_Parameters->Result[j]*Parameters->Result[k-j]*B_Zero;
				B_Zero *= Temp_Denominator_Parameters->Result[0]; // b_0
			}
			Parameters->Result[k] += B_Zero*Temp_Numerator_Parameters->Result[k];
			
		}

		B_Zero = Temp_Denominator_Parameters->Result[0]; //b_0
		
		for (int k = 0; k <= Parameters->Degree_of_Expansion; k++)
		{
			Parameters->Result[k] *= *(Parameters->Ten_Power);
			Parameters->Result[k] /= B_Zero;

			B_Zero *= Temp_Denominator_Parameters->Result[0];
			
		}	

		
		/*for (int g = 0; g <= Parameters->Degree_of_Expansion; g++)
		  {
		    cout << Parameters->Result[g] << "t^" << g << " + " ;
		  }
		cout << endl;*/
		//cout << "Expansion: DIV END" << endl;
		
		
		
		delete Temp_Numerator_Parameters;
		
	      
		delete Temp_Denominator_Parameters;
		
	}
	if (Node_Type == POLYTREE_EXP)
	{
		//cout << "Expansion: EXP" << endl;
		
	  
		unsigned int Exponent_Integer = Number_of_Children;
		unsigned int mask = 1 << (8 * (sizeof (unsigned int)) - 1);
		Taylor_Parameters *Temp_Parameters;
		
		int Bit_Position;

		

		
		for(Bit_Position = 8 * (sizeof (unsigned int)) - 1; (mask & Exponent_Integer) == 0; Bit_Position--)
			mask >>= 1;

		
		Temp_Parameters = new Taylor_Parameters;
		Temp_Parameters->Degree_of_Expansion = Parameters->Degree_of_Expansion;
		Temp_Parameters->Ten_Power = Parameters->Ten_Power;
		Temp_Parameters->Result = new ZZ[Parameters->Degree_of_Expansion + 1];
		       	

		Children[0]->Taylor_Expansion(Temp_Parameters);

		//cout << "EXP: Copying result of zero child to Parameters->Result...";
		for(int i = 0; i <= Parameters->Degree_of_Expansion; i++)
		  Parameters->Result[i] = Temp_Parameters->Result[i];
		//cout << "done." << endl;

		for(int m = 0; m < Bit_Position; m++)
		  {
		    	mask >>= 1;
		    
		    	for (int j = Parameters->Degree_of_Expansion; j >= 0; j--)
			{

			  	if(j > 0)
			    	{
					Parameters->Result[j] *= 2 * Parameters->Result[0];
					for (int k = 1; k < j; k++)
					{	
						Parameters->Result[j] += Parameters->Result[k]*Parameters->Result[j-k];
					}
			    	}
			  	else // j == 0
			    		Parameters->Result[j] *= Parameters->Result[j];

				}

		    		if(mask & Exponent_Integer)
		      		{

					for (int j = Parameters->Degree_of_Expansion; j >= 0; j--)
					{
						Parameters->Result[j] *= Temp_Parameters->Result[0];
						for (int k = 1; k <= j; k++)
						{	
							Parameters->Result[j] += Temp_Parameters->Result[k]*Parameters->Result[j-k];
						}
					}

		      		}
		     
		  	}
		

		delete Temp_Parameters->Result;
		delete Temp_Parameters;
		
		
		/*for (int g = 0; g <= Parameters->Degree_of_Expansion; g++)
		  {
		    cout << Parameters->Result[g] << "t^" << g << " + " ;
		  }
		cout << endl;
		cout << "Expansion: EXP END" << endl;
		*/
		
	}

	//Copy our results so that we do not have to calculate again for this node
	

	//cout << "Copying results for next time...";
	for (int i = 0;i <= Parameters->Degree_of_Expansion; i++)
	{		
		Taylor_Expansion_Result[i] = Parameters->Result[i];

	}
	//cout << "done." << endl;

	// Hooray its not dirty!
	Taylor_Expansion_Result_Dirty = 0;
}


void T_Node::Taylor_Expansion (Taylor_Parameters *Parameters)
{

	if (Node_Type == POLYTREE_T_NODE)
	{
		//cout << "Expansion: T_Node" << endl;

		for (int i = 0; i < Parameters->Degree_of_Expansion + 1; i++)
		{
			Parameters->Result[i] = 0;
		}

		//cout << "Expansion T_Node: results zerod" << endl;
		if (Exponent < 0 )
		{
			cerr << "Exponent of T_Node is negative.  Not supposed to happen!" << endl;
			exit (1);
		}
	
		int	Index_of_Result;
		
		conv(Index_of_Result, Exponent); 
				
		if (Exponent <= Parameters->Degree_of_Expansion)
			Parameters->Result[Index_of_Result] = Coefficient;
		
		//cout << "Result set" << endl;
		
		/*for (int g = 0; g <= Parameters->Degree_of_Expansion; g++)
		  {
		    cout << Parameters->Result[g] << "t^" << g << " + " ;
		  }
		cout << endl;
		cout << "Expansion: T_Node END" << endl;
		*/
		
	}
	

}

int PolyTree_Node::Print_Rational_Functions_to_File (ofstream &Output_File)
{
	char	Operator;
	

	if (!Output_File)
	{
		cerr << "Error opening output file in Print_Rational_Functions_to_File" << endl;
		exit (1);
	}
	if (Node_Type == POLYTREE_EXP)
	{
		Output_File << "(";
		Children[0]->Print_Rational_Functions_to_File (Output_File);
		Output_File << "^" << Number_of_Children << ")";
		return 1;
	}
	
	if (Node_Type == POLYTREE_ADD) // +
		Operator = '+';
	if (Node_Type == POLYTREE_MUL) // *
	{
		Operator = '*';
		//Check for zero
		
		/*for (int i = 0; i < Number_of_Children; i++)
		{
			
			if (Children[i]->Node_Type == POLYTREE_T_NODE);
			{
				if (Children[i]->Check_For_Zero() == 1)
					return 0;
			}
		}*/
	}

	if (Node_Type == POLYTREE_DIV)
		Operator = '/';
	
	
	
	Output_File << "(";

	for (unsigned int i = 0;i < Number_of_Children; i++)
	{
		if (Children[i]->Print_Rational_Functions_to_File (Output_File) == 1)
		{
			if (i+1 != Number_of_Children)
				//if (Children[i+1]->Check_For_Zero () == 0 || Node_Type != POLYTREE_ADD)
					Output_File << " " << Operator << " ";
		}
	}
		
	Output_File << ")";


	return 1;	
}

int T_Node::Print_Rational_Functions_to_File (ofstream &Output_File)
{
	if (!Output_File)
	{
		cerr << "Error opening output file in Print_Rational_Functions_to_File" << endl;
		exit (1);
	}
  
	if (Node_Type == POLYTREE_T_NODE)
	{
		if (Exponent != 0 )
		{
		  
			if (Coefficient != 1)
				Output_File << "(" <<  Coefficient << ")" << "*t^" << Exponent;
			else
				Output_File << "t^" << Exponent;
		}
		else 
			Output_File << "(" << Coefficient << ")";
	}
	
	return 1;
}

Node_Controller::Node_Controller (int Dimension, int Degree)
{
	//cout << "Node_Controller Constructor: Dimension = " << Dimension << " Degree = " << Degree << endl;
	Dimension_Plus_One = Dimension + 1;
	Degree_of_Expansion = Degree;

	
	PolyTree_Node_Head = new PolyTree_Node_List;

	PolyTree_Node_Head->Data = new PolyTree_Node;

	PolyTree_Node_Head->Data->Number_of_Children = Dimension_Plus_One;
	PolyTree_Node_Head->Data->Children = new PolyTree_Node*[Dimension_Plus_One];
	// Taylor_Expansion_Result_Dirty already 1 when allocated
	
	PolyTree_Node_Head->Data->Taylor_Expansion_Result = new ZZ [Degree + 1];
	
	PolyTree_Node_Head->Next = NULL;

	PolyTree_Node_Unused = PolyTree_Node_Head;

	
	T_Node_Head = 	new T_Node_List;

	T_Node_Head->Data = new T_Node;
	T_Node_Head->Next = NULL;

	T_Node_Unused = T_Node_Head;
	
}

PolyTree_Node *Node_Controller::Get_PolyTree_Node ()
{
	//cout << "Node_Controller: Get_PolyTree_Node called." << endl;
  	PolyTree_Node *Return_Value;

 	if (PolyTree_Node_Unused->Next != NULL)
    	{
      		Return_Value = PolyTree_Node_Unused->Data;

      		Return_Value->Taylor_Expansion_Result_Dirty = 1; //DIRTY!
      
      		PolyTree_Node_Unused = PolyTree_Node_Unused->Next;

      	}
  	else //PolyTree_Node_Unused->Next == NULL
    	{
      		Return_Value = PolyTree_Node_Unused->Data;
      
      		Return_Value->Taylor_Expansion_Result_Dirty = 1; //DIRTY!
      
      		PolyTree_Node_Unused->Next = new PolyTree_Node_List;
     		PolyTree_Node_Unused = PolyTree_Node_Unused->Next;

      		PolyTree_Node_Unused->Data = new PolyTree_Node;

      		PolyTree_Node_Unused->Data->Number_of_Children = Dimension_Plus_One;
      		PolyTree_Node_Unused->Data->Children = new PolyTree_Node*[Dimension_Plus_One];
      		// Taylor_Expansion_Result_Dirty already 1 when allocated
	
      		PolyTree_Node_Unused->Data->Taylor_Expansion_Result = new ZZ [Degree_of_Expansion + 1];
	
      		PolyTree_Node_Unused->Next = NULL;
     		
    }
  
  return Return_Value;
} 

T_Node *Node_Controller::Get_T_Node ()
{
	//cout << "Node_Controller: Get_T_Node called." << endl;
	T_Node *Return_Value;

  	if (T_Node_Unused->Next != NULL)
    	{
      		Return_Value = T_Node_Unused->Data;

      		T_Node_Unused = T_Node_Unused->Next;
        
    	}
  	else //T_Node_Unused->Next == NULL
    	{
      		Return_Value = T_Node_Unused->Data;
      	
      		T_Node_Unused->Next = new T_Node_List;
      		T_Node_Unused = T_Node_Unused->Next;

      		T_Node_Unused->Data = new T_Node;

	
      		T_Node_Unused->Next = NULL;
      
    	}
  
  	return Return_Value;
} 

void Node_Controller::Reset ()
{
	PolyTree_Node_Unused = PolyTree_Node_Head;
	
	T_Node_Unused = T_Node_Head;	

}

Node_Controller::~Node_Controller ()
{
	PolyTree_Node_List	*PolyTree_Node_Iterator, *Temp_PolyTree_Node_List;
	T_Node_List		*T_Node_Iterator, *Temp_T_Node_List;
	int	PolyTree_Node_Count = 0;
	int	T_Node_Count = 0;
	
	PolyTree_Node_Iterator = PolyTree_Node_Head;

	while (PolyTree_Node_Iterator != NULL)
	{
		delete [] PolyTree_Node_Iterator->Data->Taylor_Expansion_Result;
	
		delete [] PolyTree_Node_Iterator->Data->Children;
		
		delete PolyTree_Node_Iterator->Data;
		
		Temp_PolyTree_Node_List = PolyTree_Node_Iterator;

		PolyTree_Node_Iterator = PolyTree_Node_Iterator->Next;
		
		delete Temp_PolyTree_Node_List;
		PolyTree_Node_Count++;
	}

	T_Node_Iterator = T_Node_Head;

	while (T_Node_Iterator != NULL)
	{
		delete T_Node_Iterator->Data;

		Temp_T_Node_List = T_Node_Iterator;

		T_Node_Iterator = T_Node_Iterator->Next;

		delete Temp_T_Node_List;
		T_Node_Count++;
	}

	//cout << "Controller Destructor: Deleted " << PolyTree_Node_Count << " PolyTree_Node's." << endl;
	//cout << "Controller Destructor: Deleted " << T_Node_Count << " T_Node's." << endl;



}





