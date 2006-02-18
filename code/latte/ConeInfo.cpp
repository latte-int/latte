#include <stdlib.h>
#include <time.h>
#include "ConeInfo.h"
#include <iostream>

using namespace std;

Hash_Table_Entry ConeInfo::Hash_Table[HASH_TABLE_SIZE];
int	ConeInfo::Object_Count = 0;
int	*ConeInfo::Hash_Function_Coefficients;
int	ConeInfo::Time_Stamp = 1;
int	ConeInfo::Hash_Table_Initialized_Flag = 0;
Vector_Heap_Array_Node_Controller *Vector_Heap::Controller = NULL;

Vector_Heap::Vector_Heap (int Initial_Number_of_Generators)
{
	Root_Node = NULL;
	Node_Count = 0;
	Number_of_Generators = Initial_Number_of_Generators;

	if (Controller == NULL)
		Controller = new Vector_Heap_Array_Node_Controller(Initial_Number_of_Generators); 

	if (Initial_Number_of_Generators != Controller->Get_Current_Integer_Array_Size ())
		cout << "Vector_Heap Constructor: Problem!!!" << endl;
}

Vector_Heap::~Vector_Heap ()
{
	if (Root_Node != NULL)
		Delete_Sub_Tree (Root_Node);

}

void Vector_Heap::Delete_Sub_Tree(Heap_Node *Temp_Heap_Node)
{
	//delete [] Temp_Heap_Node->Sum_Vector;
	//delete Temp_Heap_Node->Total_Sum;

	// Use Controller to delete
	Controller->Recieve_ZZ (Temp_Heap_Node->Total_Sum);
	Controller->Recieve_Integer_Array (Temp_Heap_Node->Sum_Vector);
	
	if (Temp_Heap_Node->Left != NULL)
	{
		Delete_Sub_Tree (Temp_Heap_Node->Left);
		
	}

	if (Temp_Heap_Node->Right != NULL)
	{
		Delete_Sub_Tree (Temp_Heap_Node->Right);
	}
		
	delete Temp_Heap_Node;

}

void Vector_Heap::Add_Heap (int *Sum_Vector, ZZ *Total_Sum)
{
	if (Node_Count != 0)
	{
		Node_Count++;
		
		unsigned int	Temp_Node_Mask = 1 << 31;
	
		while ( !( Temp_Node_Mask & Node_Count)  )
		{
			Temp_Node_Mask =  Temp_Node_Mask >> 1;
		}
		
		Temp_Node_Mask =  Temp_Node_Mask >> 1;

		
		Heap_Node *Temp_Heap_Node = Root_Node;
		
		while ( Temp_Node_Mask > 1)
		{
			if ( (Temp_Node_Mask & Node_Count) == 0) //Left
				Temp_Heap_Node = Temp_Heap_Node->Left;
			else //right
				Temp_Heap_Node = Temp_Heap_Node->Right;
			
			Temp_Node_Mask =  Temp_Node_Mask >> 1;
		}	

		if ( (Temp_Node_Mask & Node_Count) == 0) //Add to left
		{
			//cout << "Add to left" << endl;
			Temp_Heap_Node->Left = new Heap_Node;

			Temp_Heap_Node->Left->Parent = Temp_Heap_Node;
			Temp_Heap_Node->Left->Right = Temp_Heap_Node->Left->Left = NULL;
		
			Temp_Heap_Node = Temp_Heap_Node->Left;
		}
		else //Add to Right
		{
			//cout << "Add to right" << endl;
			Temp_Heap_Node->Right = new Heap_Node;

			Temp_Heap_Node->Right->Parent = Temp_Heap_Node;
			Temp_Heap_Node->Right->Right = Temp_Heap_Node->Right->Left = NULL;
		
			Temp_Heap_Node = Temp_Heap_Node->Right;
		}
	
		Temp_Heap_Node->Sum_Vector = new int [ Number_of_Generators];
		// Use Controller
		//Temp_Heap_Node->Sum_Vector = Controller->Get_Integer_Array ();
	
		//cout << "Vector added to heap for cone: " << this << " [";
		for (int i = 0;i < Number_of_Generators; i++)
		{
			//cout << Sum_Vector[i] << " ";
			Temp_Heap_Node->Sum_Vector[i] = Sum_Vector[i];
		}
		//cout << "]" << endl;

		//Temp_Heap_Node->Total_Sum = new ZZ;
		// Use Controller
		Temp_Heap_Node->Total_Sum = Controller->Get_ZZ ();
		
		
		*(Temp_Heap_Node->Total_Sum) = *Total_Sum;

	
		//resort
		//cout << "Calling restore" << endl;
		Restore_Up (Temp_Heap_Node);
	}
	else // no nodes exist, create first. 
	{
		//cout << "No nodes.  Adding to root" << endl;
		Root_Node = new Heap_Node;

		Root_Node->Left = Root_Node->Right = NULL;

		//Root_Node->Sum_Vector = new int [Number_of_Generators];
		//Root_Node->Total_Sum = new ZZ;
		// Use Controller
		Root_Node->Sum_Vector = Controller->Get_Integer_Array ();			
		Root_Node->Total_Sum = Controller->Get_ZZ ();
		
		//cout << "Vector added to heap for cone: " << this << " [";
		for (int i = 0;i < Number_of_Generators; i++)
		{
			Root_Node->Sum_Vector[i] = Sum_Vector[i];
			//cout << Root_Node->Sum_Vector[i] << " ";
		}
		//cout << "]" << endl;

		*(Root_Node->Total_Sum) = *Total_Sum;

		Root_Node->Parent = NULL;
		Node_Count++;
			
	}
}

int Vector_Heap::Get_Top_Heap (int Sum_Vector [], ZZ *Total_Sum)
{
	if (Root_Node != NULL)
	{
		for (int i = 0; i < Number_of_Generators; i++)
		{
			Sum_Vector[i] = Root_Node->Sum_Vector[i];
		}
		
		*Total_Sum = *(Root_Node->Total_Sum);
		
		return 1; //sucessfull, that is root_node not null
	}
	else
		return 0; //failure, root_node null
	
}

int Vector_Heap::Pop_Top_Heap (int Sum_Vector [], ZZ *Total_Sum)
{
	//Print_Tree ();
	
	if (Root_Node != NULL)
	{
		
	if (Node_Count != 1)
	{
		//cout << "Pop_Top_Heap: Node_Count = " << Node_Count << ". [" << this << "]. ";
		//get the return values setup, then resort
		//cout << "[ ";
		for (int i = 0; i < Number_of_Generators; i++)
		{
			Sum_Vector[i] = Root_Node->Sum_Vector[i];
			//cout << Sum_Vector[i] << " ";
		}
		//cout << "]" << endl;
		//cout << *(Root_Node->Total_Sum);
		//cout << "For loop done. blarg" << endl;
		
		*Total_Sum = *(Root_Node->Total_Sum);

		//cout << "blarg" << endl;
		//swap values of the lowest and root.
		//find lowest node

		unsigned int	Temp_Node_Mask = 1 << 31;
		
		while ( !( Temp_Node_Mask & Node_Count)  )
			Temp_Node_Mask = Temp_Node_Mask >> 1;

		//cout << "Temp_Node_Mask = " << Temp_Node_Mask;
		Temp_Node_Mask >>= 1;

		Heap_Node *Temp_Heap_Node = Root_Node;
		
		while ( Temp_Node_Mask )
		{
			if ( (Temp_Node_Mask & Node_Count) == 0) //Left
			{	
				//cout << "Left. ";
				Temp_Heap_Node = Temp_Heap_Node->Left;
			}
			else //right
			{	
				//cout << "Right. ";
				Temp_Heap_Node = Temp_Heap_Node->Right;
			}
			Temp_Node_Mask >>= 1;
			
		}	
		

		for (int i = 0; i < Number_of_Generators; i++)
		{
			Root_Node->Sum_Vector[i] = Temp_Heap_Node->Sum_Vector[i];
		}
		
		*(Root_Node->Total_Sum) = *(Temp_Heap_Node->Total_Sum);
			
		//delete bottom node
	
		//delete Temp_Heap_Node->Total_Sum;
		//delete [] Temp_Heap_Node->Sum_Vector;
		// Use Controller
		Controller->Recieve_Integer_Array (Temp_Heap_Node->Sum_Vector);
		Controller->Recieve_ZZ (Temp_Heap_Node->Total_Sum);
		
		if (Temp_Heap_Node->Parent->Left == Temp_Heap_Node)
		{  //Left
			Temp_Heap_Node = Temp_Heap_Node->Parent;
			
			delete Temp_Heap_Node->Left;
			Temp_Heap_Node->Left = NULL;
			
		}
		else
		{
			Temp_Heap_Node = Temp_Heap_Node->Parent;
			delete Temp_Heap_Node->Right;
			Temp_Heap_Node->Right = NULL;
			
		}

		Node_Count--; //lower count

		//resort 

		Restore_Down (Root_Node);

		return 1;
	}
	else
	{	
		//cout << "Pop_Top_Heap: Node_Count = " << Node_Count << ". [" << this << "]. ";
		//get the return values setup, then resort
		//cout << "[ ";

		for (int i = 0; i < Number_of_Generators; i++)
		{
			Sum_Vector[i] = Root_Node->Sum_Vector[i];
			//cout << Sum_Vector[i] << " ";
		}
		//cout << "]" << endl;
		*Total_Sum = *(Root_Node->Total_Sum);
	
	
		//delete Root_Node->Total_Sum;
		//delete [] Root_Node->Sum_Vector;
		// Use Controller
		Controller->Recieve_Integer_Array (Root_Node->Sum_Vector);
		Controller->Recieve_ZZ (Root_Node->Total_Sum);

		
		delete Root_Node;
		Root_Node = NULL;

		
		Node_Count = 0;
		return 1;
	}
	}
	else
		return 0;


}

void Vector_Heap::Restore_Down (Heap_Node *Temp_Heap_Node)
{
	int Swap_Right = 1;
	int Swap_Left = 1;

	if (Temp_Heap_Node->Left == NULL)
		Swap_Left = 0;
	else if ( *(Temp_Heap_Node->Left->Total_Sum) <= *(Temp_Heap_Node->Total_Sum))
		Swap_Left = 0;
			
	if (Temp_Heap_Node->Right == NULL)
		Swap_Right = 0;	
	else if ( *(Temp_Heap_Node->Right->Total_Sum) <= *(Temp_Heap_Node->Total_Sum))
		Swap_Right	= 0;
	
	if (Swap_Left == 1 && Swap_Right == 1)
		if ( *(Temp_Heap_Node->Left->Total_Sum) > *(Temp_Heap_Node->Right->Total_Sum) ) //swap left and parent
			Swap_Right = 0;
		else
			Swap_Left = 0;
		
	if (Swap_Left == 1)
	{
		int	*Temp_Int;
		ZZ	*Temp_ZZ;
		
		Temp_Int = Temp_Heap_Node->Sum_Vector;
		Temp_ZZ = Temp_Heap_Node->Total_Sum;
		
		Temp_Heap_Node->Sum_Vector = Temp_Heap_Node->Left->Sum_Vector;
		Temp_Heap_Node->Total_Sum = Temp_Heap_Node->Left->Total_Sum;

		Temp_Heap_Node->Left->Sum_Vector = Temp_Int;
		Temp_Heap_Node->Left->Total_Sum = Temp_ZZ;
	
		Restore_Down (Temp_Heap_Node->Left);
	}
	if (Swap_Right == 1) //swap right and parent
	{
			
		int	*Temp_Int;
		ZZ	*Temp_ZZ;

		Temp_Int = Temp_Heap_Node->Sum_Vector;
		Temp_ZZ = Temp_Heap_Node->Total_Sum;

		Temp_Heap_Node->Sum_Vector = Temp_Heap_Node->Right->Sum_Vector;
		Temp_Heap_Node->Total_Sum = Temp_Heap_Node->Right->Total_Sum;

		Temp_Heap_Node->Right->Sum_Vector = Temp_Int;
		Temp_Heap_Node->Right->Total_Sum = Temp_ZZ;

		Restore_Down (Temp_Heap_Node->Right);
		
	}

}

int Vector_Heap::Check_Top_Heap (ZZ *Total_Sum)
{
	if ( Root_Node != NULL)
	{
		if ( *(Root_Node->Total_Sum) == *Total_Sum)
			return 1;
		else	
			return 0;
	}
	else
		return 0;
	
}

void Vector_Heap::Print_Tree ()
{
	if ( Root_Node != NULL)
	{
		Print_Sub_Tree (Root_Node);
		cout << endl;
	}
	else
		cout << "NULL" << endl;
}

void Vector_Heap::Print_Sub_Tree (Heap_Node *Temp_Heap_Node)
{
	/*cout << "Vector [";

	for (int i = 0; i < Number_of_Generators; i++)
	{
		cout << Temp_Heap_Node->Sum_Vector[i] << " ";
	}	
	cout << "] ";
	*/
	cout << "Sum: " << *(Temp_Heap_Node->Total_Sum) << "\t";

	if ( Temp_Heap_Node->Left != NULL)
	{
		//cout << "Left. ";	
		Print_Sub_Tree ( Temp_Heap_Node->Left);
	}
	if ( Temp_Heap_Node->Right != NULL)
	{
		//cout << "Right. ";
		Print_Sub_Tree ( Temp_Heap_Node->Right);
	}

}

void Vector_Heap::Restore_Up (Heap_Node *Temp_Heap_Node)
{
	if (Temp_Heap_Node->Parent != NULL)
	{
		if ( *(Temp_Heap_Node->Total_Sum) > *(Temp_Heap_Node->Parent->Total_Sum) )
		{
			//cout << "Restore_Up: swaping" << endl;			
			int	*Temp_Int;
			ZZ	*Temp_ZZ;

			Temp_Int = Temp_Heap_Node->Sum_Vector;
			Temp_ZZ = Temp_Heap_Node->Total_Sum;

			Temp_Heap_Node->Sum_Vector = Temp_Heap_Node->Parent->Sum_Vector;
			Temp_Heap_Node->Total_Sum = Temp_Heap_Node->Parent->Total_Sum;

			Temp_Heap_Node->Parent->Sum_Vector = Temp_Int;
			Temp_Heap_Node->Parent->Total_Sum = Temp_ZZ;

			Restore_Up(Temp_Heap_Node->Parent);
		}
		//else
			//cout << "No swap" << endl;
	}
	//else
		//cout << "Parent Null" << endl;
}

int Vector_Heap::Get_Node_Count ()
{
	return Node_Count;
}

ConeInfo::ConeInfo (vector *cost, listCone *listCone_pointer, int numOfVars)
{
	int NumGensPerCone = lengthListVector(listCone_pointer->rays);
	int *zero_vector = new int[NumGensPerCone];
	vector Our_Cost;
	
	Number_of_Variables = numOfVars;	
	
	Our_Cost = *cost;
	
	S_Values_Zero_Flag = 0;	

	S_Values = new ZZ[NumGensPerCone];
	signs = new int[NumGensPerCone];
	Numer_Exp = new ZZ;
	
	
	listCone_ptr = listCone_pointer;


		Coefficient = listCone_pointer->coefficient;
		
		/*if (S_Values_Zero_Flag == 1)
		{
			cout << "S_Value zero. Calculating Pertubation. " << this << endl;
			Our_Cost = Calculate_Pertubation (cost, 3);
			S_Values_Zero_Flag = 0;
		}*/
	
		*Numer_Exp = (listCone_pointer->latticePoints->first) * (Our_Cost);

		listVector *current = listCone_pointer->rays;
	
		//cout << "ConeInfo Constructor: Generators { ";	
		for(int i = 0; i < NumGensPerCone; i++)
		{
			//cout << current->first << " ";
			S_Values[i] = (Our_Cost) * (current->first);
			current = current->rest;

			if(S_Values[i] ==  0)
			{
				//cout << "ConeInfo Constructor: S_Value zero. " << this << endl;
				S_Values_Zero_Flag = 1;
			}
		
			if(S_Values[i] > 0)
			{
				*Numer_Exp -=  S_Values[i];
			
				Coefficient *= -1;
			
				S_Values[i] *= -1;
				// sign is the sign of the orginal dot product
				signs[i] = 1;
					
			}
			else
				signs[i] = -1;
					

			zero_vector[i] = 0;
		}
	//cout << " } " << endl;

	//cout << "Done computing S_Values. " << this << endl;

	
	//cout << "Coefficient after dot products " << Coefficient << endl;
	Number_of_Generators = lengthListVector(listCone_pointer->rays);	

	Order = new int[Number_of_Generators];
	

	/*cout << "ConeInfo Constructor: S_Values before sort: [ ";
	for (int i = 0; i < Number_of_Generators; i++)
	{
		cout << S_Values[i] << " ";
	}
	cout << "]" << endl;
	*/
	Sort_S_Values ();
	
	/*for (int i = 1; i < Number_of_Generators; i++)
	{
		if (S_Values[i] == S_Values[i-1])
			cout << "ConeInfo Constructor: Repeated S_Values terms." << endl;
	}*/
	zero_vector[0] = 1;
	
	/*ZZ Compare_Value;
	int All_Same = 1;
	
	Compare_Value = S_Values[0];
	for (int i = 1; i < Number_of_Generators; i++)
	{
		if (S_Values[i] != Compare_Value)
			All_Same = 0;
	}
	if (All_Same == 1)
		cout << "ConeInfo Constructor: S_Values all the same" << endl;
	*/
	Heap = new Vector_Heap(numOfVars);

	ZZ	Temp_Total_Sum;

	//Temp_Total_Sum = *Numer_Exp + S_Values[0];
	Temp_Total_Sum = S_Values[0];

	//cout << "ConeInfo Constructor: Adding 1 0 0" << endl;	
	Heap->Add_Heap(zero_vector, &Temp_Total_Sum);

	Integer_Vector_List	*New_Integer_Vector_List = new Integer_Vector_List;

	zero_vector[0] = 0;
	
	New_Integer_Vector_List->Next = NULL;
	New_Integer_Vector_List->Integer_Array = zero_vector;

	Vector_List_Head = New_Integer_Vector_List;

	Expansion_Highest_Term = new ZZ;
	*Expansion_Highest_Term = 0;

	Current_Highest_Term = new ZZ;
	*Current_Highest_Term = *Numer_Exp;



	if (Object_Count == 0) //create random array for hash table, only once.
	{
		if (Hash_Table_Initialized_Flag == 0)
		{
			//cout << "Creating hash coefficients." << endl;
			Hash_Function_Coefficients = new int[Number_of_Generators];

			srand(clock());
		
			for (int i = 0; i < Number_of_Generators; i ++)
			{
				Hash_Function_Coefficients[i] = rand();

			}
			Hash_Table_Initialized_Flag = 1;
		
	
			for (int i = 0; i < HASH_TABLE_SIZE; i++)
			{
				Hash_Table[i].Time_Stamp = 0;
			}
		}
	}

	Object_Count++;
}

ConeInfo::~ConeInfo ()
{
	Object_Count--;
	delete Heap;
	delete [] S_Values;
	delete [] signs;
	delete Numer_Exp;
	delete [] Order;
	delete Expansion_Highest_Term;
	delete Current_Highest_Term;
	
	/*if (Object_Count == 0)
	{
		//cout << "ConeInfo Destructor: Deleting hash functions coefficients." << endl;
		delete [] Hash_Function_Coefficients;

	}*/
	
	Integer_Vector_List *Temp_Integer_Vector_List, *Delete_Integer_Vector_List;

	Temp_Integer_Vector_List = Vector_List_Head;

	while (Temp_Integer_Vector_List != NULL)
	{
		Delete_Integer_Vector_List = Temp_Integer_Vector_List;
		Temp_Integer_Vector_List = Temp_Integer_Vector_List->Next;

		delete [] Delete_Integer_Vector_List->Integer_Array;
		delete Delete_Integer_Vector_List;
	}
}

int ConeInfo::Get_Object_Count ()
{
	return Object_Count;

}

int ConeInfo::Get_Time_Stamp ()
{
	return Time_Stamp;
}
void ConeInfo::Sort_S_Values ()
{
	ZZ	Temp_ZZ;
	int	Temp_int;
	
	for (int i = 0; i < Number_of_Generators; i++)
		Order[i] = i;

	// Bubble Sort!
	
	for(int i = 0; i < Number_of_Generators - 1; i++)
		for (int j = 0; j < Number_of_Generators - i - 1; j++)
	      	{
			if (S_Values[j+1] > S_Values [j])
			{
				Temp_ZZ = S_Values[j+1];
				S_Values[j+1] = S_Values[j];
				S_Values[j] = Temp_ZZ;

				Temp_int = Order[j+1];
				Order[j+1] = Order[j];
				Order[j] = Temp_int;
			}
		}	       

}


void ConeInfo::Calculate_Next_Highest_Term ()
{
	if (S_Values_Zero_Flag == 1)
	{
		cout << "Trying to dig on a S_Value zero cone! Exiting. " << this << endl;
		exit (1);
	}
	
	//cout << "CalcNextTerm: Coefficient = " << Coefficient << endl;
	/*cout << endl << "Svalues: ";
	for (int i = 0; i < Number_of_Generators; i++)
	{
		cout << S_Values[i] << ", ";

	}*/
	int First_Non_Zero;

	
	Integer_Vector_List	*Temp_Integer_Vector_List, *Delete_Node;
	
	Temp_Integer_Vector_List = Vector_List_Head;

	// Delete the old highest terms.
	while (Temp_Integer_Vector_List != NULL)
	{
		delete [] Temp_Integer_Vector_List->Integer_Array;
		
		Delete_Node = Temp_Integer_Vector_List;
		Temp_Integer_Vector_List = Temp_Integer_Vector_List->Next;
		
		delete Delete_Node;
	}

	//cout << "Calculation_Next: Data deleted." << endl;	
	
	// Retain the sign of the coefficient.
	if (Coefficient < 0)
		Coefficient = -1;
	else
		Coefficient = 1;
	
	ZZ	*New_Total_Sum;
	int	*New_Integer_Array = new int[Number_of_Generators];

	// Pop off the highest term off the top of the heap.
	if ( Heap->Pop_Top_Heap (New_Integer_Array, Expansion_Highest_Term) == 0)
		cout << "Error with Pop_Top in ConeInfo::Calculate_Next_Highest_Term " << endl;

	// Create new node.
	Temp_Integer_Vector_List = new Integer_Vector_List;
	// Store the array popped off the heap
	Temp_Integer_Vector_List->Integer_Array = New_Integer_Array;
	Temp_Integer_Vector_List->Next = NULL;

	Vector_List_Head = Temp_Integer_Vector_List;

	*Current_Highest_Term = *Numer_Exp + *Expansion_Highest_Term;
	
	if (*Expansion_Highest_Term > 0)
		//cout << "Expansion_Highest_Term Negative" << endl;
	//else
		cout << "Expansion_Highest_Term Positive." << endl;
	
	
	int	Hash_Value;
	
	// Increment the Time stamp for hashing of heap items.
	Time_Stamp++;

	if (Time_Stamp == 0)
		cout << "Calculate_Next_Highest_Term: Error, Time_Stamp == 0.  Roll Over!" << endl;	
	//debug
	//cout << Time_Stamp << endl;
	
	int *New_Integer_Array2;
	
	New_Total_Sum = new ZZ;

	//New_Integer_Array = new int[Number_of_Generators];
	
	//Allocate a New_Integer_Array2
	New_Integer_Array2 = new int[Number_of_Generators];

	for (int i = 0;i < Number_of_Generators; i++)
	{
		New_Integer_Array2[i] = Temp_Integer_Vector_List->Integer_Array[i];	
		
		//New_Integer_Array is copy used for hash table
		New_Integer_Array[i] = Temp_Integer_Vector_List->Integer_Array[i];	
	}
		
	New_Integer_Array2[0] += 1;

	*New_Total_Sum = *Expansion_Highest_Term + S_Values[0];

	// Add_Heap copies the data sent to it.  It does not retain the pointers.	
	Heap->Add_Heap(New_Integer_Array2,New_Total_Sum);

	*New_Total_Sum = *New_Total_Sum - S_Values[0];
	New_Integer_Array2[0] -= 1;

	First_Non_Zero = -1;

	for (int i = 0; i < Number_of_Generators; i++)
	{	
		if (First_Non_Zero == -1)
			if (New_Integer_Array2[i] != 0)
			{
				First_Non_Zero = i;
				break;
			}
	}
	
	//debug
	if(First_Non_Zero == Number_of_Generators)
		cout << "Calculate_Next_Highest_Term: First_Non_Zero == Number_of_Generators";
	
	if (First_Non_Zero != Number_of_Generators - 1 )
	{
		New_Integer_Array2[First_Non_Zero] -= 1;
		New_Integer_Array2[First_Non_Zero +1] += 1;
		
		*New_Total_Sum = *New_Total_Sum - S_Values[First_Non_Zero];
		*New_Total_Sum = *New_Total_Sum + S_Values[First_Non_Zero + 1];
		
		Heap->Add_Heap(New_Integer_Array2,New_Total_Sum);
	}
	//done with New_Integer_Array2, delete it, heap has its own copy.
	delete [] New_Integer_Array2;

	Hash_Value = Hash_Integer_Vector (New_Integer_Array);
	Hash_Table[Hash_Value].Time_Stamp = Time_Stamp;
	//Hash_Table stores the actual address, so dont delete New_Integer_Array;
	Hash_Table[Hash_Value].Integer_Array = New_Integer_Array;

	int Equal_Indicator;
	

	//Heap->Print_Tree ();
	while ( Heap->Check_Top_Heap(Expansion_Highest_Term) == 1)
	{
		New_Integer_Array = new int[Number_of_Generators];
		New_Total_Sum  = new ZZ;
		
		Heap->Pop_Top_Heap (New_Integer_Array, New_Total_Sum);
		Temp_Integer_Vector_List = Vector_List_Head;

		Equal_Indicator = 1;  //  1 is true   0 is false
		
	/*	while ( Temp_Integer_Vector_List != NULL)
		{
			Equal_Indicator = 1;  //  1 is true   0 is false
			
			for (int i = 0; i < Number_of_Generators; i++)
			{	
				if (Temp_Integer_Vector_List->Integer_Array[i] != New_Integer_Array[i])
				{
					Equal_Indicator = 0;
					break;
				}	
			}
			
			if (Equal_Indicator == 1)
				break;
			
			Temp_Integer_Vector_List = Temp_Integer_Vector_List->Next;
		}
	*/
		
		Hash_Value = Hash_Integer_Vector (New_Integer_Array);
		//cout << "Integer Vector hashed." << endl;	
			
		while(1)
		{
			// Linear probe the hash table for the value	
			if (Time_Stamp == Hash_Table[Hash_Value].Time_Stamp)
			{
				for (int j = 0; j < Number_of_Generators; j++)
				{	
					if (New_Integer_Array[j] != Hash_Table[Hash_Value].Integer_Array[j])
					{
						//cout << "Vectors not equal." << endl;
						Equal_Indicator = 0;
						break;
					}		
				}
				if (Equal_Indicator == 0)
				{
					Equal_Indicator = 1;
					Hash_Value = (Hash_Value + 1) % HASH_TABLE_SIZE;
				}
				else
					break;
			}
			else
			{
				Hash_Table[Hash_Value].Time_Stamp = Time_Stamp;
				Hash_Table[Hash_Value].Integer_Array = New_Integer_Array;
				Equal_Indicator = 0;
				break;
			}	
		}
		
		if (Equal_Indicator == 1)  // Vector already exists 
		{
			delete [] New_Integer_Array;
			delete New_Total_Sum;
		}
		else   // Vector does not exist.  Save it and store new items into heap.
		{
			Temp_Integer_Vector_List = new Integer_Vector_List;

			Temp_Integer_Vector_List->Integer_Array = New_Integer_Array;
			Temp_Integer_Vector_List->Next = Vector_List_Head;

			Vector_List_Head = Temp_Integer_Vector_List;
				

			if (Coefficient > 0)
				Coefficient++;
			else
				Coefficient--;


			// Allocate New_Integer_Array2	
			New_Integer_Array2 = new int[Number_of_Generators];
			
			for (int i = 0;i < Number_of_Generators; i++)
			{
				New_Integer_Array2[i] = Temp_Integer_Vector_List->Integer_Array[i];	
			}
		
			New_Integer_Array2[0] += 1;

			*New_Total_Sum = *Expansion_Highest_Term + S_Values[0];
		
			Heap->Add_Heap(New_Integer_Array2,New_Total_Sum);

			*New_Total_Sum = *New_Total_Sum - S_Values[0];
			New_Integer_Array2[0] -= 1;

			First_Non_Zero = -1;

			for (int i = 0; i < Number_of_Generators; i++)
			{	
				if (First_Non_Zero == -1)
					if (New_Integer_Array2[i] != 0)
					{
						First_Non_Zero = i;
						break;
					}
			}
			
			if (First_Non_Zero != Number_of_Generators - 1 )
			{
				New_Integer_Array2[First_Non_Zero] -= 1;
				New_Integer_Array2[First_Non_Zero +1] += 1;
		
				*New_Total_Sum = *New_Total_Sum - S_Values[First_Non_Zero];
				*New_Total_Sum = *New_Total_Sum + S_Values[First_Non_Zero + 1];
		
				Heap->Add_Heap(New_Integer_Array2,New_Total_Sum);
			}

			//Done with New_Integer_Array2, delete it and New_Total_Sum
			delete [] New_Integer_Array2;
			delete New_Total_Sum;
		}	
	
	}
	//Heap->Print_Tree ();
	//cout << "%%";	
	//Insert new vectors into the heap

	/*Temp_Integer_Vector_List = Vector_List_Head;

	New_Integer_Array = new int[Number_of_Generators];
	New_Total_Sum	= new ZZ;
	
	while (Temp_Integer_Vector_List != NULL)
	{
		for (int i = 0;i < Number_of_Generators; i++)
		{
			New_Integer_Array[i] = Temp_Integer_Vector_List->Integer_Array[i];	
		}
		
		New_Integer_Array[0] += 1;

		*New_Total_Sum = *Expansion_Highest_Term + S_Values[0];
		
		Heap->Add_Heap(New_Integer_Array,New_Total_Sum);

		*New_Total_Sum = *Expansion_Highest_Term - S_Values[0];
		New_Integer_Array[0] -= 1;

		int	First_Non_Zero = -1;

		for (int i = 0; i < Number_of_Generators; i++)
		{	
			if (First_Non_Zero == -1)
				if (New_Integer_Array[i] != 0)
				{
					First_Non_Zero = i;
					break;
				}
		}
			
		if (First_Non_Zero != Number_of_Generators - 1 )
		{
			New_Integer_Array[First_Non_Zero] -= 1;
			New_Integer_Array[First_Non_Zero +1] += 1;
		
			*New_Total_Sum = *New_Total_Sum - S_Values[First_Non_Zero];
			*New_Total_Sum = *New_Total_Sum + S_Values[First_Non_Zero + 1];
		
			Heap->Add_Heap(New_Integer_Array,New_Total_Sum);
		}
		
		Temp_Integer_Vector_List = Temp_Integer_Vector_List->Next;
	}	
	*/
	
	//delete [] New_Integer_Array;
	//delete New_Total_Sum;
	
	//cout << "CalcNextTerm Done: Coefficient: " << Coefficient << endl;
}

ConeInfo_Heap::ConeInfo_Heap ()
{
	Root_Node = NULL;
	Node_Count = 0;
}

ConeInfo_Heap::~ConeInfo_Heap ()
{
	if (Root_Node != NULL)
	{
		Delete_Sub_Tree (Root_Node);
	}
}

void ConeInfo_Heap::Clear_Tree ()
{
	if (Root_Node != NULL)
	{
		Delete_Sub_Tree (Root_Node);
		Root_Node = NULL;
		Node_Count = 0;
	}

}

void ConeInfo_Heap::Delete_Sub_Tree(ConeInfo_Heap_Node *Temp_ConeInfo_Heap_Node)
{
	//Actually deletes ConeInfo
	delete Temp_ConeInfo_Heap_Node->ConeInfo_Pointer;

	if (Temp_ConeInfo_Heap_Node->Left != NULL)
	{
		Delete_Sub_Tree (Temp_ConeInfo_Heap_Node->Left);
		
	}

	if (Temp_ConeInfo_Heap_Node->Right != NULL)
	{
		Delete_Sub_Tree (Temp_ConeInfo_Heap_Node->Right);
	}
		
	delete Temp_ConeInfo_Heap_Node;

}


void ConeInfo_Heap::Add_Heap (ConeInfo *Temp_ConeInfo)
{
	if (Node_Count != 0)
	{
		Node_Count++;
		
		unsigned int	Temp_Node_Mask = 1 << 31;
	
		while ( !( Temp_Node_Mask & Node_Count)  )
		{
			Temp_Node_Mask =  Temp_Node_Mask >> 1;
		}
		
		Temp_Node_Mask =  Temp_Node_Mask >> 1;

		
		ConeInfo_Heap_Node *ConeInfo_Temp_Heap_Node = Root_Node;
		
		while ( Temp_Node_Mask > 1)
		{
			if ( (Temp_Node_Mask & Node_Count) == 0) //Left
				ConeInfo_Temp_Heap_Node = ConeInfo_Temp_Heap_Node->Left;
			else //right
				ConeInfo_Temp_Heap_Node = ConeInfo_Temp_Heap_Node->Right;
			
			Temp_Node_Mask =  Temp_Node_Mask >> 1;
		}	

		if ( (Temp_Node_Mask & Node_Count) == 0) //Add to left
		{
			//cout << "Add to left" << endl;
			ConeInfo_Temp_Heap_Node->Left = new ConeInfo_Heap_Node;

			ConeInfo_Temp_Heap_Node->Left->Parent = ConeInfo_Temp_Heap_Node;
			ConeInfo_Temp_Heap_Node->Left->Right = ConeInfo_Temp_Heap_Node->Left->Left = NULL;
		
			ConeInfo_Temp_Heap_Node = ConeInfo_Temp_Heap_Node->Left;
		}
		else //Add to Right
		{
			//cout << "Add to right" << endl;
			ConeInfo_Temp_Heap_Node->Right = new ConeInfo_Heap_Node;

			ConeInfo_Temp_Heap_Node->Right->Parent = ConeInfo_Temp_Heap_Node;
			ConeInfo_Temp_Heap_Node->Right->Right = ConeInfo_Temp_Heap_Node->Right->Left = NULL;
		
			ConeInfo_Temp_Heap_Node = ConeInfo_Temp_Heap_Node->Right;
		}
	
		ConeInfo_Temp_Heap_Node->ConeInfo_Pointer = Temp_ConeInfo; 
	
		//resort
		//cout << "Calling restore" << endl;
		Restore_Up (ConeInfo_Temp_Heap_Node);
	}
	else // no nodes exist, create first. 
	{
		//cout << "No nodes.  Adding to root" << endl;
		Root_Node = new ConeInfo_Heap_Node;

		Root_Node->Left = Root_Node->Right = NULL;

		Root_Node->ConeInfo_Pointer = Temp_ConeInfo;
		
		Root_Node->Parent = NULL;
		Node_Count++;
			
	}
}



ConeInfo *ConeInfo_Heap::Get_Top_Heap ()
{
	if (Root_Node != NULL)
	{
		
		return Root_Node->ConeInfo_Pointer; //sucessfull, that is root_node not null
	}
	else
		return NULL; //failure, root_node null
	
}


ConeInfo *ConeInfo_Heap::Pop_Top_Heap ()
{
	//cout << "Pop_Top_Heap: Called." << endl;
	ConeInfo	*Return_Value = Root_Node->ConeInfo_Pointer;
	
	if (Root_Node != NULL)
	{
		
	if (Node_Count != 1)
	{
		//cout << "Pop_Top_Heap: Traversing tree. " << Node_Count << endl;
		//get the return values setup, then resort
		
		//swap values of the lowest and root.
		//find lowest node

		unsigned int	Temp_Node_Mask = 1 << 31;
		
		while ( !( Temp_Node_Mask & Node_Count)  )
			Temp_Node_Mask >>= 1;

		Temp_Node_Mask >>= 1;

		ConeInfo_Heap_Node *Temp_ConeInfo_Heap_Node = Root_Node;
	
		//cout << "POp_Top_Heap: Traversing tree..." << endl;	
		while ( Temp_Node_Mask )
		{
			if ( (Temp_Node_Mask & Node_Count) == 0) //Left
				Temp_ConeInfo_Heap_Node = Temp_ConeInfo_Heap_Node->Left;
			else //right
				Temp_ConeInfo_Heap_Node = Temp_ConeInfo_Heap_Node->Right;

			Temp_Node_Mask >>= 1;
			
			//if (Temp_ConeInfo_Heap_Node == NULL)
				//cout << "Pop_Top_Heap: Null found on right or left. " << Temp_Node_Mask << endl;
			
		}	
	
		//cout << "Pop_TOp_Heap: Traversing ended." << endl;	

		
		Root_Node->ConeInfo_Pointer = Temp_ConeInfo_Heap_Node->ConeInfo_Pointer;
			
		//delete bottom node
	

		//cout << "Pop_Top_Heap: Deleting bottom Node" << endl;
		if (Temp_ConeInfo_Heap_Node->Parent->Left == Temp_ConeInfo_Heap_Node)
		{  //Left
			Temp_ConeInfo_Heap_Node = Temp_ConeInfo_Heap_Node->Parent;
			
			delete Temp_ConeInfo_Heap_Node->Left;
			Temp_ConeInfo_Heap_Node->Left = NULL;
			
		}
		else
		{
			Temp_ConeInfo_Heap_Node = Temp_ConeInfo_Heap_Node->Parent;
			delete Temp_ConeInfo_Heap_Node->Right;
			Temp_ConeInfo_Heap_Node->Right = NULL;
			
		}

		Node_Count--; //lower count

		//resort 

		Restore_Down (Root_Node);

		return Return_Value;
	}
	else
	{	
		Return_Value = Root_Node->ConeInfo_Pointer;
		

		delete Root_Node;
		Root_Node = NULL;
		
		Node_Count = 0;
		return Return_Value;
	}
	}
	else
		return NULL;


}

void ConeInfo_Heap::Restore_Down (ConeInfo_Heap_Node *Temp_ConeInfo_Heap_Node)
{
	//cout << "Restore_Down: Called." << endl;
	int Swap_Right = 1;
	int Swap_Left = 1;

	if (Temp_ConeInfo_Heap_Node->Left == NULL)
		Swap_Left = 0;
	else if ( *(Temp_ConeInfo_Heap_Node->Left->ConeInfo_Pointer->Get_Current_Highest_Term () ) <= *(Temp_ConeInfo_Heap_Node->ConeInfo_Pointer->Get_Current_Highest_Term () ))
		Swap_Left = 0;
			
	if (Temp_ConeInfo_Heap_Node->Right == NULL)
		Swap_Right = 0;	
	else if ( *(Temp_ConeInfo_Heap_Node->Right->ConeInfo_Pointer->Get_Current_Highest_Term () ) <= *(Temp_ConeInfo_Heap_Node->ConeInfo_Pointer->Get_Current_Highest_Term () ))
		Swap_Right	= 0;
	
	if (Swap_Left == 1 && Swap_Right == 1)
		if ( *(Temp_ConeInfo_Heap_Node->Left->ConeInfo_Pointer->Get_Current_Highest_Term () ) > *(Temp_ConeInfo_Heap_Node->Right->ConeInfo_Pointer->Get_Current_Highest_Term () ) ) //swap left and parent
			Swap_Right = 0;
		else
			Swap_Left = 0;
		
	if (Swap_Left == 1)
	{
		ConeInfo	*Temp_ConeInfo_Pointer;
		
		Temp_ConeInfo_Pointer = Temp_ConeInfo_Heap_Node->ConeInfo_Pointer;
		
		Temp_ConeInfo_Heap_Node->ConeInfo_Pointer = Temp_ConeInfo_Heap_Node->Left->ConeInfo_Pointer;

		Temp_ConeInfo_Heap_Node->Left->ConeInfo_Pointer = Temp_ConeInfo_Pointer;
	
		Restore_Down (Temp_ConeInfo_Heap_Node->Left);
	}
	if (Swap_Right == 1) //swap right and parent
	{
			
		
		ConeInfo	*Temp_ConeInfo_Pointer;
		
		Temp_ConeInfo_Pointer = Temp_ConeInfo_Heap_Node->ConeInfo_Pointer;
		
		Temp_ConeInfo_Heap_Node->ConeInfo_Pointer = Temp_ConeInfo_Heap_Node->Right->ConeInfo_Pointer;

		Temp_ConeInfo_Heap_Node->Right->ConeInfo_Pointer = Temp_ConeInfo_Pointer;
	
		Restore_Down (Temp_ConeInfo_Heap_Node->Right);
	}
}


void ConeInfo_Heap::Restore_Up (ConeInfo_Heap_Node *Temp_ConeInfo_Heap_Node)
{
	if (Temp_ConeInfo_Heap_Node->Parent != NULL)
	{
		if ( *(Temp_ConeInfo_Heap_Node->ConeInfo_Pointer->Get_Current_Highest_Term () ) > *(Temp_ConeInfo_Heap_Node->Parent->ConeInfo_Pointer->Get_Current_Highest_Term () ) )
		{
			//cout << "Restore_Up: swaping" << endl;			
			ConeInfo	*Temp_ConeInfo;
			

			Temp_ConeInfo = Temp_ConeInfo_Heap_Node->ConeInfo_Pointer;

			Temp_ConeInfo_Heap_Node->ConeInfo_Pointer = Temp_ConeInfo_Heap_Node->Parent->ConeInfo_Pointer;

			Temp_ConeInfo_Heap_Node->Parent->ConeInfo_Pointer = Temp_ConeInfo;

			Restore_Up(Temp_ConeInfo_Heap_Node->Parent);
		}
		//else
			//cout << "No swap" << endl;
	}
	//else
		//cout << "Parent Null" << endl;
}

int ConeInfo_Heap::Check_Top_Heap(ConeInfo *Temp_ConeInfo)
{
	if (Root_Node != NULL)
	{
		if ( *(Temp_ConeInfo->Get_Current_Highest_Term ()) == *(Root_Node->ConeInfo_Pointer->Get_Current_Highest_Term () ) )
			return 1;
		else
			return 0;
	}
	else //list is empty
	{
		return 0;
	}

}

ZZ	*ConeInfo::Get_Current_Highest_Term ()
{
	return	Current_Highest_Term;
}

int	ConeInfo::Get_Coefficient ()
{
	return Coefficient;
}

int	ConeInfo::Calculate_Integral_Point (vector &Temp_Vector)
{
	if(Vector_List_Head == NULL)
		return 0;

	Temp_Vector = listCone_ptr->latticePoints->first;
	listVector *temp;

	temp = listCone_ptr->rays;

	for(int i=0; i < Number_of_Generators; i++)
	{
		if(signs[i] > 0)
		{
			Temp_Vector -= temp->first;
		}
		temp = temp->rest;
	}

	temp = listCone_ptr->rays;
	int temp_coefficients[Number_of_Generators];

	for(int i =0; i< Number_of_Generators; i++)
		temp_coefficients[Order[i]] = Vector_List_Head->Integer_Array[i];
	
	for(int i=0; i < Number_of_Generators; i++)
	{
		for(int j = 0; j < Number_of_Variables; j++)
		{
			// sign[] is the sign of the orginal dot product with the cost
			Temp_Vector[j] -= signs[i] * temp->first[j] * temp_coefficients[i];
		}
		temp = temp->rest;
		
	}	

	Integer_Vector_List *del_list = Vector_List_Head;
	Vector_List_Head = Vector_List_Head->Next;
	delete [] del_list->Integer_Array;
	delete del_list;

	return 1;
}

int	ConeInfo::Hash_Integer_Vector (int *Integer_Array)
{
	int	Return_Value = 0;	

	for (int i = 0; i < Number_of_Generators; i++)
	{
		Return_Value += Hash_Function_Coefficients[i] * Integer_Array[i];
	
	}
		
	Return_Value = Return_Value % HASH_TABLE_SIZE;
	
	if (Return_Value >= 0)
		return Return_Value;
	else
		return (-1)*Return_Value;
}

vector	Calculate_Pertubation (listCone *cones, vector *Cost, int Mod_Value, int Number_of_Variables)
{
	/* OLD CODE: 
	//cout << "Calculate_Pertubation: " << (*Cost)*(*Cost) << endl;
	vector	Return_Vector;

	Return_Vector = *Cost;

	//srand(clock());

	//unsigned long int	Integral_Factor = 1000;

	//ZZ	Normalize_Length;

	//Normalize_Length = 64000000;

	//Integral_Factor = rand ();

	ZZ	Root,N;
	
	Root = 1;

	N = Normalize_Length/((*Cost)*(*Cost));

	for (int i = 0; i < 100; i++)
	{
		Root = ((Root + N/Root)/2);
		if(Root == 0)
			Root = 1;
	}
		
	Root += 1;
	
	*/
	// NEW CODE: Dave + Pete 10/3/03
	vector	Return_Vector;
	Return_Vector = *Cost;

	listCone *tmpcone;
	listVector *tmpVector;
	ZZ	Maximum_Coordinate;
	
	Maximum_Coordinate = 0;
	
	tmpcone = cones;
	
	srand(clock());

	while (tmpcone)
	{
		tmpVector = tmpcone->rays;

		while (tmpVector)
		{
			for (int i=0; i < Number_of_Variables; i++)
			{
				if (tmpVector->first[i] >= 0 && tmpVector->first[i] > Maximum_Coordinate)
					Maximum_Coordinate = tmpVector->first[i];
				if (tmpVector->first[i] <= 0 && -1* tmpVector->first[i] > Maximum_Coordinate)
					Maximum_Coordinate = tmpVector->first[i];

		
			}
			
			tmpVector = tmpVector->rest;
		}	

		tmpcone = tmpcone->rest;
	}
	
	ZZ	M_Value;

	M_Value = Mod_Value * Maximum_Coordinate * Number_of_Variables;
	
	// Take 2M + 1
	M_Value *= 2;
	M_Value += 1;

	
	
	for (int i = 0; i < Number_of_Variables; i++)
	{
		Return_Vector[i] = M_Value * Return_Vector[i];

	}	
	
	for (int i = 0; i < Number_of_Variables; i++)
	{
		Return_Vector[i] += (rand() % Mod_Value)* (( rand () % 2)*2 - 1);

	}	
	
	return Return_Vector;
}


Vector_Heap_Array_Node_Controller::Vector_Heap_Array_Node_Controller (int Initial_Integer_Array_Size)
{
	Integer_Array_Size = Initial_Integer_Array_Size;

	// Have an Integer_Vector_List ready to go
	Head_Integer_Vector_List = new Integer_Vector_List;

	Head_Integer_Vector_List->Next = NULL;
	Head_Integer_Vector_List->Integer_Array = new int[Integer_Array_Size];
	
	// Have a ZZ ready to go
	Head_ZZ_List = new ZZ_List;

	Head_ZZ_List->Next = NULL;	
	Head_ZZ_List->ZZ_Value = new ZZ;

}

Vector_Heap_Array_Node_Controller::~Vector_Heap_Array_Node_Controller ()
{
	Delete_Lists ();
	
}

void Vector_Heap_Array_Node_Controller::Delete_Lists ()
{
	Integer_Vector_List	*Temp_Integer_Vector_List, *Delete_Integer_Vector_List;
	ZZ_List			*Temp_ZZ_List, *Delete_ZZ_List;

	Temp_Integer_Vector_List = Head_Integer_Vector_List;

	while (Temp_Integer_Vector_List != NULL)
	{
		delete [] Temp_Integer_Vector_List->Integer_Array;
		Delete_Integer_Vector_List = Temp_Integer_Vector_List;

		Temp_Integer_Vector_List = Temp_Integer_Vector_List->Next;

		delete Delete_Integer_Vector_List;
		
	}

	Temp_ZZ_List = Head_ZZ_List;

	while (Temp_ZZ_List != NULL)
	{
		delete Temp_ZZ_List->ZZ_Value;
		Delete_ZZ_List = Temp_ZZ_List;
		
		Temp_ZZ_List = Temp_ZZ_List->Next;

		delete Delete_ZZ_List;
		
	}

}
	
int	Vector_Heap_Array_Node_Controller::Get_Current_Integer_Array_Size ()
{
	return Integer_Array_Size;
}
	
int	*Vector_Heap_Array_Node_Controller::Get_Integer_Array ()
{
	int	*Return_Value;

	Return_Value = Head_Integer_Vector_List->Integer_Array;

	Head_Integer_Vector_List = Head_Integer_Vector_List->Next;

	// check if we need to create more integer arrays
	if (Head_Integer_Vector_List == NULL)
	{
		Head_Integer_Vector_List = new Integer_Vector_List;

		Head_Integer_Vector_List->Next = NULL;
		Head_Integer_Vector_List->Integer_Array = new int [Integer_Array_Size];
		
	}
		
	return Return_Value;
}

ZZ	*Vector_Heap_Array_Node_Controller::Get_ZZ ()
{
	ZZ	*Return_Value;

	Return_Value = Head_ZZ_List->ZZ_Value;

	Head_ZZ_List = Head_ZZ_List->Next;

	// check if we need to create more ZZ's
	if (Head_ZZ_List == NULL)
	{
		Head_ZZ_List = new ZZ_List;
		Head_ZZ_List->Next = NULL;

		Head_ZZ_List->ZZ_Value = new ZZ;
	}

	return Return_Value;
}

void	Vector_Heap_Array_Node_Controller::Recieve_Integer_Array (int *Used_Integer_Array)
{
	Integer_Vector_List	*Temp_Integer_Vector_List;
	
	//Find the end of the list and add to it
	
	Temp_Integer_Vector_List = new Integer_Vector_List;

	Temp_Integer_Vector_List->Next = Head_Integer_Vector_List;

	Temp_Integer_Vector_List->Integer_Array = Used_Integer_Array;

	Head_Integer_Vector_List = Temp_Integer_Vector_List;	
	
}

void	Vector_Heap_Array_Node_Controller::Recieve_ZZ (ZZ *Used_ZZ_Value)
{
	ZZ_List			*Temp_ZZ_List;
	
	Temp_ZZ_List = new ZZ_List;
	
	Temp_ZZ_List->Next = Head_ZZ_List;


	Temp_ZZ_List->ZZ_Value = Used_ZZ_Value;
	
	Head_ZZ_List = Temp_ZZ_List;
}
		

