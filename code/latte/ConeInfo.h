#ifndef CONEINFO_H
#define CONEINFO_H

#include <NTL/ZZ.h>
#include "myheader.h"
#include "ramon.h"

#define HASH_TABLE_SIZE 1247677

vec_ZZ	Calculate_Pertubation (listCone *cones, vec_ZZ *Cost, int Mod_Value, int Number_of_Variables);


class ConeInfo;

struct Hash_Table_Entry
{
	int		Time_Stamp;
	int		*Integer_Array;

};

struct ConeInfo_List
{
	ConeInfo	*ConeInfo_Pointer;
	
	ConeInfo_List	*Next;
};

struct Integer_Vector_List
{
	Integer_Vector_List	*Next;
	int	*Integer_Array;
};

struct Heap_Node
{
	Heap_Node	*Parent, *Left, *Right;
		//Sum_Vector is an array of size Number_of_Generators
	int	*Sum_Vector;
	ZZ	*Total_Sum;	
};

struct ZZ_List
{
	ZZ	*ZZ_Value;
	ZZ_List	*Next;
};

class Vector_Heap_Array_Node_Controller
{
	public:
		Vector_Heap_Array_Node_Controller (int Initial_Integer_Array_Size);	
		~Vector_Heap_Array_Node_Controller ();
		
		int	*Get_Integer_Array ();
		ZZ	*Get_ZZ ();
		void	Recieve_Integer_Array (int *Used_Integer_Array);
		void	Recieve_ZZ (ZZ *Used_ZZ_Value);
		int	Get_Current_Integer_Array_Size ();	
	private:
		void Delete_Lists ();
		Integer_Vector_List 	*Head_Integer_Vector_List;	
		ZZ_List			*Head_ZZ_List;
		int	Integer_Array_Size;
};

class Vector_Heap 
{
	public:
		Vector_Heap (int Initial_Number_of_Generators);
		~Vector_Heap ();
			//Returns what is on top of the heap WITHOUT removing it.  Places results in sum_vector, total_sum	
		int	Get_Top_Heap (int Sum_Vector [], ZZ *Total_Sum);
			//Pops off the top of the heap and returns its value
		int	Pop_Top_Heap (int Sum_Vector [], ZZ *Total_Sum);
			//Add this value to the heap
		void 	Add_Heap (int *Sum_Vector, ZZ *Total_Sum);
			//Checks if the top of the heap is equal to the input data, 1 true, 0 false
		int	Check_Top_Heap (ZZ *Total_Sum);

		int	Get_Node_Count ();	
		void	Print_Tree ();		
	private:
		void	Print_Sub_Tree (Heap_Node *Node);
		//	Restores the structure on the subtree of Node
		void	Restore_Up (Heap_Node *Node);
		void	Restore_Down (Heap_Node *Node);
		void	Delete_Sub_Tree (Heap_Node *Node);
		Heap_Node *Root_Node;
		unsigned int	Node_Count;
		int	Number_of_Generators;
		static	Vector_Heap_Array_Node_Controller *Controller;
};
		
class ConeInfo
{
	public:
		ConeInfo (vec_ZZ *cost, listCone *the_cone, int numOfVars);
		~ConeInfo ();
		ZZ	*Get_Current_Highest_Term ();
		int	Get_Coefficient ();
		void	Calculate_Next_Highest_Term ();
		int	Calculate_Integral_Point (vec_ZZ &);
	//	Vector_Heap	*Heap;

		int	S_Values_Zero_Flag; // 0 not zero, 1 some value is zero
		static int	Get_Object_Count();
		static int	Get_Time_Stamp ();
	private:
		void	Sort_S_Values (); 
		listCone 	*listCone_ptr;
		ZZ		*S_Values;
		ZZ		*Expansion_Highest_Term;	// purely in terms expansion
		ZZ		*Numer_Exp;		// factored out exponent
		ZZ		*Current_Highest_Term;		// Addition of expansion and numerator
		Integer_Vector_List	*Vector_List_Head;
		unsigned int	Current_Highest_Term_Coefficient;
		int		*signs;
		int		Coefficient;
		int		*Order;
		Vector_Heap	*Heap;
		int		Number_of_Variables;
		int		Number_of_Generators;
		
		//  HASH TABLE IS LINEAR PROBING.
		int		Hash_Integer_Vector (int *);		
		
		static int	Hash_Table_Initialized_Flag; //0 not init, 1 init
		static int	Time_Stamp;
		static int	Object_Count;
		static int	*Hash_Function_Coefficients;
		static Hash_Table_Entry	Hash_Table[HASH_TABLE_SIZE];		

};		

struct ConeInfo_Heap_Node
{
	ConeInfo_Heap_Node	*Parent, *Left, *Right;
	ConeInfo	*ConeInfo_Pointer;
};

class ConeInfo_Heap 
{
	public:
		ConeInfo_Heap ();
		~ConeInfo_Heap ();
			//Returns what is on top of the heap WITHOUT removing it.  Places results in sum_vector, total_sum	
		ConeInfo	*Get_Top_Heap ();
			//Pops off the top of the heap and returns its value
		ConeInfo	*Pop_Top_Heap ();
			//Add this value to the heap
		void 	Add_Heap (ConeInfo	*Temp_Cone_Pointer);
			//Checks if the top of the heap is equal to the input data, 1 true, 0 false
		int	Check_Top_Heap (ConeInfo *);
		int	Get_Node_Count ();	
		void	Clear_Tree ();
		//void	Print_Tree ();		
	private:
		//void	Print_Sub_Tree (ConeInfo_Heap_Node *Node);
		//	Restores the structure on the subtree of Node
		void	Restore_Up (ConeInfo_Heap_Node *Node);
		void	Restore_Down (ConeInfo_Heap_Node *Node);
		void	Delete_Sub_Tree (ConeInfo_Heap_Node *Node);
		ConeInfo_Heap_Node *Root_Node;
		unsigned int	Node_Count;
};

#endif


