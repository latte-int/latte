#ifndef POLYTREE_H
#define POLYTREE_H 1

#include <NTL/ZZ.h>
#include <fstream>
using namespace NTL_NAMESPACE;
using namespace std;

// Definitions of the types of PolyTree structures
#define POLYTREE_ADD 0
#define POLYTREE_MUL 1
#define POLYTREE_DIV 5
#define POLYTREE_T_NODE 3
#define POLYTREE_S_NODE 2
#define POLYTREE_EXP 4

// Define the types of generators of cones,  either (1-rt) or (r-t)
#define ONE_SUB_RT 0
#define R_SUB_T  1

// This is used to hold all the information of Taylor_Expansion function
struct Taylor_Parameters
{
	ZZ 	*Result;
	int	Degree_of_Expansion;
	ZZ	*Ten_Power;

};

class PolyTree_Node
{
	public:
		PolyTree_Node () { Taylor_Expansion_Result_Dirty = 1; };
		
		PolyTree_Node **Children;
		
		// 0 mean +, 1 mean *, 2 means S_Node, 3 means T_Node, 4 mean ^ , 5 means / 
		char	Node_Type;
		unsigned int	Number_of_Children;
		virtual int 	Print(); //Returns 1 if it did print anything, otherwise returns 0
		virtual int	Check_For_Zero ();  //Returns 1 if this node is equal to zero
		
		//This function will return the first Parameters->Degree_of_Expansion number of terms.  
		//It accepts  Taylor_Parameter as input with the Result part of it already
		//allocated for enough space for Degree_of_Expansion + 1 terms.
		virtual void	Taylor_Expansion (Taylor_Parameters *Parameters);
		virtual int	Print_Rational_Functions_to_File (ofstream &Output_File);
	
		//		1 indicates Taylor_Expansion_Result is invalid
		int		Taylor_Expansion_Result_Dirty;
		
		ZZ		*Taylor_Expansion_Result;
};


//Currently this class is never used
/*class S_Node : public PolyTree_Node
{
	public:
		ZZ	Coefficient;
		ZZ	Exponent;
		int	Print ();
};*/

class T_Node : public PolyTree_Node //Also used for ^ type.  Exponent holds the exponent. 
{
	public:
		virtual ~T_Node () {};
		ZZ	Coefficient;
		ZZ	Exponent;
		int	Print(); //Returns 1 if it did print anything, otherwise returns 0
		int	Check_For_Zero (); // Returns 1 if zero, else returns 0
		
		//Same as description in PolyTree_Node
		void	Taylor_Expansion (Taylor_Parameters *Parameters);
		int	Print_Rational_Functions_to_File (ofstream &Output_File);
};


// Used to hold information for each generator of a cone.
struct Generator
{
	int	Form_Type;  // 0 means (1-rt)     1 means  (r-t)
	ZZ	R_Exponent;
	ZZ	T_Exponent;

};


struct Cone_Data 
{
	int	sign;	//Holds the sign of this Cone	
	int	order;	//Holds the order of this Cone
	Generator	*Generators_of_Cone;  //The denominators
	Generator	Numerator_Generator;
	
	
	
};

struct PolyTree_Node_List
{
	PolyTree_Node	*Data;
	
	PolyTree_Node_List *Next;

};

struct T_Node_List
{
	T_Node		*Data;

	T_Node_List 	*Next;

};

class Node_Controller
{
	public:
		Node_Controller (int Dimension, int Degree);
		~Node_Controller ();
			
		//  Returns a pointer to a PolyNode for someone to use
		PolyTree_Node 	*Get_PolyTree_Node ();
		//  Returns a pointer to a T_Node for someone to use
		T_Node		*Get_T_Node ();

		// Clears all the nodes this Controller controls
		void	Reset ();

	private:
		int	Dimension_Plus_One;
		int	Degree_of_Expansion;
			

		PolyTree_Node_List	*PolyTree_Node_Head;
		PolyTree_Node_List	*PolyTree_Node_Unused;

		T_Node_List		*T_Node_Head;
		T_Node_List		*T_Node_Unused;

};



#endif
