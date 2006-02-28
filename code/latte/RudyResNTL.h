#ifndef RUDYRESNTL__H
#define RUDYRESNTL__H
#include "PolyTree.h"
#include "myheader.h"
#include "cone.h"
#include "ramon.h"

// Later we will reduce the slots in this class and use further
// subclassing for the individual computation modes.  For instance,
// Taylor_Expansion_Result is only used in the "dual" method.
// -- mkoeppe

class Single_Cone_Parameters {
 public:
	listCone	*Cone; // The master cone to be decomposed.
	ZZ		*Total_Uni_Cones;
	ZZ		*Current_Simplicial_Cones_Total;
	ZZ		*Max_Simplicial_Cones_Total;
	int		Number_of_Variables;
	unsigned int	Flags;
	virtual int ConsumeCone(listCone *cone) = 0;
};

class Standard_Single_Cone_Parameters : public Single_Cone_Parameters {
 public:
	int		Degree_of_Taylor_Expansion;
	
	ZZ		*Taylor_Expansion_Result;
	ZZ		*Random_Lambda;
	ZZ		*Ten_Power;
	ZZ		*Total_Lattice_Points;

	Node_Controller *Controller;
 public:
	virtual int ConsumeCone(listCone *cone);
};



void ResidueFunction(listCone* cones, int numOfVars, int print_flag, int degree, int output_cone);

// Returns -1 if a Dot Product is zero in the denominator, otherwise 1 if ok
int
ResidueFunction_Single_Cone (listCone *cones,
			     Standard_Single_Cone_Parameters *Residue_Parameters); 

#endif

