#ifndef RUDYRESNTL__H
#define RUDYRESNTL__H
#include "PolyTree.h"
#include "myheader.h"
#include "cone.h"
#include "ramon.h"

struct Single_Cone_Parameters
{
	listCone	*Cone;
	int		Number_of_Variables;
	int		Degree_of_Taylor_Expansion;
	unsigned int	Flags;
	
	ZZ		*Taylor_Expansion_Result;
	ZZ		*Random_Lambda;
	ZZ		*Ten_Power;
	ZZ		*Total_Lattice_Points;
	ZZ		*Total_Uni_Cones;
	ZZ		*Current_Simplicial_Cones_Total;
	ZZ		*Max_Simplicial_Cones_Total;
};


void ResidueFunction(listCone* cones, int numOfVars, int print_flag, int degree, int output_cone);

// Returns -1 if a Dot Product is zero in the denominator, otherwise 1 if ok
int	ResidueFunction_Single_Cone (Single_Cone_Parameters *Residue_Parameters, Node_Controller *Controller); 

#endif

