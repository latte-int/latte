/***********************************************************************
  Author: Ruriko Yoshida
  July 24th, 2002
  Update: Febrary 3rd, 2003
  This is a program for Barvinok's decomposition of cones.
  This is a class file.

************************************************************************/
#ifndef BARVINOK__H
#define BARVINOK__H
#include <list>

#include <fstream>
#include <cstdlib>
#include <cstring>
#include <math.h>
#include <algorithm>
#include <time.h>
#include "../myheader.h"
#include "../ramon.h"
#include "Cone.h"
#include <NTL/mat_ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_RR.h>
#include <NTL/vec_RR.h>
#include <NTL/config.h>
#include <NTL/LLL.h>
#include <NTL/HNF.h>
#include <NTL/ZZ.h>
#include "../RudyResNTL.h"
#include "../flags.h"
using namespace std;

struct Barvinok_DFS_Parameters
{
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
	
	rationalVector *vertex;
	
	Node_Controller *Controller;
};


int barvinok(mat_ZZ &, list< PtrCone > &, int&);

int barvinok_Single (mat_ZZ &, int &, Single_Cone_Parameters *, Node_Controller *, rationalVector *vertex);

int barvinok_DFS(Cone *, Barvinok_DFS_Parameters *);


#endif
