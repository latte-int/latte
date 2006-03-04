/*********************************************************** -*- C++ -*- */
#ifndef BARVINOK_DEC_H
#define BARVINOK_DEC_H

#include "../flags.h"
#include "PolyTree.h"
#include "barvinok.h"

// FIXME: Move somewhere else
listCone* readListCone(rationalVector*, int);


// The traditional LattE mode: Simply collect all subdivided cones
// into a list.
class Collecting_Single_Cone_Parameters : public Single_Cone_Parameters {
public:
  Collecting_Single_Cone_Parameters();
  listCone *Decomposed_Cones;
  virtual int ConsumeCone(listCone *cone);
};

listCone*
decomposeCones(listCone*, int, unsigned int Flags,
	       char *File_Name, int max_determinant = 1);


/* Guess a generic vector; this is simply a random vector. */
vec_ZZ
guess_generic_vector(int numOfVars);

/* Functions can throw this exception when they discover the
   passed vector was not generic. */
struct NotGenericException {};

class Generic_Vector_Single_Cone_Parameters : public Single_Cone_Parameters {
public:
  vec_ZZ generic_vector;
  virtual void InitializeComputation();
};

/* Pick a tentative generic vector by calling InitializeComputation().
   Then call barvinokDecomposition_Single on all CONES.  This results
   in ConsumeCone() being called on all resulting small cones; when
   any ConsumeCone() call returns with -1, restart the computation by
   calling InitializeComputation() and start decomposing again. */
void
barvinokDecomposition_List(listCone *cones,
			   Generic_Vector_Single_Cone_Parameters &Parameters);


// The "Memory Save" mode: Perform residue calculations immediately
// at each subdivided cone in the tree, don't store the cones.
//
// FIXME: Later we will reduce the slots in this class and use further
// subclassing for the individual computation modes.  For instance,
// Taylor_Expansion_Result is only used in the "dual" method.  -- mkoeppe

class Standard_Single_Cone_Parameters
  : public Generic_Vector_Single_Cone_Parameters {
 public:
	int		Degree_of_Taylor_Expansion;
	
	ZZ		*Taylor_Expansion_Result;
	ZZ		Ten_Power;
	ZZ		Total_Lattice_Points;

	Node_Controller *Controller;
 public:
	virtual int ConsumeCone(listCone *cone);
};

void decomposeCones_Single (listCone *, int, int degree,
			    unsigned int flags, char *File_Name);

#endif
