/* This is a -*- C++ -*- header file. */
#ifndef FLAGS_H
#define FLAGS_H 1

#include "myheader.h"   // FIXME: For ListCone
#include "PolyTree.h"  // FIXME: For NodeController

#define PRINT 0x1
#define OUTPUT 0x6
#define OUTPUT0 0x2
#define OUTPUT1 0x4
#define DUAL_APPROACH 0x8
#define DECOMPOSE 0x10
#define LOAD	0x20
#define SAVE    0x40

// FIXME: This file (flags.h) is certainly not the right place yet for
// these class definitions.  --mkoeppe

struct BarvinokParameters {
  // Whether we use the
  //   - traditional LattE monomial substitution z_i |-> (1 + s)^(lambda_i) 
  //   - or the exponential substitution         z_i |-> exp(t lambda_i)
  enum { PolynomialSubstitution, ExponentialSubstitution } substitution;
  // The maximum determinant of cones that we do not subdivide
  // further.  Set to 1 to subdivide until we reach unimodular cones
  // only.
  int max_determinant;
  // A file name that is used for constructing file names for
  // temporary and semi-temporary files.
  char *File_Name;
  // Ambient dimension.
  int Number_of_Variables;
};

class Single_Cone_Parameters : public BarvinokParameters {
public:
  // Parameters that control the computation.
  unsigned int	Flags;
public:
  // Data that are used during the computation.
  //listCone	*Cone;		// The master cone to be decomposed.
  int		Cone_Index;	/* Its index in the list of all master
				   cones; only used for naming
				   triangulation caches. */
public:
  // Statistics collected during the computation.
  ZZ		Total_Uni_Cones;
  ZZ		Current_Simplicial_Cones_Total;
  ZZ		Max_Simplicial_Cones_Total;
public:
  virtual int ConsumeCone(listCone *cone) = 0;
};

// The traditional LattE mode: Simply collect all subdivided cones
// into a list.
class Collecting_Single_Cone_Parameters : public Single_Cone_Parameters {
public:
  Collecting_Single_Cone_Parameters();
  listCone *Decomposed_Cones;
  virtual int ConsumeCone(listCone *cone);
};

// The "Memory Save" mode: Perform residue calculations immediately
// at each subdivided cone in the tree, don't store the cones.
//
// FIXME: Later we will reduce the slots in this class and use further
// subclassing for the individual computation modes.  For instance,
// Taylor_Expansion_Result is only used in the "dual" method.  -- mkoeppe

class Standard_Single_Cone_Parameters : public Single_Cone_Parameters {
 public:
	int		Degree_of_Taylor_Expansion;
	
	ZZ		*Taylor_Expansion_Result;
	ZZ		*Random_Lambda;
	ZZ		Ten_Power;
	ZZ		Total_Lattice_Points;

	Node_Controller *Controller;
 public:
	virtual int ConsumeCone(listCone *cone);
};



#endif
