/*********************************************************** -*- C++ -*-
  Author: Ruriko Yoshida
  July 24th, 2002
  Update: Febrary 3rd, 2003
  This is a program for Barvinok's decomposition of cones.
  This is a class file.

************************************************************************/
#ifndef BARVINOK__H
#define BARVINOK__H

#include "cone.h"

struct BarvinokParameters {
  // FIXME: Following does not really belong here.
  // Whether we use the
  //   - traditional LattE monomial substitution z_i |-> (1 + s)^(lambda_i) 
  //   - or the exponential substitution         z_i |-> exp(t lambda_i)
  enum { PolynomialSubstitution, ExponentialSubstitution } substitution;
  // Whether to use
  //  - traditional dual decomposition
  //  - irrational primal decomposition
  typedef enum { DualDecomposition, IrrationalPrimalDecomposition } DecompositionType;
  DecompositionType decomposition; 
  // The maximum determinant of cones that we do not subdivide
  // further.  Set to 1 to subdivide until we reach unimodular cones
  // only.  Set to 0 (special case) to not subdivide at all. 
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
  ZZ		Total_Simplicial_Cones;
  ZZ		Current_Simplicial_Cones_Total;
  ZZ		Max_Simplicial_Cones_Total;
public:
  virtual int ConsumeCone(listCone *cone) = 0;
  virtual ~Single_Cone_Parameters() {}
};

/* Do a signed decomposition, modulo lower-dimensional cones, of the
   SIMPLICIAL cone spanned by the ROW VECTORS of B with apex at
   VERTEX, until the determinants of all cones are at most
   PARAMETERS->max_determinant.

   Call PARAMETERS->ConsumeCone() for each of the small cones.
*/ 
int
barvinok_Single(mat_ZZ B, Single_Cone_Parameters *Parameters,
		rationalVector *vertex);

#endif
