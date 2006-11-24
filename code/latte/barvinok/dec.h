// This is a -*- C++ -*- header file.

/* Cone.cpp -- Barvinok's decomposition of a cone.

   Copyright 2002, 2003 Raymond Hemmecke, Ruriko Yoshida
   Copyright 2006 Matthias Koeppe

   This file is part of LattE.
   
   LattE is free software; you can redistribute it and/or modify it
   under the terms of the version 2 of the GNU General Public License
   as published by the Free Software Foundation.

   LattE is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with LattE; if not, write to the Free Software Foundation,
   Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
*/

#ifndef BARVINOK_DEC_H
#define BARVINOK_DEC_H

#include "../flags.h"
#include "PolyTree.h"
#include "barvinok.h"


// The traditional LattE mode: Simply collect all subdivided cones
// into a list.
class Collecting_Single_Cone_Parameters : public Single_Cone_Parameters {
public:
  Collecting_Single_Cone_Parameters();
  Collecting_Single_Cone_Parameters(const BarvinokParameters &params); 
  listCone *Decomposed_Cones;
  virtual int ConsumeCone(listCone *cone);
};

// Obsolete:
listCone*
decomposeCones(listCone *cones, int numOfVars, unsigned int Flags,
	       char *File_Name, int max_determinant,
	       bool dualize,
	       BarvinokParameters::DecompositionType decomposition);
// Nicer, more general interface:
listCone*
decomposeCones(listCone *cones, bool dualize,
	       BarvinokParameters &param);


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
  Generic_Vector_Single_Cone_Parameters() {};
  Generic_Vector_Single_Cone_Parameters(const BarvinokParameters &params) :
    Single_Cone_Parameters(params) {};
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
  Standard_Single_Cone_Parameters() {};
  Standard_Single_Cone_Parameters(const BarvinokParameters &params) :
    Generic_Vector_Single_Cone_Parameters(params) {};
  virtual void InitializeComputation();
  virtual int ConsumeCone(listCone *cone);
};

// Decompose the polyhedral CONES down to MAX_DETERMINANT.  Then
// perform residue calculations and print results.  When DUALIZE is
// true, the CONES are given in primal space, so dualize before
// triangulating; otherwise CONES must be given in dual space already.
void
decomposeAndComputeResidue(listCone *cones, int degree, bool dualize,
			   Standard_Single_Cone_Parameters &param);

// Likewise, deprecated interface.  
void decomposeCones_Single (listCone *cones, int numOfVars, int degree,
			    unsigned int flags, char *File_Name,
			    int max_determinant,
			    bool dualize,
			    BarvinokParameters::DecompositionType decomposition);

#endif
