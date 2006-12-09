// This is a -*- C++ -*- header file.

/* barvinok.h -- Barvinok's decomposition of a cone.

   Copyright 2002, 2003 Ruriko Yoshida
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

#ifndef BARVINOK__H
#define BARVINOK__H

#include "cone.h"
#include "timing.h"

class BarvinokParameters {
public:
  // Whether we use the
  //   - traditional LattE monomial substitution z_i |-> (1 + s)^(lambda_i) 
  //   - or the exponential substitution         z_i |-> exp(t lambda_i)
  enum { PolynomialSubstitution, ExponentialSubstitution } substitution;
  // Whether to use
  //  - traditional dual decomposition
  //  - irrational primal decomposition
  //  - irrational all-primal decomposition
  //    (i.e., both triangulation and Barvinok decomposition in primal space)
  typedef enum {
    DualDecomposition,
    IrrationalPrimalDecomposition,
    IrrationalAllPrimalDecomposition
  } DecompositionType;
  DecompositionType decomposition; 
  // The kind of triangulation to use.
  typedef enum {
    RegularTriangulationWithCdd,
    RegularTriangulationWithCddlib,
    SubspaceAvoidingRecursiveTriangulation,
    SubspaceAvoidingBoundaryTriangulation,
    PlacingTriangulationWithTOPCOM
  } TriangulationType;
  TriangulationType triangulation;
  // The kind of short vectors we use for decomposition
  typedef enum {
    LatteLLL,
    SubspaceAvoidingLLL    
  } ShortVectorType;
  ShortVectorType shortvector;
  // The maximum determinant of cones that we do not subdivide
  // further.  Set to 1 to subdivide until we reach unimodular cones
  // only.  Set to 0 (special case) to not subdivide at all. 
  int max_determinant;
  // A file name that is used for constructing file names for
  // temporary and semi-temporary files.
  char *File_Name;
  // Ambient dimension.
  int Number_of_Variables;
  // Parameters that control the computation.
  unsigned int Flags;
  // Data that are used during the computation.
  int		Cone_Index;	/* Its index in the list of all master
				   cones; only used for naming
				   triangulation caches. */
  // Timers.
  Timer total_time;
  Timer read_time;
  Timer vertices_time;
  Timer irrationalize_time;
  Timer dualize_time;
  Timer triangulate_time;
  Timer decompose_time;
  // Constructor & destructor.
  BarvinokParameters();
  virtual ~BarvinokParameters();
  virtual void print_statistics(ostream &s);
};

class Single_Cone_Parameters : public BarvinokParameters {
public:
  // Statistics collected during the computation.
  ZZ		Total_Uni_Cones;
  ZZ		Total_Simplicial_Cones;
  ZZ		Current_Simplicial_Cones_Total;
  ZZ		Max_Simplicial_Cones_Total;
  int		Current_Depth;
  int		Max_Depth;
public:
  Single_Cone_Parameters() : Current_Depth(0), Max_Depth(0) {};
  Single_Cone_Parameters(const BarvinokParameters &params) :
    BarvinokParameters(params), Current_Depth(0), Max_Depth(0) {};
  virtual int ConsumeCone(listCone *cone) = 0;
  virtual ~Single_Cone_Parameters() {}
  virtual void print_statistics(ostream &s);
};

/* Do a signed decomposition, modulo lower-dimensional cones, of the
   SIMPLICIAL cone spanned by the ROW VECTORS of B with apex at
   VERTEX, until the determinants of all cones are at most
   PARAMETERS->max_determinant.

   Call PARAMETERS->ConsumeCone() for each of the small cones.
*/ 
int
barvinok_Single(mat_ZZ B, Single_Cone_Parameters *Parameters,
		const Vertex *vertex);

/* Do a signed decomposition, modulo lower-dimensional cones, of the
   polyhedral CONE, until the determinants of all cones are at most
   PARAMETERS->max_determinant.  (When the cone is not simplicial, we
   triangulate it first.)

   Call PARAMETERS->ConsumeCone() for each of the small cones.
*/ 
int
barvinokDecomposition_Single (listCone *cone,
			      Single_Cone_Parameters *Parameters);

#endif
