// This is a -*- C++ -*- header file.

/* ReadPolyhedron.h -- Handle command-line args to read a polyhedron
	       
   Copyright 2007 Matthias Koeppe

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

#ifndef LATTE_READPOLYHEDRON_H
#define LATTE_READPOLYHEDRON_H

#include "Polyhedron.h"
#include "barvinok/barvinok.h"

class ReadPolyhedronData {
public:
  // A maze of twisty parameters, all alike
  char equationsPresent[10];
  char nonneg[127];
  char cddstyle[127];
  char Vrepresentation[127];
  char dilation[127];
  char interior[127];
  int dilation_const;
  char dualApproach[127];
  string filename;
  char Memory_Save[127];
  char inthull[127];
  char grobner[127];
  char maximum[127];
  char minimize[127];
  char taylor[127];
  char rationalCone[127];
  char assumeUnimodularCones[127];
  char Singlecone[127];
  int degree;
public:
  // Data for input of cones
  bool input_homog_cone, input_dualized, have_subcones, input_listcone_format;
  string subcones_filename;
public:
  // How to compute vertex cones.
  typedef enum {
    VertexConesWithCdd,
    VertexConesWithLrs,
    VertexConesWith4ti2
  } VertexConesType;
  VertexConesType vertexcones;
  // How to obtain a non-redundant representation.
  typedef enum {
    RedundancyCheckWithCdd,
    NoRedundancyCheck,
    RedundancyCheckWithCddlib
  } RedundancyCheckType;
  RedundancyCheckType redundancycheck;
public:
  // A maze of twisty intermediate data, all alike.
  vec_ZZ cost;
  listVector *matrix;  // Sometimes the original matrix.
  mat_ZZ AA;			// Data related
  vec_ZZ bb;			// to un-projection.
  int oldnumofvars;		// 
  listVector *templistVec; 	// 
public:
  ReadPolyhedronData();
  void show_options(ostream &stream);
  bool parse_option(const char *arg);
  Polyhedron *read_polyhedron(BarvinokParameters *params);
private:
  Polyhedron *read_polyhedron_from_homog_cone_input(BarvinokParameters *params);
  Polyhedron *read_polyhedron_hairy(BarvinokParameters *params);
public:
  bool expect_dilation_factor;
  bool expect_filename;
};

#endif
