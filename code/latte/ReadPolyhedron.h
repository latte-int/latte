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

#include <cassert>
#include "Polyhedron.h"
#include "latte_cddlib.h"
#include "barvinok/barvinok.h"
#include "LattException.h"

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
  bool input_homog_cone, input_vertex_cones, input_dualized, have_subcones, input_listcone_format;
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
    RedundancyCheckWithCddlib,
    NoRedundancyCheck,
    FullRedundancyCheckWithCddlib
  } RedundancyCheckType;
  RedundancyCheckType redundancycheck;

  //this enum is to be used by the user in describing the output of Polyhedron.
  typedef enum {
	  computePrimalCones,
	  computeVertices,
  } ReadPolyhedronOutput;
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
  Polyhedron *read_polyhedron(dd_MatrixPtr M, BarvinokParameters *params, const ReadPolyhedronOutput readPolyhedronOutput);
  listVector * read_full_rank_inequality_matrix(BarvinokParameters *params); //Returns the system Ax <= b where the polytope is full dimensional, starting from the file name.
protected:
  void matrixToVerticesOrCones(listVector * theMatrix, int numOfVars, Polyhedron *& Poly, BarvinokParameters *&params);
  Polyhedron *read_polyhedron_from_homog_cone_input(BarvinokParameters *params);
  Polyhedron *read_polyhedron_from_vertex_cone_input(BarvinokParameters *params);
  Polyhedron *read_polyhedron_hairy(BarvinokParameters *params);
  Polyhedron *PolyhedronFromHrepMatrix(dd_MatrixPtr M, BarvinokParameters *params);
  void polyhedronRedundancyCheck(RedundancyCheckType redunType, dd_MatrixPtr &M);
  mat_ZZ findLatticeBasis(dd_MatrixPtr &M, int &numOfVars);
  listVector * projectOutVariables(dd_MatrixPtr &M, int &numOfVars, Polyhedron *& Poly); //reduces the input matrix to a full-dimensional matrix.
public:
  bool expect_dilation_factor;
  bool expect_filename;
};


/**
 * Class ReadPolyhedronDataRecursive's job is to set each facet inequalities to equality
 *   and compute the reduced polytope.
 *
 *  DO NOT USE THIS CLASS: it is experimental/partly finished.
 */
class ReadPolyhedronDataRecursive: public ReadPolyhedronData
{
private:
	dd_MatrixPtr ddHrep;
	mat_ZZ latticeBasis;
	mat_ZZ latticeLeftInverse;
	ZZ     latticeLeftInverseDilation;
	ZZ dilationNum;
public:
	ReadPolyhedronDataRecursive(const ReadPolyhedronData & rpd); //consturctor.
	~ReadPolyhedronDataRecursive(); //we are in charge of freeing the memory of matrix.

	//A B C D E F G H I J K L M N O P Q R S T U V W X Y Z
	void dilatePolytope();
	int dimension();
	int getFullDimensionCount() const;

	Polyhedron * findTangentCones();
	void getFacetPolytope(int row, ReadPolyhedronDataRecursive &newMatrix, vec_ZZ & l, RationalNTL &lDotNormal);
	const mat_ZZ * getLatticeInverse() const;
	const ZZ * getLatticeInverseDilation() const;
	RationalNTL getNormalFactor() const;
	int getNumberEqualities() const;
	int getNumberRows() const;
	void latticeInverse();
	void readHrepMatrixFromFile(string fileName, BarvinokParameters *params);
	void readHrepMatrix();

	//void setInequalityMatrix(listVector *newMatrix);

	void setInequalityToEquality(int i, listVector * &newMatrix, BarvinokParameters &newParm);
	void setInequality(int row);
	void setMatrix(dd_MatrixPtr m);

	RationalNTL volumeCorrection(const RationalNTL & a) const;

};





/* Helper functions. */

/* Read a VREP file in LattE format
   and create a corresponding Polyhedron. */
Polyhedron *ReadLatteStyleVrep(const char *filename, bool homogenize); //not sure if this function exist.

/* Create a polyhedron from a vrep matrix. */
Polyhedron *PolyhedronFromVrepMatrix(dd_MatrixPtr matrix, bool homogenize);

static dd_MatrixPtr ReadCddStyleMatrix(const string &filename);
#endif
