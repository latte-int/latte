// This is a -*- C++ -*- header file.

/* convert.cpp -- Data conversions

   Copyright 2002-2004 Jesus A. De Loera, David Haws, Raymond
      Hemmecke, Peter Huggins, Jeremy Tauzer, Ruriko Yoshida
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

#ifndef CONVERT_H
#define CONVERT_H

#include "cone.h"

listVector *
transformArrayBigVectorToListVector(const mat_ZZ &A,
				    int numOfVectors,
				    int numOfVars);

/* Create a matrix whose ROWS are the ray vectors of CONE. */
mat_ZZ
createConeDecMatrix(const listCone *cone, int numOfRays, int numOfVars);

/* Create a matrix whose ROWS are the facet vectors of CONE,
   scaled in a way such that
   
   < RAY_i, SCALED_FACET_j > = -DET(RAYS) * DELTA_{i,j}.
*/
mat_ZZ
createFacetMatrix(const listCone *cone, int numOfFacets, int numOfVars);


/* Likewise, but
   < RAY_i, SCALED_FACET_j > = - |DET(RAYS)| * DELTA_{i,j}.
*/
mat_ZZ
createFacetMatrix2(const listCone *cone, int numOfFacets, int numOfVars);

/* latte to NTL conversions */
mat_ZZ
convert_listVector_to_mat_ZZ(listVector *);

#endif
