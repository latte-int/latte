// This is a -*- C++ -*- header file.

/* ConeDecom.cpp -- Barvinok's decomposition of a cone.

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

#ifndef CONEDECOM__H
#define CONEDECOM__H

#include "myheader.h"
#include "barvinok.h"

/* Do a signed decomposition, modulo lower-dimensional cones, of the
   polyhedral cone spanned by the ROW VECTORS of B with apex at
   CONE_VERTEX, until the determinants of all cones are at most
   PARAMETERS->max_determinant.  (When the cone is not simplicial, we
   triangulate it first.)

   Call PARAMETERS->ConsumeCone() for each of the small cones.
*/ 
int
barvinokDecomposition_Single (const mat_ZZ &B, rationalVector *Cone_Vertex,
			      Single_Cone_Parameters *Parameters);

/* Likewise, but the cone is given as a listCone (whose "rest" we
   ignore). */ 
int
barvinokDecomposition_Single (listCone *cone,
			      Single_Cone_Parameters *Parameters);

#endif






