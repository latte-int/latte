// This is a -*- C++ -*- header file.

/* VertexConesWith4ti2.h -- Compute vertex cones with 4ti2
	       
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

#ifndef VERTEXCONESWITH4TI2_H
#define VERTEXCONESWITH4TI2_H

#include "cone.h"
#include "cone_consumer.h"

/* Let MATRIX be a list of the vectors of an irredundant inequality
   system describing a full-dimensional polytope in dimension
   NUMOFVARS.  The vectors of MATRIX are of dimension NUMOFVARS + 1,
   the 0-th coordinate being the RHS.

   Compute the vertex cones of the polytope.
*/
listCone *
computeVertexConesWith4ti2(listVector* matrix, int numOfVars);

void
computeVertexConesWith4ti2(listVector* matrix, int numOfVars,
			   ConeConsumer &consumer);

#endif
