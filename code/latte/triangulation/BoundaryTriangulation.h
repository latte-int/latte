// This is a -*- C++ -*- header file.

/* BoundaryTriangulation.h -- Boundary triangulation
	       
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

#ifndef BOUNDARYTRIANGULATION_H
#define BOUNDARYTRIANGULATION_H

#include "cone.h"
#include "barvinok/barvinok.h"

/* Return a triangulation of a full-dimensional CONE, such that the
   resulting facets avoid the subspace where the last coordinate is 0. */
listCone *
boundary_triangulation_of_cone_with_subspace_avoiding_facets
(listCone *cone, BarvinokParameters *Parameters);

#endif