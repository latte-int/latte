// This is a -*- C++ -*- header file.

/* TriangulationWithTOPCOM.h -- Compute a triangulation using TOPCOM
	       
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

#ifndef TRIANGULATIONWITHTOPCOM_H
#define TRIANGULATIONWITHTOPCOM_H

#include "cone.h"
#include <list>

/* Return a triangulation of a full-dimensional CONE.
   CONE is consumed. */
listCone *
triangulate_cone_with_TOPCOM(listCone *cone, int numOfVars);

/* Return all triangulations. */
std::list<listCone *>
all_triangulations_of_cone_with_TOPCOM(listCone *cone, int numOfVars);

#endif
