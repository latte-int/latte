// This is a -*- C++ -*- header file.

/* RegularTriangulationWith4ti2.h -- Triangulate with 4ti2
	       
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

#ifndef REGULARTRIANGULATIONWITH4TI2_H
#define REGULARTRIANGULATIONWITH4TI2_H

#include "barvinok/barvinok.h"
#include "triangulation/RegularTriangulation.h"

/* General function to compute a regular triangulation. */
void
triangulate_cone_with_4ti2(listCone *cone,
			   BarvinokParameters *Parameters,
			   height_function_type height_function,
			   void *height_function_data,
			   ConeConsumer &consumer);

/* Triangulate a full-dimensional CONE. */
void
random_regular_triangulation_with_4ti2(listCone *cone,
				       BarvinokParameters *Parameters,
				       ConeConsumer &consumer);

void
refined_delone_triangulation_with_4ti2(listCone *cone,
				       BarvinokParameters *Parameters,
				       ConeConsumer &consumer);

#endif
