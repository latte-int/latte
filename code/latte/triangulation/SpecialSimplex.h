// This is a -*- C++ -*- header file.

/* SpecialSimplex.h -- Check for a special simplex using CPLEX
	       
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

#ifndef SPECIALSIMPLEX_H
#define SPECIALSIMPLEX_H

#include "cone.h"
#include "cone_consumer.h"
#include "barvinok/barvinok.h"
#include "triangulation/RegularTriangulation.h"

/* Find a `special' full-dimensional simplicial cone spanned by some
   of the rays of CONE, i.e., one that contains +/- the n-th unit
   vector in its interior. */
listCone *
FindSpecialSimplex(listCone *cone, int numOfVars);

/* A height function for a regular triangulation that guarantees to
   make SPECIAL_CONE (DATA) a cone of the trinagulation. */
void
special_height(mpq_t height, const vec_ZZ &ray, void *data);

/* Return a triangulation of a full-dimensional CONE, such that the
   resulting facets avoid the subspace where the last coordinate is 0.
   This fails when the cone does not have a special simplicial cone.
*/
void
special_triangulation_with_subspace_avoiding_facets
(listCone *cone, BarvinokParameters *Parameters, ConeConsumer &consumer);

#endif
