// This is a -*- C++ -*- header file.

/* dual.h -- Dualize polyhedral cones

   Copyright 2002 Raymond Hemmecke, Ruriko Yoshida
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

#ifndef DUAL_H
#define DUAL_H

#include "cone.h"
#include "barvinok/barvinok.h"

BarvinokParameters::DualizationType
dualization_type_from_name(const char *name);

bool
parse_standard_dualization_option(const char *arg,
				  BarvinokParameters *params);

/* Dualize the polyhedral CONE. */
void dualizeCone(listCone *cone, int numOfVars,
		 BarvinokParameters *params);

/* Dualize the polyhedral cones in the list. */
listCone* dualizeCones(listCone *cones, int numOfVars,
		       BarvinokParameters *params);

/* Fill the slots `determinant', `facets', and `facet_divisors' of
   CONE.  The facet vectors are made primitive.  
*/
void computeDetAndFacetsOfSimplicialCone(listCone *cone, int numOfVars);

/* Destructively dualize the simplicial cones in the list.
   The ray vectors of the resulting cones are made primitive.
*/
listCone* dualizeBackCones(listCone*, int);

/* Fill the slot `facet' of all CONES by determining which of the
   INEQUALITIES are tight at the respective vertex. */
void computeTightInequalitiesOfCones(listCone *cones,
				     listVector *inequalities,
				     int numOfVars);

#endif
