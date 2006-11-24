// This is a -*- C++ -*- header file.

/* TrivialSubst.h -- Perform a trivial monomial substitution

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

#ifndef TRIVIALSUBST_H
#define TRIVIALSUBST_H

#include <iostream>
#include "cone.h"

// We substitute (x_1,\dots,x_n) -> x_n.
void
TrivialMonomialSubstitutionMapleOutput_Single(ostream &out,
					      listCone *cone, int numOfVars);  

void
TrivialMonomialSubstitutionMapleOutput(ostream &out,
				       listCone *cones, int numOfVars);

#endif
