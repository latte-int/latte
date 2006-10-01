/* RudyResNTL.h -- Polynomial substitution and residue calculations

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

#ifndef RUDYRESNTL__H
#define RUDYRESNTL__H

#include "cone.h"
#include "barvinok/dec.h" // for Standard_Single_Cone_Parameters

void ResidueFunction(listCone* cones, int numOfVars, int print_flag, int degree, int output_cone);

// Returns -1 if a Dot Product is zero in the denominator, otherwise 1 if ok
int
ResidueFunction_Single_Cone (listCone *cones,
			     Standard_Single_Cone_Parameters *Residue_Parameters); 

#endif

