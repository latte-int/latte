// This is a -*- C++ -*- header file.

/* Residue.cpp -- Polynomial substitution and residue calculations

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

#ifndef RESIDUE_H
#define RESIDUE_H

ZZ Residue(listCone* cones, int numOfVars);

// Returns a -1 if dot product in denominator is 0, 1 otherwise
// In the latter case, consumes CONES.
int Residue_Single_Cone (listCone *cones, int numOfVars,
			 const vec_ZZ &Random_Lambda,
			 ZZ *Total_Lattice_Points,
			 ZZ *Ten_Power);

#endif
