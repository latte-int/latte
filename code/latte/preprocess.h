/* preprocess.h -- Preprocessing and projecting polytopes

   Copyright 2002 Raymond Hemmecke, Ruriko Yoshida

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

#ifndef PREPROCESS__H
#define PREPROCESS__H
int ihermite(vec_ZZ *S, vec_ZZ *U, vec_ZZ* rhs, int m, int n);
listVector* preprocessProblem(listVector*, listVector*, vec_ZZ**, int*, vec_ZZ&, mat_ZZ &, char*, int);
 listVector* TransformToDualCone(listVector* matrix, int& numOfVars);
void dilateListVector(listVector* basis, int numOfVars, int dil);
vec_ZZ transpose(vec_ZZ mat, int numOfVars, int numOfRows);
vec_ZZ ProjectingUp(mat_ZZ ProjU, vec_ZZ cost, int numOfVars);
vec_RR ProjectingUpRR(mat_RR ProjU, vec_RR cost, int numOfVars);
void Interior(listVector* basis);
#endif
