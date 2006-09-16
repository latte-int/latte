/* tmpGro.cpp -- 

   Copyright 2002-2004 Jesus A. De Loera, David Haws, Raymond
      Hemmecke, Peter Huggins, Jeremy Tauzer, Ruriko Yoshida

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

void Grobner(){
  
  
  readLatteProblem(fileName,&equations,&inequalities,equationsPresent,
		   &numOfVars, nonneg, dualApproach, grobner);
  
  numOfVars--;
  numOfAllVars=numOfVars;
  
  generators=createArrayVector(numOfVars);
  matrix=Grobner(equations,inequalities,&generators,&numOfVars, &templistVec);
  
  cones=computeVertexCones(fileName,matrix,numOfVars);
  
  cones=dualizeCones(cones,numOfVars);
  
  cones=decomposeCones(cones,numOfVars, flags, fileName);
    
  cones=dualizeBackCones(cones,numOfVars);
  
    
  tmp=cones;
  while (tmp) 
    { 
      tmp->latticePoints=pointsInParallelepiped(tmp,tmp,numOfVars);     
      tmp=tmp->rest;
    }
  
  
  int oldnumofvars;
  oldnumofvars = 0;
  
  cones = ProjectUp(cones, oldnumofvars, numOfVars, templistVec);
  numOfVars = oldnumofvars;
 
  createGeneratingFunctionAsMapleInput(fileName,cones,numOfVars);  

}
