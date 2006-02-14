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
      tmp->latticePoints=pointsInParallelepipedOfUnimodularCone(tmp->vertex,tmp->rays,numOfVars);     
      tmp=tmp->rest;
    }
  
  
  int oldnumofvars;
  oldnumofvars = 0;
  
  cones = ProjectUp(cones, oldnumofvars, numOfVars, templistVec);
  numOfVars = oldnumofvars;
 
  createGeneratingFunctionAsMapleInput(fileName,cones,numOfVars);  

}
