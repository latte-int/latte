/* ProjectUp.cpp -- Compute the injection of a cone into a higher-dimensional space
	       
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

#include "ProjectUp.h"

/* ----------------------------------------------------------------- */
listCone* ProjectUp(listCone* cone, int & oldNumOfVars, int & newNumOfVars, 
             listVector *equations){

  listCone *current_cone = cone;
  vec_ZZ newVector;

  newVector.SetLength(oldNumOfVars);
 
  listVector *temp, *temp2, *current_ray, *new_ray;
  int i;

  while(current_cone)
  {

    temp2 = equations;
    //  cout << " Here 1" << endl;
    i = 0;
    while(temp2)
      {	  
	newVector[i] = temp2->first * current_cone->latticePoints->first;
	temp2 = temp2->rest;
	i++;
      }
    //  cout << " Here 2" << endl;
    for(i = oldNumOfVars - newNumOfVars; i < oldNumOfVars; i++)
      {
	newVector[i] = current_cone->latticePoints->first[i - oldNumOfVars + newNumOfVars];
      }
    // cout << " Here 3" << endl;
    delete current_cone->latticePoints;
    current_cone->latticePoints = new listVector;
    current_cone->latticePoints->rest = NULL;

    current_cone->latticePoints->first.SetLength(oldNumOfVars);
    //  cout << " Here 4" << endl;
    for(i = 0; i < oldNumOfVars; i++)
      current_cone->latticePoints->first[i] = newVector[i];

    current_ray = current_cone->rays;
    new_ray = new listVector;
    current_cone->rays = new_ray;
    //   cout << " Here 5" << endl;
    while(current_ray)
      {
	temp2 = equations;

	i = 0;
	while(temp2)
	  {	  
	    newVector[i] = temp2->first * current_ray->first;
	    temp2 = temp2->rest;
	    i++;
	  }

	for(i = oldNumOfVars - newNumOfVars; i < oldNumOfVars; i++)
	  {
	    newVector[i] = current_ray->first[i - oldNumOfVars + newNumOfVars];
	  }
	
	temp = current_ray;
	current_ray = current_ray->rest;
	delete temp;
	
	new_ray->first.SetLength(oldNumOfVars);
	//   cout << " Here 6" << endl;
	for(i = 0; i < oldNumOfVars; i++)
	  new_ray->first[i] = newVector[i];
      

	if(current_ray != NULL)
	  {
	    new_ray->rest = new listVector;
	    new_ray = new_ray->rest;
	  }
	else
	  new_ray->rest = NULL;
      }

    current_cone = current_cone->rest;
  }
  return cone;
}

/* ----------------------------------------------------------------- */
listCone* ProjectUp2(listCone* cone, int & oldNumOfVars, int & newNumOfVars, 
             mat_ZZ AA, vec_ZZ b){

  // d =  oldNumOfVars and k = newNumOfVars
  
  listCone *current_cone = cone;
  vec_ZZ newVector;
  
  newVector.SetLength(oldNumOfVars);
  
  listVector *temp, *current_ray, *new_ray;
  int i;
  
  while(current_cone)
    {
      
      //  cout << " Here 1" << endl;
      i = 0;
      newVector = b;
      
      for(i = 0; i < oldNumOfVars; i++){	  
	newVector[i] += AA[i] * current_cone->latticePoints->first;
      }
      
      //  cout << " Here 2" << endl;
      /*    for(i = oldNumOfVars - newNumOfVars; i < oldNumOfVars; i++)
	    {
	    newVector[i] = current_cone->latticePoints->first[i - oldNumOfVars + newNumOfVars];
	    }*/
	// cout << " Here 3" << endl;
      delete current_cone->latticePoints;
      current_cone->latticePoints = new listVector;
      current_cone->latticePoints->rest = NULL;
      
      current_cone->latticePoints->first.SetLength(oldNumOfVars);
      //  cout << " Here 4" << endl;
      for(i = 0; i < oldNumOfVars; i++)
	current_cone->latticePoints->first[i] = newVector[i];
      
      current_ray = current_cone->rays;
      new_ray = new listVector;
      current_cone->rays = new_ray;
      //   cout << " Here 5" << endl;
      while(current_ray)
	{
	  i = 0;
	  for(i = 0; i < oldNumOfVars; i++)
	    {	  
	      newVector[i] = AA[i] * current_ray->first;
	    }
	  
	  // 	for(i = oldNumOfVars - newNumOfVars; i < oldNumOfVars; i++)
	  // 	  {
	  // 	    newVector[i] = current_ray->first[i - oldNumOfVars + newNumOfVars];
	  // 	  }
	  
	  temp = current_ray;
	  current_ray = current_ray->rest;
	  delete temp;
	  
	  new_ray->first.SetLength(oldNumOfVars);
	  //   cout << " Here 6" << endl;
	  for(i = 0; i < oldNumOfVars; i++)
	    new_ray->first[i] = newVector[i];
	  
	  
	  if(current_ray != NULL)
	    {
	      new_ray->rest = new listVector;
	      new_ray = new_ray->rest;
	    }
	  else
	    new_ray->rest = NULL;
	}
      
      current_cone = current_cone->rest;
    }
  return cone;
}
