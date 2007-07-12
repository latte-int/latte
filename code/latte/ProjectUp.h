// This is a -*- C++ -*- header file.

/* ProjectUp.h -- Compute the injection of a cone into a higher-dimensional space

   Extracted from ReadingFile.h, which is:
     Copyright 2002, 2003 Raymond Hemmecke, Ruriko Yoshida
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

#ifndef PROJECTUP_H
#define PROJECTUP_H

#include "cone.h"
#include "cone_consumer.h"

listCone* ProjectUp(listCone* cone, int & oldNumOfVars, int & newNumOfVars, 
		    listVector *equations);
listCone* ProjectUp2(listCone* cone, int & oldNumOfVars, int & newNumOfVars, 
		     mat_ZZ AA, vec_ZZ b);

class ProjectingUpConeTransducer : public ConeTransducer
{
public:
  ProjectingUpConeTransducer(int a_oldNumOfVars, int a_newNumOfVars,
			     const mat_ZZ &a_AA, const vec_ZZ &a_b)
    : oldNumOfVars(a_oldNumOfVars),
      newNumOfVars(a_newNumOfVars),
      AA(a_AA), b(a_b) {}
  int ConsumeCone(listCone *cone);
private:
  int oldNumOfVars;
  int newNumOfVars;
  mat_ZZ AA;
  vec_ZZ b;
};

#endif
