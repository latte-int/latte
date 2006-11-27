/* ConeDecom.cpp -- Barvinok's decomposition of a cone.

   Copyright 2002,2003 Ruriko Yoshida
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

#include <list>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <string>

#include "../ramon.h"
#include "../print.h"
#include "../cone.h"
#include "convert.h"
#include "barvinok.h"
#include "triangulation/triangulate.h"

/*
  The first step is to triangulate a cone into simplicial cones.
  Then, by using Barvinok's decomposition, we decompose each
  simplicial cone into unimodular cones.
*/

int
barvinokDecomposition_Single (listCone *cone,
			      Single_Cone_Parameters *Parameters)
{
  int status = 1;
  listCone *triang = triangulateCone(cone, Parameters->Number_of_Variables, Parameters);
  listCone *t;
  for (t = triang; t!=NULL; t=t->rest) {
    int num_rays = lengthListVector(t->rays);
    mat_ZZ B = createConeDecMatrix(t, num_rays, Parameters->Number_of_Variables);
    if ((status = barvinok_Single(B, Parameters, t->vertex)) == -1)
      goto BAILOUT;
  }
 BAILOUT:
  freeListCone(triang);
  return status;
}

