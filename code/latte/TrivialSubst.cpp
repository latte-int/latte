/* TrivialSubst.cpp -- Perform a trivial monomial substitution

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

#include "TrivialSubst.h"

void
TrivialMonomialSubstitutionMapleOutput_Single(ostream &out,
					      listCone *cone, int numOfVars)

{
  out << "(" << cone->coefficient << ")";
  out << " * (";
  listVector *lattice_point;
  bool first = true;
  for (lattice_point = cone->latticePoints;
       lattice_point != NULL;
       lattice_point = lattice_point->rest) {
    if (!first) out << " + ";
    out << "t^(" << lattice_point->first[numOfVars - 1] << ")";
    first = false;
  }
  out << ") / (";
  listVector *ray;
  first = true;
  for (ray = cone->rays; ray != NULL; ray = ray->rest) {
    if (!first) out << " * ";
    if (ray->first[numOfVars - 1] == 0) {
      cerr << "Trivial monomial substitution: Singularity encountered." << endl;
      exit(1);
    }
    out << "(1 - t^(" << ray->first[numOfVars - 1] << "))";
    first = false;
  }
  out << ")" << endl;
}

void
TrivialMonomialSubstitutionMapleOutput(ostream &out,
				       listCone *cones, int numOfVars)

{
  listCone *cone;
  bool first = true;
  for (cone = cones; cone!=NULL; cone=cone->rest) {
    if (!first) out << " + ";
    TrivialMonomialSubstitutionMapleOutput_Single(out, cone, numOfVars);
    first = false;
  }
}
  
/**************** CLASS IMPLEMENTATION ***************/

  
TrivialSubstitutionWritingConeConsumer::TrivialSubstitutionWritingConeConsumer(const string &filename)
  : genfun_stream(filename.c_str()), first_term(true)
{}

int TrivialSubstitutionWritingConeConsumer::ConsumeCone(listCone *cone) {
  if (cone->latticePoints != NULL) {
    if (!first_term)
      genfun_stream << " + ";

    TrivialMonomialSubstitutionMapleOutput_Single(genfun_stream, cone,
                                          cone->latticePoints->first.length());
    genfun_stream << endl;
    first_term = false;
  }
  freeCone(cone);
}
