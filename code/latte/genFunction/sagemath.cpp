/* 
   This modified version of /code/latte/genFunction/maple.cpp
   encodes short rational functions in a human- and machine-friendly format.
*/

/* sagemath.cpp -- Create Sagemath input

   Copyright 2002-2004 Jesus A. De Loera, David Haws, Raymond
      Hemmecke, Peter Huggins, Jeremy Tauzer, Ruriko Yoshida
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

#include <climits>
#include "sagemath.h"

/* ----------------------------------------------------------------- */
void writeTermToFile(ofstream & out, const vec_ZZ &v, int numOfVars) {
  int i;
  
  for (i = 0; i < numOfVars; ++i) {
    out << v[i];
    if (i < numOfVars - 1)
      out << " ";
  }
  out << endl;
}

/* ----------------------------------------------------------------- */
void writeTermOfGeneratingFunctionToFile(ofstream & out, listCone *cone,
					 int numOfVars) {
  listVector *tmp;

  if (cone->coefficient == 0)
    return;

  out << "{" << endl;

  out << cone->coefficient << endl;

  tmp = cone->latticePoints;
  out << lengthListVector(tmp) << endl;
  while (tmp) {
    writeTermToFile(out, tmp->first, numOfVars);
    tmp = tmp->rest;
  }

  tmp = cone->rays;
  out << lengthListVector(tmp) << endl;
  while (tmp) {
    writeTermToFile(out, tmp->first, numOfVars);
    tmp = tmp->rest;
  }

  out << "}" << endl;
}

/* ----------------------------------------------------------------- */
void createGeneratingFunctionAsSageMathInput(const char *fileName, listCone* cones, 
					  int numOfVars) {
  char sagemathInFileName[PATH_MAX];
  listCone *tmp;

  strcpy(sagemathInFileName,fileName);
  strcat(sagemathInFileName,".sage");

  ofstream out(sagemathInFileName);
  if (!out) {
    printf("Error opening output file in createGeneratingFunctionAsSageMathInput!");
    exit(1);
  }

  out << "gF:=";
  tmp=cones;
  while (tmp->rest) {
    writeTermOfGeneratingFunctionToFile(out,tmp,numOfVars);
    out << "+";
    tmp=tmp->rest;
  }
  writeTermOfGeneratingFunctionToFile(out,tmp,numOfVars);
  out << ";\n";

  out.close();
  return;
}
/* ----------------------------------------------------------------- */

void createGeneratingFunctionAsSageMathInputGrob(listCone* cones, 
					  int numOfVars, ofstream & out) {

  listCone *tmp;

  if (!out) {
    printf("Error opening output file in createGeneratingFunctionAsSageMathInput!");
    exit(1);
  }


  tmp=cones;
  while (tmp->rest) {
    writeTermOfGeneratingFunctionToFile(out,tmp,numOfVars);
    out << "+";
    tmp=tmp->rest;
  }
  writeTermOfGeneratingFunctionToFile(out,tmp,numOfVars);
   out << "+";


  return;
}
/* ----------------------------------------------------------------- */

GeneratingFunctionWritingConeConsumer::GeneratingFunctionWritingConeConsumer(const std::string &genfun_filename)
  : genfun_stream(genfun_filename.c_str()),
    first_term(true)
{
}


int GeneratingFunctionWritingConeConsumer::ConsumeCone(listCone *cone)
{
  if (cone->latticePoints != NULL) {
    writeTermOfGeneratingFunctionToFile(genfun_stream, cone,
					cone->latticePoints->first.length());
    genfun_stream << flush;
  }
  freeCone(cone);
  return 1; // means "success, please continue"
}
