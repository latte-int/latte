/* maple.cpp -- Create Maple input

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

#include "../cone.h"
#include "../ramon.h"
#include "../vertices/cdd.h"
/* ----------------------------------------------------------------- */
void writeTermToFile(ofstream & out, vec_ZZ v, int numOfVars) {
  int i,firstEntry;

  firstEntry=0;

  for (i=0;i<numOfVars;i++) {
    if (v[i]!=0) {
      if (firstEntry==1) out << "*";
      firstEntry=1;
      if (v[i]<0) out << "x[" << i << "]^(" << v[i] << ")";
      if (v[i]==1) out << "x[" << i << "]";
      if (v[i]>1) out << "x[" << i << "]^" << v[i];
    }
  }
  /* Is v=0? */
  if (firstEntry==0) {out << "1";}
 
  return ;
}
/* ----------------------------------------------------------------- */
void writeTermOfGeneratingFunctionToFile(ofstream & out, listCone *cone, 
					 int numOfVars) {
  int len;
  vec_ZZ v;
  listVector *tmp;

  if (cone->coefficient==0) return;
  if (cone->coefficient!=1) out << "(" << cone->coefficient << ")*";

  tmp=cone->latticePoints;
  len=lengthListVector(tmp);
  if (len>1) out << "(";
  while (tmp) {
    v=tmp->first;
    writeTermToFile(out,v,numOfVars);
    if (tmp->rest!=0) out << "+";
    tmp=tmp->rest;  
  }
  if (len>1) out << ")";

  out << "/";

  tmp=cone->rays;
  out << "(";
  while (tmp) {
    out << "(1-";
    v=tmp->first;
    writeTermToFile(out,v,numOfVars);
    out << ")";
    if (tmp->rest!=0) out << "*";
    tmp=tmp->rest;  
  }
  out << ")";

  return;
}
/* ----------------------------------------------------------------- */
void createGeneratingFunctionAsMapleInput(char *fileName, listCone* cones, 
					  int numOfVars) {
  char mapleInFileName[127];
  listCone *tmp;

  strcpy(mapleInFileName,fileName);
  strcat(mapleInFileName,".maple");

  ofstream out(mapleInFileName);
  if (!out) {
    printf("Error opening output file in createGeneratingFunctionAsMapleInput!");
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

void createGeneratingFunctionAsMapleInputGrob(listCone* cones, 
					  int numOfVars, ofstream & out) {

  listCone *tmp;

  if (!out) {
    printf("Error opening output file in createGeneratingFunctionAsMapleInput!");
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

