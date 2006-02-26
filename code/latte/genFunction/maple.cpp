#include "../myheader.h"
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
    exit(0);
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
    exit(0);
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

