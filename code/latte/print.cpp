#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include <list>
#include <vector>

#include "myheader.h"
#include "cone.h"
#include "ramon.h"
#include "print.h"

/* ----------------------------------------------------------------- */
void printVector(vec_ZZ v, int numOfVars) {
  int i;

//    if (v==0) {
//      cout << "[]\n";
//      return ;
//    }
  cout << "[";
  for (i=0; i<(numOfVars-1); i++) {
    cout << v[i] << " ";
  }
  cout << v[i] << "]" << endl;
  return ;
}
/* ----------------------------------------------------------------- */
void printListVector(listVector* basis, int numOfVars) {
  if (basis==0) cout << "[]\n";
  while(basis) {
    printVector(basis->first,numOfVars);
    basis = basis->rest;
  }
/*  printf("\n"); */
  return ;
}
/* ----------------------------------------------------------------- */
void printRationalVector(rationalVector *v, int numOfVars) {
  int i;

//    if (v==0) {
//      cout << "[]\n";
//      return ;
//    }

  cout << "[";
  for (i=0; i<(numOfVars-1); i++) {
    if ((v->denominator)[i]==1)
      cout << (v->enumerator)[i] << " ";
    else
      cout << (v->enumerator)[i] << "/" << (v->denominator)[i] << " ";
  }

  if ((v->denominator)[i]==1)
    cout << (v->enumerator)[i] << "]" << endl;
  else
    cout << (v->enumerator)[i] << "/" << (v->denominator)[i] << "]" << endl;
  return ;
}
/* ----------------------------------------------------------------- */
void printCone(listCone* cones, int numOfVars) {
  printConeToFile(cout, cones, numOfVars);
}
/* ----------------------------------------------------------------- */
void printListCone(listCone* cones, int numOfVars) {
  if (cones==0) cout << "No cones in list.\n";
  while(cones) {
    printCone(cones,numOfVars);
    cones = cones->rest;
  }
  cout << endl;
  return ;
}
/* ----------------------------------------------------------------- */
void printVectorToFile(ostream & out, vec_ZZ v, int numOfVars) {
  int i;

//    if (v==0) {
//      out << "[]\n";
//      return ;
//    }
  out << "[";
  for (i=0; i<(numOfVars-1); i++) {
    out << v[i] << " ";
  }
  out << v[i] << "]\n";
  return ;
}
/* ----------------------------------------------------------------- */
void printListVectorToFile(ostream & out, listVector* basis, int numOfVars) {
  if (basis==0) {
    out << "[]\n";
    return;
  }

  while(basis) {
    printVectorToFile(out,basis->first,numOfVars);
    basis = basis->rest;
  }
  return ;
}
/* ----------------------------------------------------------------- */
void printVectorToFileWithoutBrackets(ostream & out, vec_ZZ v, 
				      int numOfVars) {
  int i;

//    if (v==0) return ;

  for (i=0; i<(numOfVars-1); i++) {
    out << v[i] << " ";
  }
  out << v[i] << endl;
  return ;
}
/* ----------------------------------------------------------------- */
void printListVectorToFileWithoutBrackets(ostream & out, listVector* basis, 
					  int numOfVars) {
  if (basis==0) {
    out << numOfVars << " 0\n";
    return;
  }

  while(basis) {
    printVectorToFileWithoutBrackets(out,basis->first,numOfVars);
    basis = basis->rest;
  }
  return ;
}
/* ----------------------------------------------------------------- */
void printRationalVectorToFile(ostream & out, rationalVector *v, 
			       int numOfVars) {
  int i;

  if (v==0) {
    out << "[]\n";
    return ;
  }
  out << "[";
  for (i=0; i<(numOfVars-1); i++) {
    if ((v->denominator)[i]==1)
      out << (v->enumerator)[i] << " ";
    else
      out << (v->enumerator)[i] << "/" << (v->denominator)[i] << " ";
  }
  
  if ((v->denominator)[i]==1)
    out << (v->enumerator)[i] << "]\n";
  else
    out << (v->enumerator)[i] << "/" << (v->denominator)[i] << "]\n";
  return ;
}
/* ----------------------------------------------------------------- */
void printRationalVectorToFileWithoutBrackets(ostream & out, 
					      rationalVector *v, 
					      int numOfVars) {
  int i;

  if (v==0) {
    return ;
  }
  for (i=0; i<(numOfVars); i++) {
    if ((v->denominator)[i]==1)
      out << (v->enumerator)[i] << " "; 
    else
      out << (v->enumerator)[i] << "/" << (v->denominator)[i] << " ";
  }

  out << endl;
  return ;
}
/* ----------------------------------------------------------------- */
void printConeToFile(ostream & out,listCone* cones, int numOfVars) {
  out << "==========\n";
  out << "Cone.\n";

  out << "Coefficient: " << cones->coefficient << endl;

  out << "Vertex: ";
  printRationalVectorToFile(out,cones->vertex,numOfVars);

  out << "Extreme rays:\n";  
  printListVectorToFile(out,cones->rays,numOfVars);

  cout << "Determinant:" << cones->determinant << endl;
  
  cout << "Facets:\n";  
  printListVector(cones->facets,numOfVars);

  out << "Lattice points in parallelepiped:\n";
  printListVectorToFile(out,cones->latticePoints,numOfVars);
  out << "==========\n\n";

  return ;
}
/* ----------------------------------------------------------------- */
void printListConeToFile(char *fileName, listCone* cones, int numOfVars) {
  ofstream out(fileName);
  if (!out) {
    printf("Error opening output file for writing in printListConeToFile!");
    exit(0);
  }

  if (cones==0) out << "No cones in list.\n";

  while (cones) {
    printConeToFile(out,cones,numOfVars);
    cones = cones->rest;
  }
  out << endl;

  out.close();
  return ;
}
/* ----------------------------------------------------------------- */
void printResidueFile(char* fileName, listCone* cones, int numOfVars) {
  int numOfTerms;
  char outFileName[127];
  listVector *tmp;
  listCone *C;

  strcpy(outFileName,fileName);
  strcat(outFileName,".residue");

  ofstream out(outFileName);
  if (!out) {
    printf("Error opening output file for writing in printResidueFile!");
    exit(0);
  }
  if (cones==0) out << "No cones in list.\n";

  numOfTerms=0;

  C=cones;
  while (C) {
    numOfTerms=numOfTerms+lengthListVector(C->latticePoints);
    C=C->rest;
  }


  out << numOfVars << " " << lengthListVector(cones->rays) << " " <<
    numOfTerms << "\n\n";

  while (cones) {
    tmp=cones->latticePoints;
    while (tmp) {
      out << cones->coefficient << endl;
      printVectorToFileWithoutBrackets(out,tmp->first,numOfVars);
      printListVectorToFileWithoutBrackets(out,cones->rays,numOfVars);
      out << endl;
      tmp=tmp->rest;
    }
    cones = cones->rest;
  }
  out << endl;

  out.close();
  return ;
}

