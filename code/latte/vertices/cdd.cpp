/* cdd.cpp -- Computation of all vertex cones via CDD

   Copyright 2002 Raymond Hemmecke, Ruriko Yoshida
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

#include "config.h"
#include "cone.h"
#include "print.h"
#include "ramon.h"
#include "rational.h"
#include <list>
#include <cassert>
#include "latte_system.h"

using namespace std;
/* ----------------------------------------------------------------- */


listVector* CopyListVector(listVector* A, int numOfVars) {

  vec_ZZ v;
  listVector *L, *endL;

  v=createVector(numOfVars);
  L=createListVector(v);
  endL=L;

  while(A) {
    v = A -> first;
    endL->rest = createListVector(v);
    endL = endL->rest;
    A = A->rest;
  }

  return (L->rest);
}

/* ----------------------------------------------------------------- */
listCone* CopyListCones(listCone* RudyCones, int numOfVars, 
			rationalVector* Opt_vertex) {
  int s;
  listCone *cones, *endCones, *newCone, *tmp;

  cones=createListCone();
  endCones=cones;
  tmp = RudyCones;
  while (tmp) {

    if((tmp->vertex->vertex->numerators()==Opt_vertex->numerators())
       && (tmp->vertex->vertex->denominators()==Opt_vertex->denominators())){ 
      newCone=createListCone();
      s = tmp->coefficient;
      newCone->coefficient=s;
      newCone->rays = CopyListVector(tmp->rays, numOfVars);
      newCone->facets = CopyListVector(tmp->facets, numOfVars);
      newCone->vertex = tmp->vertex;
      endCones->rest=newCone;
      endCones=endCones->rest;
    }
    tmp = tmp ->rest; 
  }
  //  printListCone(cones->rest,numOfVars);
  return (cones->rest);
}

/* ----------------------------------------------------------------- */
listCone* CopyListCones(listCone* RudyCones, int numOfVars) {
  int s;
  listCone *cones, *endCones, *newCone, *tmp;
  listVector *latticepoint;
  cones=createListCone();
  endCones=cones;
  tmp = RudyCones;
  while (tmp) {

//    if((tmp->vertex->numerators()==Opt_vertex->numerators()) && 
//(tmp->vertex->denominators()==Opt_vertex->denominators())){ 
      newCone=createListCone();
      latticepoint = createListVector(tmp->latticePoints->first);
      newCone->latticePoints = latticepoint;
      s = tmp->coefficient;
      newCone->coefficient=s;
      newCone->rays = CopyListVector(tmp->rays, numOfVars);
      newCone->facets = CopyListVector(tmp->facets, numOfVars);
      newCone->vertex = tmp->vertex;
      endCones->rest=newCone;
      endCones=endCones->rest;
    
    tmp = tmp ->rest; 
  }
  //  printListCone(cones->rest,numOfVars);
  return (cones->rest);
}

/* ----------------------------------------------------------------- */

void createCddIneFile(listVector* matrix, int numOfVars) {
  int i;
  listVector* tmp;

  ofstream OUT;
  OUT.open("latte_cdd.ine");
  OUT << "H-representation" << endl;
  OUT << "begin " << endl;
  OUT << lengthListVector(matrix) << " " << numOfVars << " integer" << endl;
  tmp=matrix;  
  while (tmp) {
    for (i=0; i<(numOfVars); i++) OUT << (tmp->first)[i] << " ";
    OUT << endl;
    tmp=tmp->rest;
  }
  OUT << "end" << endl;
  OUT << "adjacency" << endl;
  OUT << "incidence" << endl;
  OUT.close();

  return ;
}
/* ----------------------------------------------------------------- */
void createLrsIneFile(listVector* matrix, int numOfVars) {
  int i;
  listVector* tmp;

  ofstream OUT;
  OUT.open("latte_lrs.ine");
  OUT << "H-representation" << endl;
  OUT << "begin" << endl;
  OUT << lengthListVector(matrix) << " " << numOfVars << " integer" << endl;
  tmp=matrix;  
  while (tmp) {
    for (i=0; i<(numOfVars); i++) OUT << (tmp->first)[i] << " ";
    OUT << endl;
    tmp=tmp->rest;
  }
  OUT << "end" << endl;
  OUT.close();

  return ;
}

/* ----------------------------------------------------------------- */
void createLrsIneFileToPostAnalysys(listVector* matrix, int numOfVars) {

  ifstream in;
  ofstream OUT;
  OUT.open("latte_cdd.ine");
  in.open("latte_lrs.ine");
  int x, y;
  string tmpString;
  while (tmpString!="begin") getline(in,tmpString);
  OUT << "H-representation" << endl;
  OUT << "begin" << endl;
  in >> x >> y;
  getline(in,tmpString);
  OUT << x << " " << y << " integer" << endl;
  while (tmpString!="end") {getline(in,tmpString); OUT << tmpString << endl;}

  OUT << "postanalysis " << endl;
  OUT << "adjacency" << endl;

  OUT.close();

  return ;
}

/* ----------------------------------------------------------------- */
void createLrsExtFileToPostAnalysys(listVector* matrix, int numOfVars) {
  ifstream in, in2;
  ofstream OUT;
  OUT.open("latte_cdd.ext");
  in.open("latte_lrs.ext");
  in2.open("latte_lrs.ext");
  int counter = 0;
  string tmpString;
  while (tmpString!="begin") getline(in2,tmpString);
  getline(in2,tmpString);
  while (tmpString!="end") {getline(in2,tmpString); counter++;}
  counter--;
  while (tmpString!="begin") getline(in,tmpString);  
  OUT << "V-representation" << endl;
  OUT << "begin" << endl;
  getline(in,tmpString);
  OUT << counter << " " << numOfVars << " rational" << endl;
  while (tmpString!="end") {getline(in,tmpString); OUT << tmpString << endl;}

  OUT << "hull " << endl;
 
  OUT.close();

  return ;
}


/* ----------------------------------------------------------------- */


void createCddExtFile(listVector* matrix, int numOfVars) {
  int i;
  listVector* tmp;

  ofstream OUT;
  OUT.open("latte_cdd.ext1");
  OUT << "V-representation" << endl;
  OUT << "begin" << endl;
  OUT << lengthListVector(matrix) << " " << numOfVars << " integer" << endl;
  tmp=matrix;  
  while (tmp) {
    for (i=0; i<(numOfVars); i++) OUT << (tmp->first)[i] << " ";
    OUT << endl;
    tmp=tmp->rest;
  }
  OUT << "end" << endl;
  OUT << "hull" << endl;
  OUT.close();

  return ;
}

/* ----------------------------------------------------------------- */
void createCddExtFile2(const char* filename) {
  int i, numOfVec, numOfVars;
  string tmpString;
  // listVector* tmp;
  ifstream in(filename);
  in >> numOfVec >> numOfVars;
  ofstream OUT;
    getline(in,tmpString);
  OUT.open("latte_cdd.ext");
  OUT << "V-representation" << endl;
  OUT << "begin" << endl;
  OUT << numOfVec << " " << numOfVars << " integer" << endl;
  for(i = 0; i < numOfVec; i++){
    getline(in,tmpString);
    OUT << tmpString << endl;
  }
  OUT << "end" << endl;
  OUT << "hull" << endl;
  OUT.close();

  return ;
}
/* ----------------------------------------------------------------- */
void createCddIneLPFile(listVector* matrix, int numOfVars, vec_ZZ & cost) {
  int i;
  listVector* tmp;
  // cout << cost << " " << numOfVars << endl;
  ofstream OUT;

  OUT.open("LP.ine");
  OUT << "H-representation" << endl;
  OUT << "begin " << endl;
  OUT << lengthListVector(matrix) << " " << numOfVars << " integer" << endl;
  tmp=matrix;  
  while (tmp) {
    for (i=0; i<(numOfVars); i++) OUT << (tmp->first)[i] << " ";
    OUT << endl;
    tmp=tmp->rest;
  }
  OUT << "end" << endl;
  OUT << "maximize" << endl;
  OUT << 0 << " ";
  for(i = 0; i < numOfVars - 1; i++)
    OUT << cost[i] << " ";
  OUT << endl;
  OUT.close();

  return ;
}

/* ----------------------------------------------------------------- */

rationalVector* ReadLpsFile(int numOfVars, bool verbose = true)
{
  ifstream in("LP.lps");
  string tmpString;
  ZZ x, y;
  if (verbose) {
    cout << "Reading .lps file...";
    cout.flush();
  }
  rationalVector* OptVector;
  OptVector = createRationalVector(numOfVars);
  if(!in){
    cerr << "Cannot open input file in ReadLpsFile." << endl;
    exit(1);
  }

  while (tmpString!="begin") getline(in,tmpString);
  in >> tmpString;
  for(int i = 0; i < numOfVars; i++)
    {
      in >> tmpString >> tmpString;

      x=0;
      y=0;
      ReadCDD(in,x,y);
      OptVector->set_entry(i, x, y);
      
    }
  if (verbose) {
    cout <<"done." << endl;
  }
  return OptVector;
}

/* ----------------------------------------------------------------- */


listVector* createListOfInequalities(listVector* matrix, int numOfVars) {
  int i,j;
  ZZ g;
  vec_ZZ v;
  listVector *tmp, *inequalities, *endInequalities;

  /* Copying equality constraints as Az<=b, -Az<=-b. */

  tmp=matrix;  
  inequalities=createListVector(createVector(numOfVars));
  endInequalities=inequalities;

  while (tmp) {
    v=createVector(numOfVars);
    for (i=0; i<(numOfVars); i++) v[i]=-(tmp->first)[i+1];

/* Normalize vector by removing gcd of entries. */
    g=v[0];
    for(i=1; i<numOfVars; i++) g=GCD(g,v[i]);
    g=abs(g);
    if (g!=1) {
      for (i=0; i<numOfVars; i++) 
	v[i]=(v[i])/g;
    }

    endInequalities->rest=createListVector(v);
    endInequalities=endInequalities->rest;

/*    for (i=0; i<(numOfVars); i++) v[i]=(tmp->first)[i+1]);
    endInequalities->rest=createListVector(v);
    endInequalities=endInequalities->rest; */
    tmp=tmp->rest;
  }

  /* Writing non-negativity constraints. */
  for (i=0; i<(numOfVars); i++) {   
    v=createVector(numOfVars);
    for (j=0; j<(numOfVars); j++) {
      if (i==j) v[i]=1;
      else v[j]=0;
    }
    endInequalities->rest=createListVector(v);
    endInequalities=endInequalities->rest;
  }    

  return(inequalities->rest);
}
/* ----------------------------------------------------------------- */
listCone* readCddExtFile(int &numOfVars) {
  int i,j,numOfVertices;
  ZZ x,y;
  char cddInFileName[127];
  rationalVector *v;
  listCone *cones,*endCones,*c;
  string tmpString;

  cout << "Reading .ext file..."; cout.flush();

  strcpy(cddInFileName,"latte_cdd.ext");

  ifstream in(cddInFileName);
  if(!in){
    cerr << "Cannot open input file in readCddExtFile." << endl;
    exit(1);
  }

  while (tmpString!="begin") getline(in,tmpString);

  in >> numOfVertices >> numOfVars >> tmpString;

  cones=createListCone();
  endCones=cones;
  if(numOfVertices == 0){
    cout << "Empty Polytope." << endl;
    ofstream OUT("numOfLatticePoints");
    OUT << 0 << endl;
    exit(0);}

  if(numOfVertices == 1){
    char read = 'a';
    int flag = 0;
    ofstream OUT("numOfLatticePoints");
    // cout << tmpString << endl;
    //getline(in , tmpString);
    in.get(read); 
    //cout << tmpString << endl;
    
    while(((read == '\n' || read == '\r') || read == ' ') || read == '\t')
      {in.get( read); 
      if(read == '0') {	cout << "\n\nUnbounded polytope!" << endl << endl;
      exit(0);}}
    
    while(read != '\n' && read != '\r')
      {
	if(read == '/')
	  flag = 1;

	in.get( read); 
      }
    
    if(flag)
      {
	cout << "Integrally empty Polytope." << endl;
	OUT << 0 << endl;

      }
    
    else
      {
	cout << "\n\n*****  Total number of lattice points: " << 1 << " ****" << endl << endl;
	OUT << 1 << endl;

      }

    exit(0);}

  for (i=0; i<numOfVertices; i++) {
    v=createRationalVector(numOfVars-1);
    for (j=0; j<numOfVars; j++) {
      x=0;
      y=0;
      ReadCDD(in,x,y);
      if (j>0) {
	v->set_entry(j-1, x, y);
      } else {
	if (x==0) {
	  cout << "\n\nGiven polyhedron is unbounded!!!\n\n";
	  ofstream Empty("numOfLatticePoints");
	  Empty << 0 << endl;
	  exit(0);
	}
      }
    }
    c=createListCone();
    c->vertex = new Vertex(v);
    endCones->rest=c;
    endCones=endCones->rest;
  }

  in.close();

  cout << "done.\n";

  listCone *result = cones->rest;
  freeCone(cones);
  return result;
}
/* ----------------------------------------------------------------- */
listCone* readCddEadFile(listCone* cones, int numOfVars) {
  int i,j,k,numOfVertices,numOfRays;
  char cddInFileName[127];
  vec_ZZ v;
  rationalVector **vertices;
  listVector *rays, *endRays;
  listCone *tmp;
  string tmpString;

  cout << "Reading .ead file..."; cout.flush();

  strcpy(cddInFileName,"latte_cdd.ead");

  ifstream in(cddInFileName);
  if(!in){
    cerr << "Cannot open input file in readCddEadFile." << endl;
    exit(1);
  }

  while (tmpString!="begin") getline(in,tmpString);

  in >> numOfVertices;

  vertices=createArrayRationalVector(numOfVertices);

  tmp=cones;
  for (i=0; i<numOfVertices; i++) {
    vertices[i]=tmp->vertex->vertex;
    tmp=tmp->rest;
  }
  tmp=cones;
  for (i=0; i<numOfVertices; i++) {
    in >> k;
    if (i!=(k-1)) {
      cout << "Vertex numbering in file latte_cdd.ead is not increasing!\n";
      system_with_error_check("rm -f latte_cdd.*"); 
      exit(1);
    }

    in >> numOfRays;
    in >> tmpString;

    rays=createListVector(createVector(numOfVars));
    endRays=rays;

    for (j=0; j<numOfRays; j++) {
      in >> k;
      v=constructRay(vertices[i],vertices[k-1],numOfVars-1);
      endRays->rest=createListVector(v);
      endRays=endRays->rest;
    }

    tmp->rays=rays->rest;
    delete rays; // only deletes the dummy head
    tmp=tmp->rest;
  }
  delete[] vertices;

  in.close();
  cout << "done.\n";

  return (cones);
}

/* ----------------------------------------------------------------- */
void CreatExtEadFile(){
  char cddInFileName[127];
  string tmpString;

  strcpy(cddInFileName,"latte_cdd.out");
  ifstream in(cddInFileName);
  if(!in){
    cerr << "Cannot open input file in readCddEadFile." << endl;
    exit(1);
  }

  while (tmpString!="end") getline(in,tmpString);
  getline(in,tmpString);

  ofstream outExt("latte_cdd.ext");
  while (tmpString!="end"){
    getline(in,tmpString);
    outExt << tmpString << endl;}
  getline(in,tmpString);

  ofstream outEad("latte_cdd.ead");
  while (tmpString!="end"){
    getline(in,tmpString);
    outEad << tmpString << endl;}  

}

/* ----------------------------------------------------------------- */
listCone* readCddEadFileFromVrep(listCone* cones, int numOfVars) {
  int i,j,k,numOfVertices,numOfRays, counter = 0;
  char cddInFileName[127];
  vec_ZZ v;
  rationalVector **vertices;
  listVector *rays, *endRays;
  listCone *tmp;
  string tmpString;

  cout << "Reading .ead file..."; cout.flush();

  strcpy(cddInFileName,"latte_cdd.ead");
  int tmp_int;
  ifstream in(cddInFileName);
  if(!in){
    cerr << "Cannot open input file in readCddEadFile." << endl;
    exit(1);
  }

  while (tmpString!="begin") {getline(in,tmpString); counter++;  if(counter > 10) {cerr << "Redundant vertices!" << endl; exit(1);}}

  in >> numOfVertices >> tmp_int;

  vertices=createArrayRationalVector(numOfVertices);

  tmp=cones;
  for (i=0; i<numOfVertices; i++) {
    vertices[i]=tmp->vertex->vertex;
    tmp=tmp->rest;
  }
  tmp=cones;
  for (i=0; i<numOfVertices; i++) {
    in >> k;
    if (i!=(k-1)) {
      cout << "Vertex numbering in file latte_cdd.ead is not increasing!\n";
      system_with_error_check("rm -f latte_cdd.*"); 
      exit(1);
    }

    in >> numOfRays;
    in >> tmpString;

    rays=createListVector(createVector(numOfVars));
    endRays=rays;

    for (j=0; j<numOfRays; j++) {
      in >> k;
      v=constructRay(vertices[i],vertices[k-1],numOfVars-1);
      endRays->rest=createListVector(v);
      endRays=endRays->rest;
    }

    tmp->rays=rays->rest;
    tmp=tmp->rest;
  }

  in.close();
  cout << "done.\n";

  return (cones);
}


/* ----------------------------------------------------------------- */
listCone* computeVertexCones(char* fileName, listVector* matrix, 
			     int numOfVars) {
  char cddOutFileName[127], command[127];
  listCone *cones;

/* Compute vertices and edges with cdd. */

  createCddIneFile(matrix,numOfVars+1);

  cout << "Computing vertices and edges with cdd...";
  cout.flush();
  system_with_error_check(CDD_PATH " latte_cdd.ine > latte_cdd.out");
  cout << "done." << endl;

  strcpy(command,"cp latte_cdd.ext ");
  strcat(command,fileName);
  strcat(command,".ext");
  system_with_error_check(command);

  strcpy(command,"cp latte_cdd.ead ");
  strcat(command,fileName);
  strcat(command,".ead");
  system_with_error_check(command);

  {
    int ext_numOfVars;
    cones=readCddExtFile(ext_numOfVars);
    assert(ext_numOfVars == numOfVars+1);
  }
  cones=readCddEadFile(cones,numOfVars+1);
  system_with_error_check("rm -f latte_cdd.*"); 

  strcpy(cddOutFileName,fileName);
  strcat(cddOutFileName,".cdd");
  printListConeToFile(cddOutFileName,cones,numOfVars);

  return(cones);
}

/* ----------------------------------------------------------------- */

listCone* computeVertexConesViaLrs(char* fileName, listVector* matrix, 
			     int numOfVars) {

  char cddOutFileName[127], command[127];
  listCone *cones;

/* Compute vertices with lrs. */

  createLrsIneFile(matrix,numOfVars+1);

  cout << "Computing vertices with lrs...";
  system_with_error_check(LRS_PATH " latte_lrs.ine > latte_lrs.ext");
  cout << "done.\n\n";

  createLrsIneFileToPostAnalysys(matrix, numOfVars + 1);
  createLrsExtFileToPostAnalysys(matrix, numOfVars + 1);

  cout << "Computing edges with cdd...";
  system_with_error_check(CDD_PATH " latte_cdd.ine > latte_cdd.out");
  cout << "done.\n\n";

  strcpy(command,"cp latte_cdd.ext ");
  strcat(command,fileName);
  strcat(command,".ext");
  system_with_error_check(command);

  strcpy(command,"cp latte_cdd.ead ");
  strcat(command,fileName);
  strcat(command,".ead");
  system_with_error_check(command);

  {
    int ext_numOfVars;
    cones=readCddExtFile(ext_numOfVars);
    assert(ext_numOfVars == numOfVars+1);
  }
  cones=readCddEadFile(cones,numOfVars+1);
  system_with_error_check("rm -f latte_cdd.* latte_lrs.*"); 

  strcpy(cddOutFileName,fileName);
  strcat(cddOutFileName,".cdd");
  printListConeToFile(cddOutFileName,cones,numOfVars);

  return(cones);
}

/* ----------------------------------------------------------------- */
listCone* computeVertexConesFromVrep(const char* fileName, int &numOfVars) {
  char cddOutFileName[127], command[127];
  listCone *cones;

/* Compute vertices and edges with cdd. */

 createCddExtFile2(fileName);

  cout << "Computing vertices and edges with cdd...";
  system_with_error_check(COMPUTEADJACENCY_PATH " latte_cdd.ext > latte_cdd.jnk 2>&1");
  cout << "done.\n\n";
  //  CreatExtEadFile();
   strcpy(command,"cp latte_cdd.ext ");
  strcat(command,fileName);
  strcat(command,".ext");
  system_with_error_check(command);

  strcpy(command,"cp latte_cdd.ead ");
  strcat(command,fileName);
  strcat(command,".ead");
  system_with_error_check(command);

  {
    int ext_numOfVars;
    cones=readCddExtFile(ext_numOfVars);
    numOfVars = ext_numOfVars - 1;
  }
  cones=readCddEadFileFromVrep(cones,numOfVars+1);
  system_with_error_check("rm -f latte_cdd.*"); 

  strcpy(cddOutFileName,fileName);
  strcat(cddOutFileName,".cdd");
  printListConeToFile(cddOutFileName,cones,numOfVars);

  return(cones);
}

/* ----------------------------------------------------------------- */
rationalVector* LP(listVector* matrix, vec_ZZ& cost, int numOfVars,
		   bool verbose = true)
{
  rationalVector* Opt_vector;
  createCddIneLPFile(matrix,numOfVars+1,cost);
  if (verbose) {
    cout << "Computing LP... "; cout.flush();
  }
  system_with_error_check(CDD_PATH " LP.ine > LP.out");
  if (verbose) {
    cout << "done."; cout.flush();
  }
  Opt_vector = ReadLpsFile(numOfVars, verbose);
  //  cout << Opt_vector->numerators() << " " << Opt_vector -> denominator << endl;
  system_with_error_check("rm -f LP.*"); 

  return(Opt_vector);
}
/* ----------------------------------------------------------------- */
