/* ----------------------------------------------------------------- */
/*                                                                   */
/* LattE (Lattice Point Enumeration)                                 */
/*                                                                   */
/* Computing all lattice points in fundamental parallelepiped        */
/*                                                                   */
/* Author     : Raymond Hemmecke, Ruriko Yoshida                     */
/*                                                                   */
/* Created    : 07-JUN-02                                            */
/* Last Update: 20-DEC-02 by Rudy                                    */
/*                                                                   */
/* ----------------------------------------------------------------- */
#include "config.h"
#include "../myheader.h"
#include "../cone.h"
#include "../print.h"
#include "../ramon.h"
#include "../rational.h"
#include "../vertices/cdd.h"
#include <NTL/ZZ.h>
#include <string>

using namespace std;
/* ----------------------------------------------------------------- */
vector movePoint(vector x, rationalVector *coeffsX, 
		 rationalVector *coeffsVertex, vector *matrix, 
		 int numOfRays, int numOfVars) {
  int i,j;
  vector z, movement;
  rationalVector *difference;

  difference=subRationalVector(coeffsX,coeffsVertex,numOfVars);

  movement=createVector(numOfRays);
  for (i=0; i<numOfRays; i++) {
    if (difference->denominator[i]==1) 
      movement[i]=-difference->enumerator[i];
    else {
      movement[i]=-(difference->enumerator[i]/difference->denominator[i]);
      if (difference->enumerator[i]<0) movement[i]=movement[i]+1;
    }
  }

  z=copyVector(x,numOfVars);

  for (i=0; i<numOfRays; i++) {
    if (movement[i]!=0) {
      for (j=0; j<numOfVars; j++) 
/* Note that each COLUMN of matrix constitutes an extreme ray, not
   each row! */
	z[j]=z[j]+movement[i]*matrix[j][i];
    }
  }

  return (z);
}
/* ----------------------------------------------------------------- */
listVector* pointsInParallelepiped(rationalVector *vertex, listVector *rays, 
				   listVector *facets, int numOfVars) {
  int i,j,k,ii,counter,numOfVertices,numOfRays,nextIndex,tmpInt;
  rationalVector *c, *w, *coeffsVertex;
  rationalVector **coeffs;
  vector z;
  vector *points, *matrix, *originalMatrix;
  listVector *tmp, *endFacets, *hilbertBasis, *interiorHilbertBasis,
      *listOfPoints, *endListOfPoints;
  char cddInFileName[127],mlpInFileName[127];
  string tmpString;
  ZZ x,y;

  setbuf(stdout,0);

  if (facets==0) {
    strcpy(cddInFileName,"latte_cdd.ine");
    ofstream out(cddInFileName);
    if(!out){
      cerr << "Cannot open output file latte_cdd.ine in pointsInParallelepiped." << endl;
      exit(1);
    }

    out << "H-representation\n";
    out << "begin\n";
    out << lengthListVector(rays) << " " << numOfVars+1 << " integer\n";  
    
    tmp=rays;  
    while (tmp) {
      out << "0 ";
      for (i=0; i<(numOfVars); i++) out << (tmp->first)[i] << " ";
      out << endl;
      tmp=tmp->rest;
    }

    out << "end\n";
    out.close();

//      printf("Computing facets with cdd...");
    system(CDD_PATH " latte_cdd.ine > latte_cdd.out");
//      printf("done.\n");

    strcpy(cddInFileName,"latte_cdd.ext");
    ifstream in(cddInFileName);
    if(!in){
      cerr << "Cannot open input file in pointsInParallelepiped." << endl;
      exit(1);
    }

    while (tmpString!="begin") getline(in,tmpString);

    in >> numOfVertices >> tmpInt >> tmpString;

    z=createVector(numOfVars);
    facets=createListVector(z);
    endFacets=facets;

    for (i=0; i<numOfVertices; i++) {
      w=createRationalVector(numOfVars);
      for (j=0; j<numOfVars+1; j++) {
	x=0;
	y=0;
	ReadCDD(in,x,y);
	if (j>0) {
	  w->enumerator[j-1]=x;
	  w->denominator[j-1]=y;
	}
      }
      w=normalizeRationalVector(w,numOfVars);
      endFacets->rest=createListVector(w->enumerator);
      endFacets=endFacets->rest;
    }

    facets=facets->rest;
    in.close();
    system("rm latte_cdd.*");
  }

/* Facets are all known at this points. */

  strcpy(mlpInFileName,"latte_mlp");
  ofstream out(cddInFileName);
  if (!out) {
    printf("Cannot open output file latte_mlp in pointsInParallelepiped!");
    exit(0);
  }

  out << numOfVars << " " << lengthListVector(facets) << endl;
  printListVectorToFileWithoutBrackets(out,facets,numOfVars);
  out.close();

  printf("Computing Hilbert basis with mlp...");
  system("./mlp dual latte_mlp >latte_mlp.out");
  printf("done.\n");

  strcpy(mlpInFileName,"latte_mlp.dual.hil");
  printf("Reading Hilbert basis from file...");
  hilbertBasis=readListVectorMLP(mlpInFileName,&ii);
  printf("done.\n");
/*   rays=readListVector("latte_mlp.dual.ray"); */

/*    system("rm latte_mlp*"); */

  numOfRays=lengthListVector(rays);

  matrix=createArrayVector(numOfRays);
  for (i=0; i<numOfVars; i++) matrix[i]=createVector(numOfRays);

  k=0;
  tmp=rays;
  while (tmp) {
    for (i=0; i<numOfVars; i++) matrix[i][k]=(tmp->first)[i];
    k++;
    tmp=tmp->rest;
  }

  points=createArrayVector(2000000);
  coeffs=createArrayRationalVector(2000000);

  tmp=hilbertBasis;
  counter=0;

  printf("rays:\n");
  printListVector(rays,numOfVars);

  while (tmp) {
    if (isVectorInListVector(tmp->first,rays,numOfVars)==0) {
      points[counter]=tmp->first;
      coeffs[counter]=solveLinearSystem(matrix,tmp->first,numOfVars,numOfRays);
      counter++;
    }
    tmp=tmp->rest;
  }

  originalMatrix=createArrayVector(numOfVars);
  for (i=0; i<numOfRays; i++) 
    originalMatrix[i]=copyVector(matrix[i],numOfRays);

  for (i=0; i<numOfVars; i++) {
    for (j=0; j<numOfRays; j++) {
      matrix[i][j]=matrix[i][j]*(vertex->denominator[i]);
    }
  }

  coeffsVertex=solveLinearSystem(matrix,vertex->enumerator,numOfVars,
				 numOfRays);

  if (counter==0) {
/* No interior HB element -> move origin as unique single point. */
    z=createVector(numOfVars);
    for (i=0;i<numOfVars;i++) z[i]=0;
    return(createListVector(movePoint(z,0,coeffsVertex,originalMatrix,
				      numOfRays,numOfVars)));
  }

  printf("Enumerating lattice points...");

/* Compute all interior points. */
  listOfPoints=transformArrayVectorToListVector(points,counter);
  interiorHilbertBasis=transformArrayVectorToListVector(points,counter);
  tmp=listOfPoints;
  while (tmp->rest) tmp=tmp->rest;
  endListOfPoints=tmp;

  printf("interior HB elements = %d / %d\n",lengthListVector(interiorHilbertBasis),lengthListVector(hilbertBasis));

/* Counter gives the number of interior HB elements. */

  nextIndex=counter;
  counter=0;

  printf("\nStart while loop.\n");
/*    printf("\nRays.\n"); */
/*    printListVector(rays,numOfVars); */
/*    printf("\nHilbert basis.\n"); */
/*    printListVector(hilbertBasis,numOfVars); */

  while (counter<nextIndex) {
    tmp=interiorHilbertBasis;
    k=0;
    while (tmp) {
      z=copyVector(tmp->first,numOfVars);
      z=addVector(z,points[counter],numOfVars);
      if (isVectorInListVector(z,listOfPoints,numOfVars)==0) {
	c=addRationalVectorsWithUpperBoundOne(coeffs[k],coeffs[counter],
					      numOfRays);
	if (c!=0) {
	  points[nextIndex]=z;
	  coeffs[nextIndex]=c;
	  endListOfPoints->rest=createListVector(z);
	  endListOfPoints=endListOfPoints->rest;
	  nextIndex++;
	  if (nextIndex==500*(nextIndex/500)) {
	    printf("counter = %d / %d.\n",counter,nextIndex);
	  }
	}
      }else z.kill();
      k++;
      tmp=tmp->rest;
    }
    counter++;
    if (counter==10*(counter/10)) {
      printf("counter = %d / %d.\n",counter,nextIndex);
    }
  }

/* Add the origin to list of points in parallelepiped. */
  z=createVector(numOfVars);
  for (i=0;i<numOfVars;i++) z[i]=0;
  points[nextIndex]=z;
  coeffs[nextIndex]=0;
  endListOfPoints->rest=createListVector(z);
  endListOfPoints=endListOfPoints->rest;

  printf("done.\n");

/* Schon mal nicht schlecht, aber jetzt muessen die Punkte noch 
   verschoben werden... */

  printf("Moving lattive points...");

  k=0;
  tmp=listOfPoints;
  while (tmp) {
      tmp->first=movePoint(tmp->first,coeffs[k],coeffsVertex,
			   originalMatrix,numOfRays,numOfVars);
      k++;
      tmp=tmp->rest;
  }

  printf("done.\n");

  return(listOfPoints);
}
/* ----------------------------------------------------------------- */
listVector* pointsInParallelepipedOfUnimodularCone(rationalVector *vertex, 
						  listVector *rays, 
						  int numOfVars) {
  int i,j,k,numOfRays;
  ZZ a,b;
  vector lambda,w,z;
  vector *matrix, *originalMatrix;
  listVector *points, *tmp;
  rationalVector *coeffs;

/*  printf("Computing point in parallelepiped of unimodular cone.\n"); */

/*  printRationalVector(vertex,numOfVars); */

  //points=createListVector(createVector(numOfVars));

  numOfRays=lengthListVector(rays);
  if (numOfRays!=numOfVars) {
    printf("Cone is NOT simplicial!\n");
    exit(0);
  }

/*  printf("numOfRays = %d, numOfVars = %d\n", numOfRays, numOfVars); */

  matrix=createArrayVector(numOfVars);
  originalMatrix=createArrayVector(numOfVars);

  for (i=0; i<numOfVars; i++) matrix[i]=createVector(numOfRays);

  k=0;
  tmp=rays;
  while (tmp) {
    for (i=0; i<numOfVars; i++) matrix[i][k]=(tmp->first)[i];
    k++;
    tmp=tmp->rest;
  }

 for (i=0; i<numOfRays; i++) 
    originalMatrix[i]=copyVector(matrix[i],numOfRays);

  for (i=0; i<numOfVars; i++) {
    for (j=0; j<numOfRays; j++) {
      matrix[i][j]=matrix[i][j]*(vertex->denominator[i]);
    }
  }
  w=copyVector(vertex->enumerator,numOfVars);

  /* Now we have to solve: matrix * x = w */

  coeffs=solveLinearSystem(matrix,w,numOfVars,numOfRays);

  lambda=createVector(numOfRays);

/* Note that we make heavy use of numOfVars=numOfRays! That is, we
   assume the cone to be simplicial! */

  for (i=0; i<numOfVars; i++) {
    a=coeffs->denominator[i];
    b=coeffs->enumerator[i];

    if (((b/a)*a)==b) {
      lambda[i]=b/a;
    } else {
      if (((b<0) && (a>0)) || ((b>0) && (a<0))) {
//          lambda[i]=-abs(b)/abs(a); 
         lambda[i]=abs(b)/abs(a); 
         lambda[i]=-lambda[i]; 
     } else 
	lambda[i]=b/a;
      if (b>0) lambda[i]=lambda[i]+1;
    }
  }

//  cout << "lambda = ";
//  printVector(lambda,numOfVars);

  z=createVector(numOfVars);
  for (i=0; i<numOfVars; i++) z[i]=0;

  for (i=0; i<numOfRays; i++) {
    for (j=0; j<numOfVars; j++) {
      z[j]=z[j]+lambda[i]*originalMatrix[j][i];
    }
  }

//  cout << "point = ";
//  printVector(z,numOfVars);

  //points->rest=createListVector(z);
  
  coeffs->denominator.kill ();
  coeffs->enumerator.kill ();

  delete coeffs;
  
  points=createListVector(z);
  delete [] matrix;
  delete [] originalMatrix;
  w.kill();
  z.kill();
  lambda.kill();

  return(points);
}
/* ----------------------------------------------------------------- */
