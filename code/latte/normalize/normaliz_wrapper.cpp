/*
4ti2 -- A software package for algebraic, geometric and combinatorial
problems on linear spaces.

Copyright (C) 2006 4ti2 team.
Main author(s): Raymond Hemmecke.

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA. 
*/

/* ----------------------------------------------------------------- */
/*                                                                   */
/* Deals with MK's normaliz under symmetry                           */
/*                                                                   */
/* Author   : Raymond Hemmecke                                       */
/*                                                                   */
/* ----------------------------------------------------------------- */
#include <stdio.h>
#include <string.h>

// Stuff from lib4ti2util
extern "C" {
#include "util/myheader.h"
#include "util/print.h"
#include "util/orbit.h"
#include "util/output.h"
#include "util/vector.h"
}

// Silly interface to MK's normaliz
int normalize_commandline(char *command);

listVector* candidates;
int callsToNormaliz;

/* ----------------------------------------------------------------- */
listVector* myReadListVector(int *numOfVars, char *fileName) {
  int numOfVectors;
  listVector *basis, *endBasis;
  vector b;
  FILE *in;

  setbuf(stdout,0);
  if (!(in = fopen(fileName,"r"))) {
    printf("File \"%s\" not found for reading!\n",fileName);
    return(0);
  }

  fscanf(in,"%d",&numOfVectors);
  fscanf(in,"%d",numOfVars);

  if (numOfVectors==0) {
    fclose(in);
    return (0);
  }

  b=createVector(*numOfVars);
  for (int j=0; j<(*numOfVars); j++) fscanf(in,"%d",&b[j]);
  basis = createListVector(b);
  endBasis = basis;

  for (int i=1; i<numOfVectors; i++) {
    b=createVector(*numOfVars);
    for (int j=0; j<(*numOfVars); j++) fscanf(in,"%d",&b[j]);
    endBasis = updateBasis(createListVector(b), endBasis);
  }
  fclose(in);
  return(basis);
}
/* ----------------------------------------------------------------- */
void myPrintListVectorToFile(char* fileName, listVector* basis, int numOfVars) {
  int len;
  FILE* out;

  if (!(out = fopen(fileName,"w"))) {
    printf("Error opening output file!");
    exit (0);
  }
  if (basis==0) {
    fprintf(out,"0 %d\n",numOfVars);
    fclose(out);
    return;
  }

  len=lengthListVector(basis);
  fprintf(out,"%d %d\n",len,numOfVars);
  while(basis) {
    printVectorToFile(out,basis->first,numOfVars);
    basis = basis->rest;
  }
  fprintf(out,"\n");

  fclose(out);
  return ;
}
/* ----------------------------------------------------------------- */
listVector* extractSmallCones(listVector **mainCones, int threshold,
			      int numOfVars) {
  int norm;
  listVector *smallCones, *bigCones, *tmp, *tmp2;

  smallCones=0;
  bigCones=0;
  tmp=*mainCones;
  while (tmp) {
    norm=normOfVector(tmp->first,numOfVars);
    if (norm<threshold) {
      tmp2=tmp;
      tmp=tmp->rest;
      tmp2->rest=smallCones;
      smallCones=tmp2;
    } else {
      tmp2=tmp;
      tmp=tmp->rest;
      tmp2->rest=bigCones;
      bigCones=tmp2;
    }
  }
  *mainCones=bigCones;

  return (smallCones);
}
/* ----------------------------------------------------------------- */
listVector* extractSimplicialCones(listVector *simplicialCones, 
				   listVector **smallCones, int dimension,
				   int numOfVars) {
  int norm;
  listVector *smallConesLeft, *tmp, *tmp2;

  printf("numOfVars = %d\n",numOfVars);
  printf("dimension = %d\n",dimension);

  smallConesLeft=0;
  tmp=*smallCones;
  while (tmp) {
    norm=normOfVector(tmp->first,numOfVars);
    if (norm==dimension) {
      tmp2=tmp;
      tmp=tmp->rest;
      tmp2->rest=simplicialCones;
      simplicialCones=tmp2;
    } else {
      tmp2=tmp;
      tmp=tmp->rest;
      tmp2->rest=smallConesLeft;
      smallConesLeft=tmp2;
    }
  }
  *smallCones=smallConesLeft;

  return (simplicialCones);
}
/* ----------------------------------------------------------------- */
void runNormaliz(char *inFileName, char *outFileName, char *normaliz, 
		 char *raysFileName, int rayToBePulled, int trivial) {
  char command[1000];
  int retval;

/*    strcpy(command,"/home/mkoeppe/w/latte/build-gcc411-64/dest/bin/"); */

  if (normaliz[0])
    strcpy(command,normaliz);
  else
    strcpy(command, "dummy");
  
  strcat(command," --quiet --no-triang-file --only-triangulate");
  strcat(command," --triangulation-pull-rays=");
  sprintf(command,"%s%d",command,rayToBePulled);
  strcat(command," --subcones=");
  strcat(command,inFileName);
  strcat(command," --output-subcones=");
  strcat(command,outFileName);
  if (trivial==1) {
    strcat(command," --output-trivial-subcones=");
    strcat(command,outFileName);
    strcat(command,".trivial");
  }
  strcat(command," ");
  strcat(command,raysFileName);

  if (normaliz[0]) {
    callsToNormaliz++;
    printf("callsToNormaliz = %d\n",callsToNormaliz);
    strcat(command," >out.tmp");
    /*    printf("%s\n",command); */
    retval = system(command);
  }
  else
    retval = normalize_commandline(command);

  if (retval != 0) {
    fprintf(stderr, "Normaliz returned nonzero exit status.\n");
    exit(1);
  }
  return;
}
/* ----------------------------------------------------------------- */
void checkCones(char *simplicialConesFileName,char *checkFileName, 
		char* raysFileName, char *normaliz) {
  int numOfVars;
  listVector *HB, *tmp;
  char hilFileName[127],command[1000];
  int retval;

  printf("Checking simplicial cones.\n");

  if (normaliz[0])
    strcpy(command,normaliz);
  else
    strcpy(command, "dummy");

  strcat(command," --quiet --no-initial-triangulation");
  strcat(command," --no-triang-file --reduction=cplex");
  strcat(command," --max-determinant-for-enumeration=10000");
  strcat(command," --reduction-rays-file=");
  strcat(command,raysFileName);
  strcat(command,".check");
  strcat(command," --subcones=");
  strcat(command,simplicialConesFileName);
  strcat(command," ");
  strcat(command,raysFileName);

  printf("%s\n",command);
  if (normaliz[0]) {
    callsToNormaliz++;
    /*    printf("callsToNormaliz = %d\n",callsToNormaliz); */
    strcat(command," >out.tmp");
    /*    printf("%s\n",command); */
    do {
      retval = system(command);
      /*      printf("retval = %d\n",retval); */
      if (retval != 0) {
	fprintf(stderr, "No license, wtf?\n");
      }
    } while (retval != 0);
  } 
  else
    retval = normalize_commandline(command);
  if (retval != 0) {
    fprintf(stderr, "Normaliz returned nonzero exit status.\n");
    exit(1);
  }
  
  printf("CPLEX done.\n");

  strcpy(hilFileName,raysFileName);
  strcat(hilFileName,"--subcones-");
  strcat(hilFileName,simplicialConesFileName);
  strcat(hilFileName,".hil");

  HB=myReadListVector(&numOfVars,hilFileName);
  if (HB) {
    printf("New candidates for HB:\n");
    printListVector(HB,numOfVars);
  }

  if (HB) {
    tmp=HB;
    while (tmp->rest) tmp=tmp->rest;
    tmp->rest=candidates;
    candidates=HB;
  }

  strcpy(hilFileName,raysFileName);
  strcat(hilFileName,".hilbert");
  myPrintListVectorToFile(hilFileName,candidates,numOfVars);

  return;
}
/* ----------------------------------------------------------------- */
listVector* locallyPullingRay(listVector* smallCones, listVector* mainOrbits, 
			      char* smallConesInFileName, char* smallConesOutFileName, 
			      char* trivialSmallConesOutFileName, char* simplicialConesFileName, 
			      char* normaliz, char* raysFileName, int localRayToBePulled, 
			      int dimension, int numOfVars) {

  int i,k,len;
  listVector *trivialSmallCones, *newSmallCones, *allSmallCones, *simplicialCones, *tmp, *tmp2;

  simplicialCones=0;
  allSmallCones=0;
  trivialSmallCones=0;
  len=lengthListVector(smallCones);
   
  k=100000;
  while (smallCones) {
    printf("Remaining smallCones to be checked = %d\n",len);    
    if (len<k+1) {
      tmp=smallCones;
      smallCones=0;
      len=0;
    } else {
      tmp=smallCones;
      tmp2=smallCones;
      for (i=1;i<k;i++) tmp2=tmp2->rest;
      smallCones=tmp2->rest;
      tmp2->rest=0;
      len=len-k;
    }
    
    myPrintListVectorToFile(smallConesInFileName,tmp,numOfVars);

    if (tmp) {
      freeAllOfListVector(tmp);
      tmp=0;
    }
  
    printf("Locally pulling ray = %d\n",localRayToBePulled);
    runNormaliz(smallConesInFileName,smallConesOutFileName,normaliz,
                raysFileName,localRayToBePulled,1);
    newSmallCones=myReadListVector(&numOfVars,smallConesOutFileName);
    trivialSmallCones=myReadListVector(&numOfVars,trivialSmallConesOutFileName);
    if (trivialSmallCones) {
      tmp=trivialSmallCones;
      while (tmp->rest) tmp=tmp->rest;
      tmp->rest=allSmallCones;
      allSmallCones=trivialSmallCones;
      trivialSmallCones=0;
    }

    if (newSmallCones) {
      /*      printf("trivial new smallCones = %d, ",lengthListVector(trivialSmallCones));
	      printf("non-trivial new smallCones = %d -> ", lengthListVector(newSmallCones)); */
      simplicialCones=extractSimplicialCones(simplicialCones,&newSmallCones,dimension,numOfVars);
      newSmallCones=extractNonDominatedVectors(newSmallCones,mainOrbits,numOfVars);
      /*      printf("uncovered = %d -> ",lengthListVector(newSmallCones)); */

      if (newSmallCones) {
        tmp=newSmallCones;
        while (tmp->rest) tmp=tmp->rest;
        tmp->rest=allSmallCones;
        allSmallCones=newSmallCones;
        newSmallCones=0;
      }

      /*      printf("simplicial = %d\n",lengthListVector(simplicialCones)); */
      if (simplicialCones) { 
        myPrintListVectorToFile(simplicialConesFileName,simplicialCones,numOfVars);
        checkCones(simplicialConesFileName,raysFileName,raysFileName,normaliz);
        freeAllOfListVector(simplicialCones);
	simplicialCones=0;
      } 
    }
  }

  return (allSmallCones);
}
/* ----------------------------------------------------------------- */
listVector* pullOneRay(char* simplicialConesFileName, char* mainConesInFileName, 
            char* mainConesOutFileName, char* normaliz, char* raysFileName, 
            listVector* mainCones, listVector* symmGroup,
            int rayToBePulled, int numOfVars, int dimension) {
  int maxNorm,threshold,localRayToBePulled;
  listVector *simplicialCones, *smallCones, *mainOrbits;
  char mainConesInFileNameNumbered[127], mainConesOutFileNameNumbered[127],
    smallConesInFileName[127],smallConesOutFileName[127],trivialSmallConesOutFileName[127];
  
  if (mainCones==0) return (mainCones);

  if (dimension==0) {
    printf("Dimension of cone not specified!!!\n");
    exit(1);
  }

  mainOrbits=0;

  printf("\n=======================================================\n");
  printf("\nPulling ray = %d\n",rayToBePulled);
  printf("\n=======================================================\n\n");
  maxNorm=maximalNormInListVector(mainCones,numOfVars);
  if (maxNorm==dimension) {
    simplicialCones=mainCones;
    mainCones=0;
    myPrintListVectorToFile(simplicialConesFileName,simplicialCones,numOfVars);
    if (simplicialCones) { 
      checkCones(simplicialConesFileName,raysFileName,raysFileName,normaliz);
      freeAllOfListVector(simplicialCones);
      simplicialCones=0;
      return (mainCones);
    }
  } else {
    printListVector(mainCones,numOfVars);

    sprintf(mainConesInFileNameNumbered,"%s.%d",mainConesInFileName,rayToBePulled);
    sprintf(mainConesOutFileNameNumbered,"%s.%d",mainConesOutFileName,rayToBePulled);
    myPrintListVectorToFile(mainConesInFileNameNumbered,mainCones,numOfVars);
    freeAllOfListVector(mainCones);

    runNormaliz(mainConesInFileNameNumbered,mainConesOutFileNameNumbered,
    		    normaliz,raysFileName,rayToBePulled,0);
    mainCones=myReadListVector(&numOfVars,mainConesOutFileNameNumbered);

    threshold=maximalNormInListVector(mainCones,numOfVars);
    smallCones=extractSmallCones(&mainCones,threshold,numOfVars);
    if (smallCones==0) printf("No small cones.\n");

    if (smallCones) {
      printf("main cones = %d, small cones = %d\n",lengthListVector(mainCones),
            lengthListVector(smallCones));

      if (mainOrbits) freeAllOfListVector(mainOrbits);
      mainOrbits=expandRepresentativeIntoFullOrbits(mainCones,symmGroup,numOfVars,10);
      printf("mainOrbits = %d,   ",lengthListVector(mainOrbits));
          
      smallCones=extractNonDominatedVectors(smallCones,mainOrbits,numOfVars);
      printf("uncovered smallCones = %d -> ",lengthListVector(smallCones));
    	
      simplicialCones=extractSimplicialCones(simplicialCones,&smallCones,
    						                 dimension,numOfVars);
      printf("simplicial = %d\n",lengthListVector(simplicialCones));
      if (simplicialCones) { 
        checkCones(simplicialConesFileName,raysFileName,raysFileName,normaliz);
        freeAllOfListVector(simplicialCones);
        simplicialCones=0;
      }
        
      /* Replace this by a vector hasBeedPulled to make order of pulling more flexible. */

      localRayToBePulled=rayToBePulled;
      strcpy(smallConesInFileName,mainConesInFileName);
      strcat(smallConesInFileName,".smallcones.in");
      strcpy(smallConesOutFileName,mainConesInFileName);
      strcat(smallConesOutFileName,".smallcones.out");
      strcpy(trivialSmallConesOutFileName,mainConesInFileName);
      strcat(trivialSmallConesOutFileName,".smallcones.out.trivial");
          
      while (smallCones) {
        localRayToBePulled++;
        smallCones=locallyPullingRay(smallCones,mainOrbits,smallConesInFileName,
				     smallConesOutFileName,trivialSmallConesOutFileName,
				     simplicialConesFileName,normaliz,raysFileName,localRayToBePulled,
				     dimension,numOfVars);
      }
    }
  }

  if (mainOrbits) freeAllOfListVector(mainOrbits);

  return (mainCones);
}
/* ----------------------------------------------------------------- */
static void usage()
{
  fprintf(stderr, "new version\n");
  fprintf(stderr, "usage: normaliz_wrapper [OPTIONS...] FILENAME\n");
}
/* ----------------------------------------------------------------- */
int main(int argc, char *argv[]) {
  int i,rayToBePulled,localRayToBePulled,dimension,numOfVars,threshold,
    maxNorm,trivialPulling;
  vector v;
  listVector *mainCones, *symmGroup, *smallCones, *trivialSmallCones, 
    *simplicialCones, *tmp;
  char raysFileName[127],symFileName[127],mainConesInFileName[127],
    mainConesInFileNameNumbered[127],mainConesOutFileName[127],
    mainConesOutFileNameNumbered[127],smallConesInFileName[127],
    smallConesOutFileName[127],trivialSmallConesOutFileName[127],
    simplicialConesFileName[127],reductionRaysFileName[127],
    action[127],normaliz[127];

  if (argc < 2) {
    usage();
    exit(1);
  }

  callsToNormaliz=0;  

  normaliz[0] = '\0'; /* initialize... --mkoeppe */
  
  setbuf(stdout,0);

  strcpy(action, "pullall"); /* default */
  strcpy(raysFileName,argv[argc-1]);
  strcpy(symFileName,argv[argc-1]);
  strcat(symFileName,".sym.full");
  strcpy(simplicialConesFileName,argv[argc-1]);
  strcat(simplicialConesFileName,".simplicial");
  strcpy(action,"pullall");
  rayToBePulled=0;
  threshold=-1;
  dimension=0;
  strcpy(reductionRaysFileName,"");
  strcpy(mainConesInFileName,"");
  strcpy(mainConesOutFileName,"");
  strcpy(mainConesInFileNameNumbered,"");
  strcpy(mainConesOutFileNameNumbered,"");

  for (i=1;i<argc-1;i++) {
    if (strncmp(argv[i], "--symmetry-file",15) == 0) {
	strcpy(symFileName,argv[i]+16);
      } else if (strncmp(argv[i], "--main-cones-in-file",20) == 0) {
	strcpy(mainConesInFileName,argv[i]+21);
      } else if (strncmp(argv[i], "--main-cones-out-file",21) == 0) {
	strcpy(mainConesOutFileName,argv[i]+22);
      } else if (strncmp(argv[i], "--simplicial-cones-file",23) == 0) {
	strcpy(simplicialConesFileName,argv[i]+24);
      } else if (strncmp(argv[i], "--reduction-rays-file",21) == 0) {
	strcpy(reductionRaysFileName,argv[i]+22);
      } else if (strncmp(argv[i], "--triangulation-pull-rays",25) == 0) {
	rayToBePulled=atoi(argv[i]+26);
      } else if (strncmp(argv[i], "--dimension",11) == 0) {
	dimension=atoi(argv[i]+12);
      } else if (strncmp(argv[i], "--main-cones-threshold",22) == 0) {
	threshold=atoi(argv[i]+23);
      } else if (strncmp(argv[i], "--action",8) == 0) {
	strcpy(action,argv[i]+9);
      } else if (strncmp(argv[i], "--normaliz",10) == 0) {
	strcpy(normaliz,argv[i]+11);
      } 
  }
  symmGroup=myReadListVector(&numOfVars,symFileName);
  if (symmGroup==0) {
    v=createVector(numOfVars);
    for (i=0;i<numOfVars;i++) v[i]=i;
    symmGroup=createListVector(v);
  }
  if (mainConesInFileName[0]=='\0') {
    strcpy(mainConesInFileName,argv[argc-1]);
    strcat(mainConesInFileName,".mainCones.in");
    sprintf(mainConesInFileNameNumbered,"%s.%d",mainConesInFileName,1);
    v=createVector(numOfVars);
    for (i=0;i<numOfVars;i++) v[i]=1;
    myPrintListVectorToFile(mainConesInFileNameNumbered,createListVector(v),
			  numOfVars);
  } else {
    sprintf(mainConesInFileNameNumbered,"%s.%d",mainConesInFileName,1);
  }
  if (mainConesOutFileName[0]=='\0') {
    strcpy(mainConesOutFileName,argv[argc-1]);
    strcat(mainConesOutFileName,".mainCones.out");
  }
  printf("symFileName = %s\n",symFileName);
  printf("mainConesInFileName = %s\n",mainConesInFileName);
  printf("mainConesOutFileName = %s\n",mainConesOutFileName);
  printf("simplicialConesFileName = %s\n",simplicialConesFileName);
  printf("ray to be pulled = %d\n",rayToBePulled);
  printf("dimension = %d\n",dimension);
  printf("threshold = %d\n",threshold);
  printf("action = %s\n",action);
  printf("normaliz = %s\n",normaliz);

  mainCones=myReadListVector(&numOfVars,mainConesInFileNameNumbered);
  smallCones=0;
  trivialSmallCones=0;
  simplicialCones=0;
  candidates=0;
  trivialPulling=0;

  if (strncmp(action,"pullRay",7)==0) {
    mainCones=pullOneRay(simplicialConesFileName,mainConesInFileName, 
                         mainConesOutFileName,normaliz,raysFileName,mainCones, 
                         symmGroup,rayToBePulled,numOfVars,dimension);
  } else
  if (strncmp(action,"pullall",7)==0) {
    rayToBePulled=0;
    while (mainCones) {
      rayToBePulled++;
      mainCones=pullOneRay(simplicialConesFileName,mainConesInFileName, 
                           mainConesOutFileName,normaliz,raysFileName,mainCones, 
                           symmGroup,rayToBePulled,numOfVars,dimension);
    }
  } else {
    fprintf(stderr, "--action=%s not handled.\n", action);
    exit(1);
  }

  return(0);
}
/* ----------------------------------------------------------------- */
