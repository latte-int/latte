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
listVector* checkCones(listVector *candidates,char *simplicialConesFileName,
		       char *checkFileName, char* raysFileName, 
		       char *normaliz) {
  int numOfVars;
  listVector *HB, *tmp;
  char hilFileName[127],command[1000];
  int retval;

  if (normaliz[0])
    strcpy(command,normaliz);
  else
    strcpy(command, "dummy");

  strcat(command," --quiet --no-initial-triangulation");
  strcat(command," --no-triang-file --reduction=cplex");
  strcat(command," --max-determinant-for-enumeration=10000");
  strcat(command," --reduction-rays-file=");
  strcat(command,raysFileName);
  strcat(command," --subcones=");
  strcat(command,simplicialConesFileName);
  strcat(command," ");
  strcat(command,raysFileName);

  if (normaliz[0]) {
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
  
  strcpy(hilFileName,raysFileName);
  strcat(hilFileName,"--subcones-");
  strcat(hilFileName,simplicialConesFileName);
  strcat(hilFileName,".hil");

  HB=readListVector(&numOfVars,hilFileName);

  if (HB) {
    tmp=HB;
    while (tmp->rest) tmp=tmp->rest;
    tmp->rest=candidates;
    candidates=HB;
  }

  strcpy(hilFileName,raysFileName);
  strcat(hilFileName,".hilbert");
  printListVectorToFile(hilFileName,candidates,numOfVars);

  return (candidates);
}

static void usage()
{
  fprintf(stderr, "usage: normaliz_wrapper [OPTIONS...] FILENAME\n");
}


/* ----------------------------------------------------------------- */
int main(int argc, char *argv[]) {
  int i,rayToBePulled,localRayToBePulled,dimension,numOfVars,threshold,
    maxNorm,trivialPulling;
  vector v;
  listVector *mainCones, *symmGroup, *smallCones, *trivialSmallCones, 
    *mainOrbits, *simplicialCones, *tmp, *candidates;
  char raysFileName[127],symFileName[127],mainConesInFileName[127],
    mainConesOutFileName[127],smallConesInFileName[127],
    smallConesOutFileName[127],trivialSmallConesOutFileName[127],
    simplicialConesFileName[127],action[127],normaliz[127];

  if (argc < 2) {
    usage();
    exit(1);
  }
  
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
  strcpy(mainConesInFileName,"");
  strcpy(mainConesOutFileName,"");

  for (i=1;i<argc-1;i++) {
    if (strncmp(argv[i], "--symmetry-file",15) == 0) {
	strcpy(symFileName,argv[i]+16);
      } else if (strncmp(argv[i], "--main-cones-in-file",20) == 0) {
	strcpy(mainConesInFileName,argv[i]+21);
      } else if (strncmp(argv[i], "--main-cones-out-file",21) == 0) {
	strcpy(mainConesOutFileName,argv[i]+22);
      } else if (strncmp(argv[i], "--simplicial-cones-file",23) == 0) {
	strcpy(simplicialConesFileName,argv[i]+24);
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
  symmGroup=readListVector(&numOfVars,symFileName);
  if (symmGroup==0) {
    v=createVector(numOfVars);
    for (i=0;i<numOfVars;i++) v[i]=i;
    symmGroup=createListVector(v);
  }
  if (mainConesInFileName[0]=='\0') {
    strcpy(mainConesInFileName,argv[argc-1]);
    strcat(mainConesInFileName,".mainCones.in");
    v=createVector(numOfVars);
    for (i=0;i<numOfVars;i++) v[i]=1;
    printListVectorToFile(mainConesInFileName,createListVector(v),numOfVars);
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

  mainCones=readListVector(&numOfVars,mainConesInFileName);
  mainOrbits=0;
  smallCones=0;
  trivialSmallCones=0;
  simplicialCones=0;
  candidates=0;
  trivialPulling=0;

  rayToBePulled=0;
  if (strncmp(action,"pullall",7)==0) {
    if (dimension==0) {
      printf("Dimension of cone not specified!!!\n");
      exit(1);
    }
    while (mainCones) {
      rayToBePulled++;
      printf("\n=======================================================\n");
      printf("\nPulling ray = %d\n",rayToBePulled);
       printf("\n=======================================================\n\n");
      maxNorm=maximalNormInListVector(mainCones,numOfVars);
      if (maxNorm==dimension) {
	tmp=mainCones;
	while (tmp->rest) tmp=tmp->rest;
	tmp->rest=simplicialCones;
	simplicialCones=mainCones;
	mainCones=0;
	printListVectorToFile(simplicialConesFileName,simplicialCones,
			      numOfVars);
	if (simplicialCones) { 
	  candidates=checkCones(candidates,simplicialConesFileName,
				raysFileName,raysFileName,normaliz);

	  freeAllOfListVector(simplicialCones);
	  simplicialCones=0;
	}
      } else {
	printListVector(mainCones,numOfVars);
	runNormaliz(mainConesInFileName,mainConesOutFileName,normaliz,
		    raysFileName,rayToBePulled,0);
	mainCones=readListVector(&numOfVars,mainConesOutFileName);
	threshold=maximalNormInListVector(mainCones,numOfVars);
	smallCones=extractSmallCones(&mainCones,threshold,numOfVars);
	printListVectorToFile(mainConesInFileName,mainCones,numOfVars);
	if (smallCones) {
	  printf("main cones = %d, small cones = %d\n",
		 lengthListVector(mainCones),lengthListVector(smallCones));
	  
	  if (strncmp(mainConesInFileName,"346",3)==0) {
	    tmp=mainCones;
	    while (tmp) {
	      (tmp->first)[2]=1;
	      (tmp->first)[5]=1;
	      (tmp->first)[14]=1;
	      (tmp->first)[17]=1;
	      (tmp->first)[26]=1;
	      tmp=tmp->rest;
	    }
	  }
	  if (strncmp(mainConesInFileName,"355",3)==0) {
	    tmp=mainCones;
	    while (tmp) {
	      (tmp->first)[2]=1;
	      (tmp->first)[5]=1;
	      (tmp->first)[8]=1;
	      (tmp->first)[17]=1;
	      tmp=tmp->rest;
	    }
	  }
	  if (mainOrbits) freeAllOfListVector(mainOrbits);
	  mainOrbits=expandRepresentativeIntoFullOrbits(mainCones,symmGroup,
							numOfVars,10);
	  if (strncmp(mainConesInFileName,"346",3)==0) {
	    tmp=mainCones;
	    while (tmp) {
	      (tmp->first)[2]=0;
	      (tmp->first)[5]=0;
	      (tmp->first)[14]=0;
	      (tmp->first)[17]=0;
	      (tmp->first)[26]=0;
	      tmp=tmp->rest;
	    }
	  }
	  if (strncmp(mainConesInFileName,"355",3)==0) {
	    tmp=mainCones;
	    while (tmp) {
	      (tmp->first)[2]=0;
	      (tmp->first)[5]=0;
	      (tmp->first)[8]=0;
	      (tmp->first)[17]=0;
	      tmp=tmp->rest;
	    }
	  }
	  printf("mainOrbits = %d,   ",lengthListVector(mainOrbits));
      
	  smallCones=extractNonDominatedVectors(smallCones,mainOrbits,
						numOfVars);
	  printf("uncovered smallCones = %d -> ",lengthListVector(smallCones));
	
	  simplicialCones=extractSimplicialCones(simplicialCones,&smallCones,
						 dimension,numOfVars);
	  printf("simplicial = %d\n",lengthListVector(simplicialCones));
    
	  /* Replace by a vector hasBeedPulled to make order of pulling 
	     more flexible. */
	  localRayToBePulled=rayToBePulled;
	  strcpy(smallConesInFileName,mainConesInFileName);
	  strcat(smallConesInFileName,".smallcones.in");
	  strcpy(smallConesOutFileName,mainConesInFileName);
	  strcat(smallConesOutFileName,".smallcones.out");
	  strcpy(trivialSmallConesOutFileName,mainConesInFileName);
	  strcat(trivialSmallConesOutFileName,".smallcones.out.trivial");
      
	  while (smallCones) {
/*  	    if (trivialPulling==0) */
	      printListVectorToFile(smallConesInFileName,smallCones,
				    numOfVars);
	    if (smallCones) {
	      freeAllOfListVector(smallCones);
	      smallCones=0;
	    }
	    localRayToBePulled++;
	    if (localRayToBePulled==rayToBePulled) localRayToBePulled++;
	    
	    printf("Locally pulling ray = %d\n",localRayToBePulled);
	    runNormaliz(smallConesInFileName,smallConesOutFileName,normaliz,
			raysFileName,localRayToBePulled,1);
	    smallCones=readListVector(&numOfVars,smallConesOutFileName);
	    trivialSmallCones=readListVector(&numOfVars,
					     trivialSmallConesOutFileName);
	    if (trivialSmallCones) {
	      trivialPulling=1;
	    } else {
	      trivialPulling=0;
	    }
	    if (smallCones) {
	      printf("trivial smallCones = %d, ",
		     lengthListVector(trivialSmallCones));
	      printf("non-trivial smallCones = %d -> ",
		     lengthListVector(smallCones));
	      simplicialCones=extractSimplicialCones(simplicialCones,
						     &smallCones,dimension,
						     numOfVars);
	      smallCones=extractNonDominatedVectors(smallCones,mainOrbits,
						    numOfVars);
	      printf("uncovered = %d -> ",lengthListVector(smallCones));
	      printf("nonsimplicial = %d and ",lengthListVector(smallCones));
	      printf("simplicial = %d\n",lengthListVector(simplicialCones));
	      printListVectorToFile(simplicialConesFileName,simplicialCones,
				    numOfVars);
	      if (simplicialCones) { 
		candidates=checkCones(candidates,simplicialConesFileName,
				      raysFileName,raysFileName,normaliz);
		freeAllOfListVector(simplicialCones);
		simplicialCones=0;
	      }
	    } else {
	      printf("No new cones. Trivial smallCones = %d\n",
		     lengthListVector(trivialSmallCones));
	    }
	    if (smallCones) {
	      tmp=smallCones;
	      while (tmp->rest) tmp=tmp->rest;
	      tmp->rest=trivialSmallCones;
	      trivialSmallCones=0;
	    } else {
	      smallCones=trivialSmallCones;
	      trivialSmallCones=0;
	    }
	  }
	}
      }
    }  
  }
  else {
    fprintf(stderr, "--action=%s not handled.\n", action);
    exit(1);
  }

  return(0);
}
/* ----------------------------------------------------------------- */
