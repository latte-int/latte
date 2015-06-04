/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2014 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cppmain.cpp
 * @brief  main file for C++ compilation
 * @author Tobias Achterberg
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <string>
#include <iostream>

//regular scip headers
#include "scip/scip.h"
#include "scip/scipshell.h"
#include "scip/scipdefplugins.h"
#include "scip/scipshell.h"
#include "scip/message_default.h"
#include "objscip/objrelax.h"

//my scip headers
#include "LatteSummationRelaxor.h"
#include "LatteIntegrationRelaxor.h"
#include "LatteBranchRule.h"

//latte headers
#include "barvinok/barvinok.h"
#include "ReadPolyhedron.h"
#include "LattException.h"
#include "nonlinearOptimization/WeightedExponentialSubs.h"
#include "nonlinearOptimization/BoxOptimization.h"
#include "print.h"
#include "dual.h"
#include "rational.h"
#include "integration/burstTrie.h"
#include "integration/PolyTrie.h"

using namespace std;


/** creates a SCIP instance with default plugins, evaluates command line parameters, runs SCIP appropriately,
 *  and frees the SCIP instance
 */
SCIP_RETCODE SCIPrunMyShell(
   int                   argc,               /**< number of shell parameters */
   char**                argv,               /**< array with shell parameters */
   const char*           defaultsetname      /**< name of default settings file */
   )
{
	SCIP* scip = NULL;
	char * fileName = NULL;
	int type = 0;

	for(int i = 0; i < argc; ++i)
		if ( strcmp(argv[i], "-f") == 0 && i + 1 < argc)
			fileName = argv[i+1];

	cout << "Usage  : " << argv[0] <<" <--int | --sum | --none for integration or summation> <--bb | --none for branch and bound> -f filename" << endl;
	cout << "Example: " << argv[0] << " --sum --none -f file.pip" <<endl;

	if ( fileName == NULL )
	{
		cout << "pip file missing. use ./scip -f file.pip where the 2nd line has a latte-style polynomial" << endl;
		exit(1);
	}

   /*********
    * Setup *
    *********/

   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) );

   /* include default SCIP plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* include latte plugins */
   ///*
   if (strcmp(argv[1], "--sum") == 0)
	   SCIP_CALL( SCIPincludeObjRelax(scip, new LatteSummationRelaxor(scip, fileName), TRUE) );
   else if (strcmp(argv[1], "--int") == 0)
	   SCIP_CALL( SCIPincludeObjRelax(scip, new LatteIntegrationRelaxor(scip, fileName), TRUE) );
   else
   {
	   cout << "Latte NOT being used" <<endl;
   }
   //*/

   /* include latte branch rule */
   if (strcmp(argv[2], "--bb") == 0)
   {
	   SCIP_CALL( SCIPincludeObjBranchrule(scip, new LatteBranchRule(scip), TRUE) );
       cout << "Latte BB used" << endl;
   }

   /**********************************
    * Process command line arguments *
    **********************************/

   //remove the --sum/--int arg.
   SCIP_CALL( SCIPprocessShellArguments(scip, argc-2, argv+2, defaultsetname) );


   /********************
    * Deinitialization *
    ********************/

   SCIP_CALL( SCIPfree(&scip) );

   BMScheckEmptyMemory();

   return SCIP_OKAY;
}

/** main method starting SCIP */
int main(
   int                        argc,          /**< number of arguments from the shell */
   char**                     argv           /**< array of shell arguments */
   )
{

/*************************/
   SCIP_RETCODE retcode;

   /* run interactive shell */
   retcode = SCIPrunMyShell(argc, argv, "scip.set");

   /* evaluate retrun code of the SCIP process */
   if( retcode != SCIP_OKAY )
   {
      /* write error back trace */
      SCIPprintError(retcode);
      return -1;
   }

   return 0;
}
