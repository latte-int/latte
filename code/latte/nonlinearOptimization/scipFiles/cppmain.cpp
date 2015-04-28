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

	if ( strcmp(argv[1], "--int") == 0)
		type = 2;
	else if (strcmp(argv[1], "--sum") == 0)
		type = 1;
	else
	{
		cout << "usage: " << argv[0] << " --int | --sum -f filename\n"
				"Where --int is for integration and --sum is for summation" <<endl;
		exit(1);
	}

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
   if (type == 1)
	   SCIP_CALL( SCIPincludeObjRelax(scip, new LatteSummationRelaxor(scip, fileName), TRUE) );
   else if (type == 2)
	   SCIP_CALL( SCIPincludeObjRelax(scip, new LatteIntegrationRelaxor(scip, fileName), TRUE) );
   else
   {
	   cout << "this should never happen" <<endl;
   }

   /**********************************
    * Process command line arguments *
    **********************************/

   //remove the --sum/--int arg.
   SCIP_CALL( SCIPprocessShellArguments(scip, argc-1, argv+1, defaultsetname) );


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
