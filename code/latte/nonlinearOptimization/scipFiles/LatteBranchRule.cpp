/*
 * LatteBranchRule.cpp
 *
 *  Created on: May 8, 2015
 *      Author: bedutra
 */

#include <assert.h>
#include <iostream>
#include "LatteBranchRule.h"

using namespace std;

LatteBranchRule::LatteBranchRule(SCIP* scip): ObjBranchrule(scip, "LatteBranchRule", "branch on box variables", 4000, -1, 1)
{


}

LatteBranchRule::~LatteBranchRule() {
	// TODO Auto-generated destructor stub
}




/*
 * Local methods
 */

/* put your local methods here, and declare them static */


/*
 * Callback methods of branching rule
 */


/** destructor of branching rule to free user data (called when SCIP is exiting)
 *
 *  @see SCIP_DECL_BRANCHFREE(x) in @ref type_branch.h
 */
SCIP_DECL_BRANCHFREE(LatteBranchRule::scip_free)
{  /*lint --e{715}*/
	SCIPerrorMessage("method LatteBranchRule::scip_free not implemented yet\n");
	return SCIP_OKAY;
}

/** initialization method of branching rule (called after problem was transformed)
 *
 *  @see SCIP_DECL_BRANCHINIT(x) in @ref type_branch.h
 */
SCIP_DECL_BRANCHINIT(LatteBranchRule::scip_init)
{  /*lint --e{715}*/
	SCIPerrorMessage("method LatteBranchRule::scip_init not implemented yet\n");
   return SCIP_OKAY;
}

/** deinitialization method of branching rule (called before transformed problem is freed)
 *
 *  @see SCIP_DECL_BRANCHEXIT(x) in @ref type_branch.h
 */
SCIP_DECL_BRANCHEXIT(LatteBranchRule::scip_exit)
{  /*lint --e{715}*/
	SCIPerrorMessage("method LatteBranchRule::scip_exit not implemented yet\n");
   return SCIP_OKAY;
}

/** solving process initialization method of branching rule (called when branch and bound process is about to begin)
 *
 *  @see SCIP_DECL_BRANCHINITSOL(x) in @ref type_branch.h
 */
SCIP_DECL_BRANCHINITSOL(LatteBranchRule::scip_initsol)
{  /*lint --e{715}*/
	SCIPerrorMessage("method LatteBranchRule::scip_initsol not implemented yet\n");
   return SCIP_OKAY;
}

/** solving process deinitialization method of branching rule (called before branch and bound process data is freed)
 *
 *  @see SCIP_DECL_BRANCHEXITSOL(x) in @ref type_branch.h
 */
SCIP_DECL_BRANCHEXITSOL(LatteBranchRule::scip_exitsol)
{  /*lint --e{715}*/
	SCIPerrorMessage("method LatteBranchRule::scip_exitsol not implemented yet\n");
   return SCIP_OKAY;
}

/** branching execution method for fractional LP solutions
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - branchrule      : the branching rule itself
 *  - allowaddcons    : is the branching rule allowed to add constraints to the current node in order to cut off the
 *                      current solution instead of creating a branching?
 *  - result          : pointer to store the result of the branching call
 *
 *  possible return values for *result (if more than one applies, the first in the list should be used):
 *  - SCIP_CUTOFF     : the current node was detected to be infeasible
 *  - SCIP_CONSADDED  : an additional constraint (e.g. a conflict constraint) was generated; this result code must not be
 *                      returned, if allowaddcons is FALSE
 *  - SCIP_REDUCEDDOM : a domain was reduced that rendered the current LP solution infeasible
 *  - SCIP_SEPARATED  : a cutting plane was generated
 *  - SCIP_BRANCHED   : branching was applied
 *  - SCIP_DIDNOTFIND : the branching rule searched, but did not find a branching
 *  - SCIP_DIDNOTRUN  : the branching rule was skipped
 */
SCIP_DECL_BRANCHEXECLP(LatteBranchRule::scip_execlp)
{  /*lint --e{715}*/
	//SCIPerrorMessage("method LatteBranchRule::scip_execlp not implemented yet\n");
   assert(result != NULL);
   *result = SCIP_DIDNOTRUN;
   return SCIP_OKAY;
}

/** branching execution method for external candidates
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - branchrule      : the branching rule itself
 *  - allowaddcons    : is the branching rule allowed to add constraints to the current node in order to cut off the
 *                      current solution instead of creating a branching?
 *  - result          : pointer to store the result of the branching call
 *
 *  possible return values for *result (if more than one applies, the first in the list should be used):
 *  - SCIP_CUTOFF     : the current node was detected to be infeasible
 *  - SCIP_CONSADDED  : an additional constraint (e.g. a conflict constraint) was generated; this result code must not be
 *                      returned, if allowaddcons is FALSE
 *  - SCIP_REDUCEDDOM : a domain was reduced that rendered the current pseudo solution infeasible
 *  - SCIP_BRANCHED   : branching was applied
 *  - SCIP_DIDNOTFIND : the branching rule searched, but did not find a branching
 *  - SCIP_DIDNOTRUN  : the branching rule was skipped
 */
SCIP_DECL_BRANCHEXECEXT(LatteBranchRule::scip_execext)
{  /* Main branch function */

	SCIP_VAR** orgVars = SCIPgetOrigVars(scip);
	int orgVarsNum = SCIPgetNOrigVars(scip);
	double maxWidth = 0;
	double width;
	double midPoint;
	int maxIndex = -1;


   SCIP_Real brpoint;
   SCIP_NODE* downchild;
   SCIP_NODE* eqchild;
   SCIP_NODE* upchild;

	for(int i = 0; i < orgVarsNum ; ++i)
	{
		if (orgVars[i]->name[0] == 'x')
		{
			width = SCIPvarGetUbLocal(orgVars[i]) - SCIPvarGetLbLocal(orgVars[i]);
			if ( width > maxWidth)
			{
				maxWidth = width;
				maxIndex = i;
				midPoint = (SCIPvarGetUbLocal(orgVars[i]) + SCIPvarGetLbLocal(orgVars[i]))/2.0;
			}
		}
	}

	if ( maxIndex == -1)
	{
		*result = SCIP_DIDNOTFIND;
		return SCIP_OKAY;
	}

	brpoint = SCIPgetBranchingPoint(scip, orgVars[maxIndex], midPoint);


	/* perform the branching */
	SCIP_CALL( SCIPbranchVarVal(scip, orgVars[maxIndex], brpoint, &downchild, &eqchild, &upchild) );

	if( downchild != NULL || eqchild != NULL || upchild != NULL )
	{
	  *result = SCIP_BRANCHED;
	  //cout << "Branched on x" << maxIndex << " at " << midPoint << "\n";
	}
	else
	{
	  /* if there are no children, then variable should have been fixed by SCIPbranchVarVal */
	  assert(SCIPisEQ(scip, SCIPvarGetLbLocal(orgVars[maxIndex]), SCIPvarGetUbLocal(orgVars[maxIndex])));
	  *result = SCIP_REDUCEDDOM;
	}

	return SCIP_OKAY;

}

/** branching execution method for not completely fixed pseudo solutions
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - branchrule      : the branching rule itself
 *  - allowaddcons    : is the branching rule allowed to add constraints to the current node in order to cut off the
 *                      current solution instead of creating a branching?
 *  - result          : pointer to store the result of the branching call
 *
 *  possible return values for *result (if more than one applies, the first in the list should be used):
 *  - SCIP_CUTOFF     : the current node was detected to be infeasible
 *  - SCIP_CONSADDED  : an additional constraint (e.g. a conflict constraint) was generated; this result code must not be
 *                      returned, if allowaddcons is FALSE
 *  - SCIP_REDUCEDDOM : a domain was reduced that rendered the current pseudo solution infeasible
 *  - SCIP_BRANCHED   : branching was applied
 *  - SCIP_DIDNOTFIND : the branching rule searched, but did not find a branching
 *  - SCIP_DIDNOTRUN  : the branching rule was skipped
 */
SCIP_DECL_BRANCHEXECPS(LatteBranchRule::scip_execps)
{  /*lint --e{715}*/
	//SCIPerrorMessage("method LatteBranchRule::scip_execps not implemented yet\n");
   assert(result != NULL);
   *result = SCIP_DIDNOTRUN;
   return SCIP_OKAY;
}


