/*
 * LatteBranchRule.h
 *
 *  Created on: May 8, 2015
 *      Author: bedutra
 */

#ifndef LATTEBRANCHRULE_H_
#define LATTEBRANCHRULE_H_


#include "objscip/objbranchrule.h"

class LatteBranchRule:  public  scip::ObjBranchrule
{
public:
	LatteBranchRule(SCIP* scip);
	virtual ~LatteBranchRule();



	   /** destructor of branching rule to free user data (called when SCIP is exiting)
	    *
	    *  @see SCIP_DECL_BRANCHFREE(x) in @ref type_branch.h
	    */
	   virtual SCIP_DECL_BRANCHFREE(scip_free);

	   /** initialization method of branching rule (called after problem was transformed)
	    *
	    *  @see SCIP_DECL_BRANCHINIT(x) in @ref type_branch.h
	    */
	   virtual SCIP_DECL_BRANCHINIT(scip_init);

	   /** deinitialization method of branching rule (called before transformed problem is freed)
	    *
	    *  @see SCIP_DECL_BRANCHEXIT(x) in @ref type_branch.h
	    */
	   virtual SCIP_DECL_BRANCHEXIT(scip_exit);

	   /** solving process initialization method of branching rule (called when branch and bound process is about to begin)
	    *
	    *  @see SCIP_DECL_BRANCHINITSOL(x) in @ref type_branch.h
	    */
	   virtual SCIP_DECL_BRANCHINITSOL(scip_initsol);

	   /** solving process deinitialization method of branching rule (called before branch and bound process data is freed)
	    *
	    *  @see SCIP_DECL_BRANCHEXITSOL(x) in @ref type_branch.h
	    */
	   virtual SCIP_DECL_BRANCHEXITSOL(scip_exitsol);

	   /** branching execution method for fractional LP solutions
	    *
	    *  @see SCIP_DECL_BRANCHEXECLP(x) in @ref type_branch.h
	    */
	   virtual SCIP_DECL_BRANCHEXECLP(scip_execlp);

	   /** branching execution method for external candidates
	    *
	    *  @see SCIP_DECL_BRANCHEXECEXT(x) in @ref type_branch.h
	    */
	   virtual SCIP_DECL_BRANCHEXECEXT(scip_execext);

	   /** branching execution method for not completely fixed pseudo solutions
	    *
	    *  @see SCIP_DECL_BRANCHEXECPS(x) in @ref type_branch.h
	    */
	   virtual SCIP_DECL_BRANCHEXECPS(scip_execps);

};


SCIP_RETCODE SCIPincludeBranchrule_LatteBranchRule(
   SCIP*                 scip                /**< SCIP data structure */
);

#endif /* LATTEBRANCHRULE_H_ */
