/*
 * LatteIntegrationRelaxor.cpp
 *
 *  Created on: Apr 28, 2015
 *      Author: bedutra
 */

#include "LatteIntegrationRelaxor.h"
#include "scip/primal.h"
#include "scip/struct_mem.h"
#include <fstream>


/*
 * Local methods
 */



LatteIntegrationRelaxor::LatteIntegrationRelaxor(SCIP* scip, char * fileName): ObjRelax(scip, "LatteIntRelaxor", "Latte Integration Relaxor", 1000, 1)
{
	pipFileName = fileName;
}

LatteIntegrationRelaxor::~LatteIntegrationRelaxor() {
	//do nothing. scip_free has job of freeing stuff
}


/** destructor of relaxator to free user data (called when SCIP is exiting)	*/
SCIP_DECL_RELAXFREE(LatteIntegrationRelaxor::scip_free)
{
	SCIPerrorMessage("method LatteIntegrationRelaxor::scip_free not implemented yet\n");

	cout << "Latte Integration Relaxor Stats" << endl;



	return SCIP_OKAY;
}


/** initialization method of relaxator (called after problem was transformed)*/
SCIP_DECL_RELAXINIT(LatteIntegrationRelaxor::scip_init)
{


	if (pipFileName == NULL) {
		cerr << "polynomial file name is missing" << endl;
		THROW_LATTE(LattException::ue_FileNameMissing);
	}

	//get polynomial objective
	ifstream polynomialFile(pipFileName);
	string polyStr;
	getline(polynomialFile, polyStr);
	getline(polynomialFile, polyStr);
	cout << "poly str: " << polyStr.c_str() << endl;
	loadMonomials(originalPolynomial, polyStr.c_str());
	polynomialFile.close();


	int k = 3;



	return SCIP_OKAY;
}


/** deinitialization method of relaxator (called before transformed problem is freed)*/
SCIP_DECL_RELAXEXIT(LatteIntegrationRelaxor::scip_exit)
{
	//SCIPerrorMessage("method LatteSummationRelaxor::scip_exit not implemented yet\n");
   return SCIP_OKAY;
}


/** solving process initialization method of relaxator (called when branch and bound process is about to begin)*/
SCIP_DECL_RELAXINITSOL(LatteIntegrationRelaxor::scip_initsol)
{
	//SCIPerrorMessage("method LatteSummationRelaxor::scip_initsol not implemented yet\n");
   return SCIP_OKAY;
}


/** solving process deinitialization method of relaxator (called before branch and bound process data is freed)*/
SCIP_DECL_RELAXEXITSOL(LatteIntegrationRelaxor::scip_exitsol)
{
	//SCIPerrorMessage("method LatteSummationRelaxor::scip_exitsol not implemented yet\n");
   return SCIP_OKAY;
}


/** execution method of relaxator*/
/** execution method of relaxator
 *
 *  The method is called in the node processing loop. It solves the current subproblem's relaxation.
 *  Like the LP relaxation, the relaxator should only operate on COLUMN variables.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - relax           : the relaxator itself
 *  - lowerbound      : pointer to store a lowerbound for the current node
 *  - result          : pointer to store the result of the relaxation call
 *
 *  possible return values for *result (if more than one applies, the first in the list should be used):
 *  - SCIP_CUTOFF     : the node is infeasible in the variable's bounds and can be cut off
 *  - SCIP_CONSADDED  : an additional constraint was generated, and the relaxator should not be called again on the
 *                      same relaxation
 *  - SCIP_REDUCEDDOM : a variable's domain was reduced, and the relaxator should not be called again on the same
 *                      relaxation
 *  - SCIP_SEPARATED  : a cutting plane was generated, and the relaxator should not be called again on the same relaxation
 *  - SCIP_SUCCESS    : the relaxator solved the relaxation and should not be called again on the same relaxation
 *  - SCIP_SUSPENDED  : the relaxator interrupted its solving process to wait for additional input (e.g. cutting
 *                      planes); however, it is able to continue the solving in order to improve the dual bound
 *  - SCIP_DIDNOTRUN  : the relaxator was skipped
 */
//SCIP_DECL_RELAXEXEC(x) SCIP_RETCODE x (SCIP* scip, SCIP_RELAX* relax, SCIP_Real* lowerbound, SCIP_RESULT* result)
SCIP_DECL_RELAXEXEC(LatteIntegrationRelaxor::scip_exec)
{

	bool isImprovement = false;
	SCIP_VAR** orgVars = SCIPgetOrigVars(scip);
	int orgVarsNum = SCIPgetNOrigVars(scip);
	vec_RR lowerBound, upperBound;
	lowerBound.SetLength(orgVarsNum - 1);
	upperBound.SetLength(orgVarsNum - 1);


	for(int i = 0; i < orgVarsNum - 1; ++i)
	{
		lowerBound[i] = SCIPvarGetLbLocal(orgVars[i]);
		upperBound[i] = SCIPvarGetUbLocal(orgVars[i]);

	}

	//> Note that SCIPgetLowerbound() and SCIPgetUpperbound() do the same as
	//> SCIPgetDualbound() and SCIPgetPrimalbound(), but with respect to the
	//> transformed problem space.

	box1.setBounds(lowerBound, upperBound);



	//********************************************
	//New local lower bound in transformed problem
	//********************************************
	// sense minimize -poly
	if (SCIPgetLocalLowerbound(scip) < -to_double(box1.U))
	{
		numBasicLocal++;
		isImprovement = true;
	}

	SCIPupdateLocalLowerbound(scip, -to_double(box1.U));


	//*********************************************
	//New global upper bound in transformed problem
	//*********************************************

	SCIP_Real newUpperBound  = -to_double(box1.L);


	//note that there are two SCIPprimalSetUpperbound on the online doc, but only one with this signature.
	if ( newUpperBound < scip->primal->upperbound)
	{
		SCIPprimalSetUpperbound(scip->primal, scip->mem->probmem, scip->set, scip->stat, scip->eventqueue, scip->transprob, scip->tree, scip->lp, newUpperBound);
		numBasicGlobal++;
		isImprovement = true;
	}


	//*********************************************************************
	//eval some points. Gives new global upper bound in transformed problem
	//*********************************************************************
	//TODO

	//*******************
	//end of cheap tricks
	//*******************
	if ( box1.N > 500000 || box1.N < 5000  || isImprovement)
		return SCIP_OKAY;

	if ( numFailedToImprove < 0)
	{
		numFailedToImprove++;
		return SCIP_OKAY;
	}

	//box1.printNumberOfPoints();
	//SCIPretransformObj(scip, 2);
	//SCIPsetObjlimit(scip, 2);



	box1.decomposePoly(BoxOptimization::naturalSummation);
	box1.findSPolynomial(BoxOptimization::naturalSummation, lowerBound, upperBound);
	box1.findRange(10);


	box2.setBounds(lowerBound, upperBound);
	box2.decomposePoly(BoxOptimization::naturalSummation);
	box2.findSPolynomial(BoxOptimization::naturalSummation, lowerBound, upperBound);
	box2.findRange(10);



	if (SCIPgetLocalLowerbound(scip) < -to_double(box2.maximumUpperbound()) )
	{	numKnormLocal++; isImprovement = true;}
	if (SCIPgetLocalLowerbound(scip) < -to_double(box2.U) )
	{	numKnormRangeLocal++; isImprovement = true;}
	if ( scip->primal->upperbound >  -to_double(box2.maximumLowerBound()) )
	{	numKnormGlobal++; isImprovement = true;}
	if ( scip->primal->upperbound >  -to_double(box2.L) )
	{	numKnormRangeGlobal++; isImprovement = true;}
	if ( scip->primal->upperbound >  -to_double(box2.L + box2.currentMap.eval(-box2.L) / box1.currentMap.eval(-box2.L)) )
	{	numRatioGlobal++; isImprovement = true;}
	numLatteCalled++;

	if ( isImprovement == false)
		numFailedToImprove++;
	if (numFailedToImprove >= 1)
		numFailedToImprove = -1000;


	SCIPupdateLocalLowerbound(scip, min(-to_double(box2.maximumUpperbound()), -to_double(box2.U))); //Maybe not better than -box1.U, but the scip function will take the best of the two.


	newUpperBound = -to_double(box2.L + box2.currentMap.eval(-box2.L) / box1.currentMap.eval(-box2.L));
	newUpperBound = min(newUpperBound, -to_double(box2.maximumLowerBound()) );
	newUpperBound = min(newUpperBound, -to_double(box2.L));
	if ( newUpperBound < scip->primal->upperbound)
	{
		SCIPprimalSetUpperbound(scip->primal, scip->mem->probmem, scip->set, scip->stat, scip->eventqueue, scip->transprob, scip->tree, scip->lp, newUpperBound);
		cout << "new global bound by summation" << endl;
	}





	*result = SCIP_SUCCESS;
	return SCIP_OKAY;
}



