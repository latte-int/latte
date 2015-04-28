/*
 * LatteSummationRelaxor.cpp
 *
 *  Created on: Nov 4, 2014
 *      Author: bedutra
 */

#include "LatteSummationRelaxor.h"
#include "scip/primal.h"
#include "scip/struct_mem.h"
#include <fstream>


/*
 * Local methods
 */



LatteSummationRelaxor::LatteSummationRelaxor(SCIP* scip, char * fileName): ObjRelax(scip, "LatteSumRelaxor", "Latte Summation Relaxor", 1000, 1)
{
	pipFileName = fileName;
}

LatteSummationRelaxor::~LatteSummationRelaxor() {
	//do nothing. scip_free has job of freeing stuff
}


/** destructor of relaxator to free user data (called when SCIP is exiting)	*/
SCIP_DECL_RELAXFREE(LatteSummationRelaxor::scip_free)
{
	SCIPerrorMessage("method LatteSummationRelaxor::scip_free not implemented yet\n");

	cout << "Latte Summation Relaxor Stats" << endl;


	cout << "     BASIC BOUNDS" << endl;
	cout << "       numBasicLocal " << numBasicLocal << endl;
	cout << "       numBasicGlobal " << numBasicGlobal << endl;
	cout << "       numRandPointImprove " << numRandPointImprove << endl;
	cout << "     K NORM RANGE" << endl;
	cout << "       numKnormRangeLocal " << numKnormRangeLocal << endl;
	cout << "       numKnormRangeGlobal " << numKnormRangeGlobal << endl;
	cout << "     K NORM BOUNDS " << endl;
	cout << "       numKnormLocal " << numKnormLocal << endl;
	cout << "       numKnormGlobal " << numKnormGlobal << endl;
	cout << "       numRatioGlobal " << numRatioGlobal << endl;
	cout << "     STATS" << endl;
	cout << "       numCalled " << numCalled << endl;
	cout << "       numLatteCalled " << numLatteCalled << endl;

	return SCIP_OKAY;
}


/** initialization method of relaxator (called after problem was transformed)*/
SCIP_DECL_RELAXINIT(LatteSummationRelaxor::scip_init)
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

	box1.setPolynomial(originalPolynomial);
	box1.setPower(k);

	box2 = box1;

	//box2.setPolynomial(originalPolynomial);
	box2.setPower(k+1);

	return SCIP_OKAY;
}


/** deinitialization method of relaxator (called before transformed problem is freed)*/
SCIP_DECL_RELAXEXIT(LatteSummationRelaxor::scip_exit)
{
	//SCIPerrorMessage("method LatteSummationRelaxor::scip_exit not implemented yet\n");
   return SCIP_OKAY;
}


/** solving process initialization method of relaxator (called when branch and bound process is about to begin)*/
SCIP_DECL_RELAXINITSOL(LatteSummationRelaxor::scip_initsol)
{
	//SCIPerrorMessage("method LatteSummationRelaxor::scip_initsol not implemented yet\n");
   return SCIP_OKAY;
}


/** solving process deinitialization method of relaxator (called before branch and bound process data is freed)*/
SCIP_DECL_RELAXEXITSOL(LatteSummationRelaxor::scip_exitsol)
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
SCIP_DECL_RELAXEXEC(LatteSummationRelaxor::scip_exec)
{

	bool isImprovement = false;
	SCIP_VAR** orgVars = SCIPgetOrigVars(scip);
	int orgVarsNum = SCIPgetNOrigVars(scip);
	vec_ZZ lowerBound, upperBound;
	lowerBound.SetLength(orgVarsNum - 1);
	upperBound.SetLength(orgVarsNum - 1);

	numCalled++;
	//NOTES:
	//outfile << "node" << SCIPgetNNodes(scip) << " depth " << SCIPgetDepth(scip)  <<endl;
	//for(int i = 0; i < orgVarsNum; ++i)
	//{
	//	cout << "node" << SCIPgetNNodes(scip) << " " << SCIPvarGetLbLocal(orgVars[i]) << "<= " << orgVars[i]->name << "<= " << SCIPvarGetUbLocal(orgVars[i]) << endl;
	//}


	for(int i = 0; i < orgVarsNum - 1; ++i)
	{
		lowerBound[i] = SCIPvarGetLbLocal(orgVars[i]);
		upperBound[i] = SCIPvarGetUbLocal(orgVars[i]);

	}

	//> Note that SCIPgetLowerbound() and SCIPgetUpperbound() do the same as
	//> SCIPgetDualbound() and SCIPgetPrimalbound(), but with respect to the
	//> transformed problem space.

	box1.setBounds(lowerBound, upperBound);


	//cout << SCIPgetNNodes(scip) << " min scip  " << SCIPgetLocalLowerbound(scip) << " <= " << SCIPgetUpperbound(scip) << endl;
	//cout << SCIPgetNNodes(scip) << " min my b  " << -to_double(box1.U)             << " <= " << -to_double(box1.L) << endl;


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


	//cout << "SCIPgetObjlimit " << SCIPgetObjlimit(scip) << endl;
	//cout << "SCIPretransformObj " << SCIPretransformObj(scip,newUpperBound) << endl;
	//SCIPsetObjlimit(scip, SCIPretransformObj(scip,newUpperBound));
    //SCIP_CALL( SCIPprimalUpdateObjlimit(scip->primal, scip->mem->probmem, scip->set, scip->stat, scip->eventqueue, scip->transprob, scip->origprob, scip->tree, scip->lp) );

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
	/*
	vec_ZZ point;
	point.SetLength(orgVarsNum-1);
	for(int i = 0; i < orgVarsNum-1; ++i)
	{
		point[i] = (upperBound[i]+lowerBound[i])/2;
	}

	SCIP_SOL* newsol;
	SCIPcreateOrigSol(scip, &newsol,NULL);
	for(int i = 0; i < orgVarsNum-1; ++i)
	{
		SCIPsetSolVal(scip, newsol, orgVars[i], to_double(point[i]));
	}

	SCIPsetSolVal(scip, newsol, orgVars[orgVarsNum-1], to_double(box1.sampleLowerBound(originalPolynomial,point)) );
	unsigned int success;
	//for(int i = 0; i < orgVarsNum; ++i)
	//	cout << "var" << i << " " << SCIPgetSolVal(scip, newsol, orgVars[i]) << endl;

	SCIPtrySol(scip, newsol, TRUE, TRUE, TRUE, FALSE, &success);
	if (success)
	{
		numRandPointImprove++;
		isImprovement = true;
	}
	*/

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

	box1.printNumberOfPoints();
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

	// |f|_{k+1}/|f|_{k} \leq max f -L

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


