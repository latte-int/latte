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


	int k = 4;

	box1.setPolynomial(originalPolynomial);
	box1.setPower(k);

	box2 = box1;

	//box2.setPolynomial(originalPolynomial);
	box2.setPower(k+1);



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
	int numLatteVars = 0;
	vec_RR lowerBound, upperBound;
	static vec_RR prevLB, prevUB;
	static std::vector<int> varMap;


	*result = SCIP_SUCCESS;
	if (SCIPgetLocalLowerbound(scip) + 0.001> scip->primal->upperbound)
	{
		*result = SCIP_DIDNOTRUN;
		return SCIP_OKAY;
	}
	numCalled++;

	for(int i = 0; i < orgVarsNum; ++i)
		if (orgVars[i]->name[0] == 'x')
			numLatteVars++;

	lowerBound.SetLength(numLatteVars);
	upperBound.SetLength(numLatteVars);

	if ( !varMap.size())
		varMap.resize(numLatteVars);

	int var;
	RR boxLen;
	for(int i = 0; i < orgVarsNum ; ++i)
	{
		if (orgVars[i]->name[0] == 'x')
		{
			var = atoi(orgVars[i]->name + 1) -1;
			//cout << "var = " << var << endl;
			lowerBound[var] = SCIPvarGetLbLocal(orgVars[i]);
			upperBound[var] = SCIPvarGetUbLocal(orgVars[i]);
			boxLen += upperBound[var] - lowerBound[var];
			varMap[var] = i;
		}
		//cout << "node" << SCIPgetNNodes(scip) << " index "<< i << " " << SCIPvarGetLbLocal(orgVars[i]) << "<= " << orgVars[i]->name << "<= " << SCIPvarGetUbLocal(orgVars[i])  << endl;
	}




	bool newProblem = false;
	if (numCalled <= 1)
	{
		prevLB = lowerBound;
		prevUB = upperBound;
		newProblem = true;
	}
	else
	{
		for (int i = 0; i < numLatteVars && !newProblem; ++i)
		{
			if ( lowerBound[i] != prevLB[i] || upperBound[i] != prevUB[i])
			{
				//cout << SCIPgetNNodes(scip) << "old: " << prevLB[i] << " <=x[" << i << "]<= " << prevLB[i] << " new: "   << lowerBound[i] << " <=x[" << i << "]<= " << upperBound[i] << "\n";
				newProblem = true;
			}
		}

		prevLB = lowerBound;
		prevUB = upperBound;
	}

	if (! newProblem )
	{
		//cout << " same box on " << SCIPgetNNodes(scip) << "; ";
		//cout << "attempted to branch on x" << maxIndex << endl;
		//SCIPaddExternBranchCand(scip, orgVars[varMap[maxIndex]], -1000000, SCIPvarGetAvgSol(orgVars[varMap[maxIndex]]));
		 //  SCIP*                 scip,               /**< SCIP data structure */
		  // SCIP_VAR*             var,                /**< variable to insert */
		  // SCIP_Real             score,              /**< score of external candidate, e.g. infeasibility */
		  // SCIP_Real             solval              /**< value of the variable in the current solution */
		  // )
		numCalled--;
		return SCIP_OKAY;
	}


	box1.setBounds(lowerBound, upperBound); //k

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
	if ( isImprovement )//|| boxLen > 10)
		return SCIP_OKAY;


	if ( numFailedToImprove < 0)
	{
		numFailedToImprove++;
		return SCIP_OKAY;
	}


	cout << "**before scip bound " <<  SCIPgetLocalLowerbound(scip) << " <=min -f <= " <<scip->primal->upperbound  <<endl;

	box1.findSPolynomial(BoxOptimizationContinuous::naturalSummation, lowerBound, upperBound);


	box2.setBounds(lowerBound, upperBound);
	cout << lowerBound <<endl;
	cout <<upperBound <<endl;
	box2.printStats();
	box2.lipschitz = 10000  ;
	box2.findSPolynomial(BoxOptimizationContinuous::naturalSummation, lowerBound, upperBound);
	box2.printSpolynomial();

	//see if we can get a better range bound for f(x).
	box2.U = min(box2.U, -to_RR(SCIPgetLocalLowerbound(scip)));
	box1.U = min(box1.U, -to_RR(SCIPgetLocalLowerbound(scip)));
	//cout << "range update is off b/c lipschitz is invalid" <<endl;

	box2.setFmaxLowerbound(-SCIPgetLocalLowerbound(scip));
	//bool fr = box2.findRange(100);
	//if (fr)
	//	cout << fr << "=range update improved box bound <-----------------------------------------------------------------------------------" << endl;





	/*
	int d = box2.originalPolynomial.varCount -1;
	RR p1, p2;
	p1 = d;
	p1 /= d+box2.currentPower;
	p2 = box2.currentPower;
	p2 /= d+box2.currentPower;

	//p1 = d/(d+k)
	//p2 = k/(d+k)

	RR f1, f2, f3,newL;

	int sign = 1;
	if (box2.currentPower % 2 )
		sign = -1;

	cout << sign << " cur pow" << box2.currentPower << endl;
	cout << "box2.evalSpoly(-box2.U)*|sign|=" <<box2.evalSpoly(-box2.U) << endl;

	f1 = pow(box2.evalSpoly(-box2.U)*sign/box2.V, inv(to_RR(d+box2.currentPower)));
	f2 = pow( (box2.M*box2.lipschitz * (d+box2.currentPower)) * inv(to_RR(d)), p1);
	f3 = pow( to_RR(d+box2.currentPower)/to_RR(box2.currentPower), p2);


	newL = box2.U - f1*f2*f3;

	cout << "wtf2" << endl;
	*/

	/*
	{
		int d = box2.originalPolynomial.varCount -1;
		RR p1, p2;
		p1 = d;
		p1 /= d+box2.currentPower;
		p2 = box2.currentPower;
		p2 /= d+box2.currentPower;

		//p1 = d/(d+k)
		//p2 = k/(d+k)

		RR f1, f2, f3;

		f1 = pow(box2.evalSpoly(-box2.L)/box2.V, inv(to_RR(d+box2.currentPower)));
		f2 = pow( (box2.M*box2.lipschitz * (d+box2.currentPower)) * inv(to_RR(d)), p1);
		f3 = pow( to_RR(d+box2.currentPower)/to_RR(box2.currentPower), p2);

		cout << "s-poly at " << -box2.L << " is " << box2.evalSpoly(-box2.L) << endl;
		cout << f1 << " * " << f2 << " *  " << f3 << " = " << f1*f2*f3 <<endl;
		cout << "u()="  << box2.L + f1*f2*f3 << endl;
	}
	*/



	if (SCIPgetLocalLowerbound(scip) < -to_double(box2.maximumUpperbound()) )
	{	numKnormLocal++; isImprovement = true;}
	if (SCIPgetLocalLowerbound(scip) < -to_double(box2.U) )
	{	numKnormRangeLocal++; isImprovement = true;}
	if ( scip->primal->upperbound >  -to_double(box2.maximumLowerBound()) )
	{	numKnormGlobal++; isImprovement = true;}
	if ( scip->primal->upperbound >  -to_double(box2.L) )
	{	numKnormRangeGlobal++; isImprovement = true;}
	if ( scip->primal->upperbound >  -to_double(box2.L + box2.evalSpoly(-box2.L) / box1.evalSpoly(-box2.L)) )
	{	numRatioGlobal++; isImprovement = true;}
	numLatteCalled++;


	//cout << "node " << SCIPgetNNodes(scip) << " k " << box1.currentPower <<  endl;

	//cout << "k-nomr " <<  -to_double(box2.maximumUpperbound()) << " <=min -f <= " << -to_double(box2.maximumLowerBound()) << endl;
	//cout <<  box2.maximumLowerBound() << " <= max f <= " << box2.maximumUpperbound() << endl;
	//cout << "Ratio " <<  " min -f <= "<< -to_double(box2.L + box2.evalSpoly(-box2.L) / box1.evalSpoly(-box2.L)) << endl;
	//cout << "      " << "L=" <<box2.L << " b2=" << box2.evalSpoly(-box2.L) << " b1=" << box1.evalSpoly(-box2.L) << endl;



	cout << "box len " << boxLen << endl;


	numFailedToImprove = -50;
//	if ( isImprovement == false)
//		numFailedToImprove++;
//	if (numFailedToImprove >= 1)
//		numFailedToImprove = -1000;




	cout << "lower bound s" << SCIPgetLocalLowerbound(scip) << ", -U=" << -to_double(box2.U) << ", -u()=" << -to_double(box2.maximumUpperbound()) << endl;
	SCIPupdateLocalLowerbound(scip, max(-to_double(box2.maximumUpperbound()), -to_double(box2.U))); //Maybe not better than -box1.U, but the scip function will take the best of the two.


	newUpperBound = -to_double(box2.L + box2.evalSpoly(-box2.L) / box1.evalSpoly(-box2.L));
	newUpperBound = min(newUpperBound, -to_double(box2.maximumLowerBound()) );
	newUpperBound = min(newUpperBound, -to_double(box2.L));

	cout << "*latte bound " <<  -to_double(box2.maximumUpperbound()) << " <=min -f <= " << newUpperBound <<endl;

	//cout << "new global bound by summation: new= " << newUpperBound << " old=" << scip->primal->upperbound << endl;
	if ( newUpperBound < scip->primal->upperbound)
	{
		//cout << "new global bound by summation: new= " << newUpperBound << " old=" << scip->primal->upperbound << endl;
		SCIPprimalSetUpperbound(scip->primal, scip->mem->probmem, scip->set, scip->stat, scip->eventqueue, scip->transprob, scip->tree, scip->lp, newUpperBound);

	}


	cout << "**after scip bound " <<  SCIPgetLocalLowerbound(scip) << " <=min -f <= " <<scip->primal->upperbound  <<endl;
	cout << "isImprovement" << isImprovement << endl;
	if ( isImprovement)
		cout << "#<---------------------------------------------------------------------------------------------------------------------------------------\n";

	/*
	cout << "isImprov = " << isImprovement << endl;
	if (isImprovement == false)
	{
		double maxWidth = 0;
		double midPoint, width;
		int maxIndex;
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

		brpoint = SCIPgetBranchingPoint(scip, orgVars[maxIndex], midPoint);


		// perform the branching
		SCIP_CALL( SCIPbranchVarVal(scip, orgVars[maxIndex], brpoint, &downchild, &eqchild, &upchild) );

		if( downchild != NULL || eqchild != NULL || upchild != NULL )
		{
		 // *result = SCIP_BRANCHED;
		  cout << "Branched on x" << maxIndex << " at " << midPoint << endl;
		}
		else
		{
		  // if there are no children, then variable should have been fixed by SCIPbranchVarVal
		  assert(SCIPisEQ(scip, SCIPvarGetLbLocal(orgVars[maxIndex]), SCIPvarGetUbLocal(orgVars[maxIndex])));
		  cout << "var " << maxIndex << " fixed" << endl;
		  // *result = SCIP_REDUCEDDOM;
		}
	}
	*/




	return SCIP_OKAY;
}



