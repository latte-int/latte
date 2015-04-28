/*
 * LatteSummationRelaxor.h
 *
 *  Created on: Nov 4, 2014
 *      Author: bedutra
 */

#ifndef LATTESUMMATIONRELAXOR_H_
#define LATTESUMMATIONRELAXOR_H_

#include <string>

#include "objscip/objrelax.h"

#include "banner.h"
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


class LatteSummationRelaxor: public scip::ObjRelax
{
private:
	char * pipFileName;
	monomialSum originalPolynomial;
	BoxOptimization box1;
	BoxOptimization box2;


	ZZ numFailedToImprove;


	ZZ numBasicLocal;
	ZZ numBasicGlobal;
	ZZ numRandPointImprove;
	ZZ numKnormLocal;
	ZZ numKnormRangeLocal;
	ZZ numKnormGlobal;
	ZZ numKnormRangeGlobal;
	ZZ numRatioGlobal;
	ZZ numCalled;
	ZZ numLatteCalled;
public:


	LatteSummationRelaxor(SCIP* scip, char * fileName);   /**< SCIP data structure */
	virtual ~LatteSummationRelaxor();


	/** destructor of relaxator to free user data (called when SCIP is exiting)	*/
	virtual SCIP_DECL_RELAXFREE(scip_free);

	/** initialization method of relaxator (called after problem was transformed)*/
	virtual SCIP_DECL_RELAXINIT(scip_init);

	/** deinitialization method of relaxator (called before transformed problem is freed)*/
	virtual SCIP_DECL_RELAXEXIT(scip_exit);

	/** solving process initialization method of relaxator (called when branch and bound process is about to begin)*/
	virtual SCIP_DECL_RELAXINITSOL(scip_initsol);

	/** solving process deinitialization method of relaxator (called before branch and bound process data is freed)*/
	virtual SCIP_DECL_RELAXEXITSOL(scip_exitsol);

	/** execution method of relaxator*/
	virtual SCIP_DECL_RELAXEXEC(scip_exec);

};






#endif /* LATTESUMMATIONRELAXOR_H_ */
