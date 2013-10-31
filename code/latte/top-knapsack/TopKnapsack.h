/*
 * TopKnapsack.h
 *
 *  Created on: Apr 23, 2013
 *      Author: bedutra
 */

#ifndef TOPKNAPSACK_H_H
#define TOPKNAPSACK_H_H

#include "latte_ntl.h"
#include "LattException.h"
#include "cone.h"
#include "integration/burstTrie.h"
#include "integration/PolyTrie.h"
#include "integration/GeneralMonomialSum.h"
#include "PeriodicFunction.h"
#include <vector>

using namespace std;

class TopKnapsack;

class MobiusPair{
public:
	ZZ gcd;
	ZZ mu;
	bool mobiusValid;

	MobiusPair();
	MobiusPair(const ZZ& g, const ZZ& m);
};

class MobiusList {
private:
	void computeMobius(int i);
public:
	vector<MobiusPair> list;
	MobiusList();
	void insertGCD(const ZZ& v);
	void computeMobius();
	void print() const;
};


class MobiusSeriesList: public MobiusList
{
public:
	vector<GeneralMonomialSum<PeriodicFunction, int> *> unweightedSeries; //[i] is a series corresponding to list[i]. Still needs to by build by list[i].mu and -1^i/i!
	void computeMobius();
};

class BernoulliFirstKind
{
private:
	vector<RationalNTL> B;
public:
	void setBernoulli(int k); //find first k Bernoulli numbers of the first kind.
	const RationalNTL& operator[](int i) const;
};




class TopKnapsack {
public:
	vec_ZZ alpha;
	int N;
	int order; //k as in N-k
	int largestSubsetComputed;
	MobiusSeriesList gcds;
	BernoulliFirstKind bernoulli;
	vector<PeriodicFunction> coeffsNminusk;


	void everyGCD(int k);


	void E(const int fIndex);
	void findLatticeBasis(mat_ZZ & latticeBasis, const vector<ZZ> & fnDivAlpha, const ZZ & f) const;
	void findLatticeBasisInv(mat_ZZ & invLatticeBasisScaled, ZZ & invLatticeBasisD, mat_ZZ & invLatticeBasis,const mat_ZZ & latticeBasis) const;
	void findVertex(vec_ZZ & tVertex, const ZZ &f, const vector<ZZ> &fnDivAlpha) const;
	listCone*  findUnimodularCones(const mat_ZZ & invLatticeBasis)const;
	void findResidue(GeneralMonomialSum<PeriodicFunction, int> & fSeries, const vec_ZZ & tVertex, const listCone *uniCones,
			const mat_ZZ & latticeBasis, const mat_ZZ & invLatticeBasis,
			const ZZ &invLatticeBasisD, const vector<ZZ> &fnDivAlpha, const vector<ZZ> & fDivAlpha);

	void expandNonperiodicPart(GeneralMonomialSum<PeriodicFunction, int> &a, const vector<ZZ>  & fDivAlpha);
	void expandPeriodicPart(ZZ & bottomCoeffPeriodicPart, GeneralMonomialSum<PeriodicFunction, int> & a, const int numPoles, const vector<ZZ> & expa, const vector<ZZ> & expe);
	void expandExponentialPart(GeneralMonomialSum<PeriodicFunction, int> & exExpansion, const int numPoles, const vec_ZZ &a, const vec_ZZ & e, const vector<RationalNTL> & f);

	void expandF1Case(GeneralMonomialSum<PeriodicFunction, int> & expansion);

	void packageAnswer();
	void printAnswer(ostream & out);
/*
	void expandOneTerm(monomialSum & oneExpansion, const ZZ & expx, const ZZ & expe);
	void expandOneNonPerturbedTerm(monomialSum & oneExpansion, const ZZ & a);
	void expandPeriodicExponential(GeneralMonomialSum<PeriodicFunction, int> & pxeProduct,
			const mat_ZZ &B, const mat_ZZ & invB,
			const ZZ & invBd, const listCone * oneCone,
			const vec_ZZ & s, const vec_ZZ & alpha, const vec_ZZ & beta);

	void expandPerturbationTerms(GeneralMonomialSum<PeriodicFunction, int> & ans,
			const GeneralMonomialSum<PeriodicFunction, int> & p1, const monomialSum & p2);
*/
	static ZZ binomial(int n, int k);
	TopKnapsack();
	~TopKnapsack();

	void set(const vec_ZZ& list);

	void coeff_NminusK(int k);

	static void printMatrix(const mat_ZZ &A);//for debugging.

};



#endif /* TOPKNAPSACK_H_ */
