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
	vector<MobiusPair> list;
	void computeMobius(int i);
public:
	MobiusList();
	void insertGCD(const ZZ& v);
	void computeMobius();
	void print() const;

friend class TopKnapsack;
};


class BernoulliSecondKind
{
private:
	vector<RationalNTL> B;
public:
	void setBernoulli(int k); //find first k Bernoulli numbers of the second kind.
	const RationalNTL& operator[](int i) const;
};

class TopKnapsack {
private:
	vec_ZZ alpha;
	int N;
	int order; //k as in N-k
	int largestSubsetComputed;
	MobiusList gcds;
	BernoulliSecondKind bernoulli;

	void everyGCD(int k);


	void E(ZZ f);
	void findLatticeBasis(mat_ZZ & latticeBasis, const vector<ZZ> & fnDivAlpha, const ZZ & f) const;
	void findVertex(vec_ZZ & tVertex, const ZZ &f, const vector<ZZ> &fnDivAlpha) const;
	listCone*  findUnimodularCones(const mat_ZZ & invLatticeBasis)const;
	void findResidue(const vec_ZZ & tVertex, const listCone *uniCones,
			const mat_ZZ & latticeBasis, const mat_ZZ & invLatticeBasis,
			const ZZ &invLatticeBasisD, const vector<ZZ> &fnDivAlpha);
	void expandOneTerm(monomialSum & oneExpansion, const ZZ & expx, const ZZ & expe);
	void expandPeriodicExponential(GeneralMonomialSum<PeriodicFunction, int> & pxeProduct,
			const mat_ZZ &B, const mat_ZZ & invB,
			const ZZ & invBd, const listCone * oneCone,
			const vec_ZZ & s, const vec_ZZ & alpha, const vec_ZZ & beta);

	void expandPerturbationTerms(GeneralMonomialSum<PeriodicFunction, int> & ans,
			const GeneralMonomialSum<PeriodicFunction, int> & p1, const monomialSum & p2);
public:
	static ZZ binomial(int n, int k);
	TopKnapsack();
	~TopKnapsack();

	void set(const vec_ZZ& list);
	void findAll();
	void coeff_NminusK(int k);

	static void printMatrix(const mat_ZZ &A);//for debugging.

};



#endif /* TOPKNAPSACK_H_ */
