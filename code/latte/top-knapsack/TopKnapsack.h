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


/**
 * This class computes the top k terms of the Ehrhart polynomial for knapsack polytopes.
 * This program implements the algorithms in "Coefficients of the Denumerant" by Velleda Baldoni, Nicole Berline, Jes\'us A. De Loera, Brandon E. Dutra, Matthias K\"oppe, and Mich\`ele Vergne
 */

using namespace std;

class TopKnapsack;

/**
 * Holds a gcd value and it mu(f) value.
 */
class MobiusPair{
public:
	ZZ gcd;
	ZZ mu;
	bool mobiusValid; //!< if true, then the mu value for this gcd has already been computed.

	MobiusPair();
	MobiusPair(const ZZ& g, const ZZ& m);
};

/**
 * Holds all the gcd values with their mu values. Base class
 */
class MobiusList {
private:
	void computeMobius(int i);
public:
	vector<MobiusPair> list;
	MobiusList();
	virtual ~MobiusList();

	void insertGCD(const ZZ& v);
	void computeMobius();
	void print() const;
};

/**
 * Holds the unweighted polynomial series for each gcd value.
 * Derived class.
 */
class MobiusSeriesList: public MobiusList
{
public:
	vector<GeneralMonomialSum<PeriodicFunction, int> *> unweightedSeries; //!< [i] is a polynomial series corresponding to gcd list[i].
	void computeMobius();
	virtual ~MobiusSeriesList();
};

/**
 * See http://en.wikipedia.org/wiki/Bernoulli_number
 * Used for doing the series expansion.
 */
class BernoulliFirstKind
{
private:
	vector<RationalNTL> B;
public:
	void setBernoulli(int k); //find first k Bernoulli numbers of the first kind.
	const RationalNTL& operator[](int i) const;
};



/**
 * Main class to compute the top Ehrhart polynomial for knapsacks.
 *
 * Example usage:
 * 		TopKnapsack tk;
 * 		tk.set(alpha);
 * 		tk.coeff_NminusK(k); //find the coeff of t^{N-k}
 * 		or ..
 * 		tk.coeff_topK(k); // find the coeff of t^N, t^{N-1}, ..., t^{N-k}
 * 		ofstream f(outFile.c_str());
 * 		tk.printAnswer(f); //best to save to file....answers are BIG
 *
 * 	Todo:
 * 		--find the gcd's faster
 * 		--can this class be used twice? like tk.coeff_topK(k); tk.coeff_topK(k+1);
 */
class TopKnapsack {
private:
	vec_ZZ alpha;								//!< list of coefficients [a_1, a_1, ...., a_{N+1}]
	int N;
	int order; 									//!< the k, as in N-k
	int largestSubsetComputed; //to delete ?
	bool topKTerms;								//!< if false, we are only computing coeff of t^{N-k}, else we are also computing the higher coefficients.
	MobiusSeriesList gcds;						//!< unscaled series expansions with gcd mu values. If topKTerms=true, gcds.unweightedSeries[i] contains terms that contribute to t^N, ..., t^{N-k} for the ith gcd only. Else, it just has terms that contribute to t^{N-k}
	BernoulliFirstKind bernoulli;				//!< store the bernoulli numbers so we do not have to always recompute.
	vector<PeriodicFunction> coeffsNminusk;		//!< Output vector. [i] = coeff of t^{N-i}.



	void coeff(int k); //!< start of computation

	void everyGCD(int k);


	void E(const int fIndex); //!< computes the series expansion of F(\alpha, gcd[f-index], T)
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

	void packageAnswer(); //!< computes the final coefficients and saves them in the coeffsNminusk vector.






	static void printMatrix(const mat_ZZ &A);//for debugging.

public:
	TopKnapsack();
	~TopKnapsack();
	static ZZ binomial(int n, int k);


	void set(const vec_ZZ& list);

	void coeff_NminusK(int k);
		void coeff_topK(int k);
		void printAnswer(ostream & out);
};



#endif /* TOPKNAPSACK_H_ */
