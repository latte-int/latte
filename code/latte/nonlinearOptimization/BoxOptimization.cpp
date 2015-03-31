/*
 * BoxOptimizationParameters.cpp
 *
 *  Created on: Apr 23, 2014
 *      Author: bedutra
 */

#include "BoxOptimization.h"
#include "ramon.h"
#include "nonlinearOptimization/WeightedExponentialSubs.h"

#include "LattException.h"

#include "top-knapsack/TopKnapsack.h"


ZZ oneVarPositiveSum(map<int, vector<RationalNTL> > & summationFormulas, const ZZ & lb, const ZZ & ub, int deg)
{

	//cout << "pos deg " << deg << " [" << lb << ", " << ub << "] " << endl;

	assert(lb >= 0);

	map<int, vector<RationalNTL> >::iterator it = summationFormulas.find(deg);

	if(it == summationFormulas.end())
	{
		//cout << "building formula for deg " << deg << endl;
		BernoulliFirstKind b;
		b.setBernoulli(deg+1);

		vector<RationalNTL> sum1n;
		sum1n.resize(deg+2);

		for(int k = 0; k <= deg ; ++k)
		{
			sum1n[deg+1-k] = b[k] * TopKnapsack::binomial(deg+1, k);
			sum1n[deg+1-k] /= (deg + 1);
			if ( k % 2 == 1)
				sum1n[deg+1-k].changeSign();
		}

		//cout << "sum1^" << deg << " =  ";
		//for(int k = 0; k < sum1n.size(); ++k)
		//	cout << sum1n[k] << "n^" << k << ",  ";
		//cout << endl;

		summationFormulas[deg] = sum1n;

		it = summationFormulas.find(deg);

		assert(it != summationFormulas.end());
	}//build the summation formula for \sum_1^n m^deg if this is the first time we needed this.







	RationalNTL eu, el;
	ZZ lpow,upow;
	lpow = lb - 1;
	upow = ub;
	for(int k = 1; k <= deg + 1; ++k)
	{
		eu += it->second[k]*upow;
		el += it->second[k]*lpow;
		lpow *= (lb-1);
		upow *= ub;
	}


	RationalNTL ans;
	ans = eu-el;

	assert(ans.getDenominator() == 1);
	return ans.getNumerator();
}

//fix me:
//--break into pos, neg, and zero parts?
//--save the computation of the powers.
ZZ oneVarSum(map<int, vector<RationalNTL> > & summationFormulas, const ZZ & lb, const ZZ & ub, int deg)
{

	//cout << "deg " << deg << " [" << lb << ", " << ub << "] " << endl;
	ZZ ans;

	if (deg == 0)
			ans = (ub - lb + 1);
	else if ( 0 <= lb)
	{
		ans = oneVarPositiveSum(summationFormulas, lb, ub, deg);
	}
	else if ( ub <= 0)
	{
		if (deg % 2 == 0)
			ans = oneVarPositiveSum(summationFormulas, -ub, -lb, deg);
		else
			ans = oneVarPositiveSum(summationFormulas, -ub, -lb, deg)*(-1);
	}
	else
	{
		ZZ one;
		one = 1;
		ans = oneVarPositiveSum(summationFormulas, one, ub, deg);
		ans += oneVarPositiveSum(summationFormulas, one, -lb, deg)*(deg % 2 == 0 ? 1 : -1);
		//0^deg = 0, for deg!=0
	}


	return ans;
}








BoxOptimization::BoxOptimization()
{
	U = 0;
	L = 0;
	N = 0;
	currentPower = 0;
	//theTrie = NULL;
	cacheWeights = NULL;
	summationFormulas = NULL;
}

BoxOptimization::~BoxOptimization()
{
	destroyMonomials(originalPolynomial);
	destroyMonomials(currentPolynomial);

	WeightedExponentialTable *t;
	while (	cacheWeights )
	{
		t = cacheWeights;
		cacheWeights = cacheWeights->next;
		delete t;
	}

	if ( summationFormulas)
		delete summationFormulas;
}

BoxOptimization & BoxOptimization::operator=(const BoxOptimization & rhs)
{
	U = rhs.U;
	L = rhs.L;
	N = rhs.N;

	BTrieIterator<RationalNTL, int>* itr =	new BTrieIterator<RationalNTL, int> ();
	term<RationalNTL, int>* term;

	originalPolynomial.termCount = 0;
	originalPolynomial.varCount = rhs.originalPolynomial.varCount;
	itr->setTrie(rhs.originalPolynomial.myMonomials,	rhs.originalPolynomial.varCount);
	itr->begin();
	for (term = itr->nextTerm(); term; term = itr->nextTerm())
		insertMonomial(term->coef, term->exps, originalPolynomial);

	currentPolynomial.termCount = 0;
	currentPolynomial.varCount = rhs.currentPolynomial.varCount;
	itr->setTrie(rhs.currentPolynomial.myMonomials,	rhs.currentPolynomial.varCount);
	itr->begin();
	for (term = itr->nextTerm(); term; term = itr->nextTerm())
		insertMonomial(term->coef, term->exps, currentPolynomial);

	delete itr;

	currentPower = rhs.currentPower;

	lowerBound = rhs.lowerBound;
	upperBound = rhs.upperBound;

	currentMap = rhs.currentMap;
	currentMapPower = rhs.currentMapPower;

	if ( rhs.cacheWeights)
	{
		THROW_LATTE_MSG(LattException::bug_NotImplementedHere,"cannot copy BoxOptimization with a cache table yet");
	}
	else
	{
		cacheWeights = NULL;
	}


	if ( rhs.summationFormulas)
		summationFormulas = new map<int, vector<RationalNTL> >(*rhs.summationFormulas);
	else
		summationFormulas = NULL;

}


void BoxOptimization::setPolynomial(const monomialSum & poly)
{
/*
	vector<bool> isSame;
	int ambDim = lowBound.length();
	int affineDim = 0;
	int i;
	int currentIndex;

	isSame.resize(ambDim);

	for(i = 0; i < ambDim; ++i)
		if ( lowBound[i] == upBound[i])
			isSame[i] = true;
		else
		{
			isSame[i] = false;
			++affineDim;
		}

	//copy lower/upperbound info over only for the full-dimension variables.
	lowerBound.SetLength(affineDim);
	upperBound.SetLength(affineDim);
	for( i = 0, currentIndex = 0; i < ambDim; ++i)
		if ( lowBound[i] != upBound[i])
		{
			lowerBound[currentIndex] = lowBound[i];
			upperBound[currentIndex] = upBound[i];
			++currentIndex;
		}


	N = 1;
	vec_ZZ maxBound;
	U = 0;
	maxBound.SetLength(ambDim);
	for( i = 0; i < ambDim; ++i)
	{
		maxBound[i] = max(abs(lowBound[i]), abs(upBound[i]));
		N *= (upBound[i] - lowBound[i] + 1);
	}
	cout << "Number of lattice points: " << N << endl;

	originalPolynomial.varCount = affineDim +1;
	originalPolynomial.termCount = 0;
*/
	originalPolynomial.varCount = poly.varCount +1;
	originalPolynomial.termCount = 0;


	BTrieIterator<RationalNTL, int>* pItr =	new BTrieIterator<RationalNTL, int> ();
	pItr->setTrie(poly.myMonomials,	poly.varCount);
	pItr->begin();


	term<RationalNTL, int>* term;
	int * exp;
	exp = new int[poly.varCount + 1];
	exp[0] = 0;

	/*
	isTrivial = true;
	for (term = pItr->nextTerm(); term; term = pItr->nextTerm())
	{
		ZZ boundTerm;
		boundTerm = 1;
		RationalNTL newCoef(term->coef);

		for (i = 0, currentIndex = 0; i < ambDim; ++i)
		{
			if ( isSame[i] == false)
			{
				exp[currentIndex + 1] = term->exps[i];
				++currentIndex;
				boundTerm *= power(maxBound[i], term->exps[i]);
				isTrivial = false; //inserted a monomial
			}
			else
			{
				newCoef *= power(lowBound[i], term->exps[i]); //lowBound[i] == upBound[i]
			}
			//cout << maxBound[i] << " ?? " << term->exps[i] << endl;


		}
		insertMonomial(newCoef, exp, originalPolynomial);
		//cout << newCoef << " " << boundTerm << " " << sign(newCoef) << endl;
		U += to_RR((newCoef) * (boundTerm * sign(newCoef)));
	}
	*/
	for (term = pItr->nextTerm(); term; term = pItr->nextTerm())
	{
		ZZ boundTerm;
		boundTerm = 1;
		for (int i = 0; i < poly.varCount; ++i)
		{
			exp[i + 1] = term->exps[i];
		}
		insertMonomial(term->coef, exp, originalPolynomial);

	}
	RationalNTL one;
	one = 1;
	for(int i = 1; i <poly.varCount + 1; ++i)
		exp[i] = 0;
	exp[0] = 1;

	insertMonomial(one, exp, originalPolynomial);
	delete [] exp;
}


void BoxOptimization::setBounds(const vec_ZZ &lowBound, const vec_ZZ &upBound)
{
	vec_ZZ maxBound;
	int ambDim = originalPolynomial.varCount - 1;

	//copy the bounds.
	lowerBound = lowBound;
	upperBound = upBound;


	N = 1;
	maxBound.SetLength(ambDim);
	for(int i = 0; i < ambDim; ++i)
	{
		maxBound[i] = max(abs(lowBound[i]), abs(upBound[i]));
		N *= (upBound[i] - lowBound[i] + 1);
	}
	//cout << "Number of lattice points: " << N << endl;



	BTrieIterator<RationalNTL, int>* pItr =	new BTrieIterator<RationalNTL, int> ();
	pItr->setTrie(originalPolynomial.myMonomials,	originalPolynomial.varCount);
	pItr->begin();


	term<RationalNTL, int>* term;
	int * exp;
	exp = new int[ambDim + 1];
	exp[0] = 0;


/*
	U = 0;
	ZZ boundTerm;
	for (term = pItr->nextTerm(); term; term = pItr->nextTerm())
	{
		if ( term->exps[0] != 0)
			continue;
		boundTerm = 1;

		RationalNTL newCoef(term->coef);

		for (int i = 0; i < ambDim; ++i)
		{
			boundTerm *= power(maxBound[i], term->exps[i+1]);
		}
		U += to_RR((newCoef) * (boundTerm * sign(newCoef)));
	}
	L = U;
	L *= -1;
	cout << "Simple bounds " << L << "<= f(x) <=" << U << endl;
*/
	//starta
	pItr->begin();
	ZZ u,l;
	u=0;
	l=0;
	for(term = pItr->nextTerm(); term; term = pItr->nextTerm())
	{
		if (term->exps[0] != 0)
			continue;

		ZZ uu,ll;
		uu = term->coef.getNumerator();
		assert(term->coef.getDenominator() == 1);
		ll = uu;
		for(int i = 0; i < ambDim; ++i)
		{
			ZZ a,b;
			if ( term->exps[i+1] == 0 )
				continue;
			if( term->exps[i+1] % 2 == 0 && upBound[i] > 0 && lowBound[i] < 0)
			{
				a = power(max(abs(upBound[i]), abs(lowBound[i])), term->exps[i+1]); //power(max(abs(upBound[i]), abs(lowBound[i]), term->exps[i]));
				b = 0;
			}
			else
			{
				a = power(upBound[i], term->exps[i+1]);
				b = power(lowBound[i], term->exps[i+1]);

			}
			//cout << "uu=" << uu << ", ll=" << ll << endl;
			//cout << "a=" << a << ", b=" << b << endl;
			ZZ tempUU;
			tempUU = max(uu*a, max(uu*b, max(ll*a,ll*b)));
			ZZ tempLL;
			tempLL = min(uu*a, min(uu*b, min(ll*a,ll*b)));
			uu = tempUU;
			ll = tempLL;
			//cout << "nuu=" << uu << ", nll=" << ll << endl;

		}
		u += uu;
		l += ll;

	}
	U = to_RR(u);
	L = to_RR(l);
	//cout << "Smarter bounds " << l << "<= f(x) <=" << u << endl;



	delete [] exp;

	//cout << "updated poly:" << originalPolynomial.varCount << ", " << originalPolynomial.termCount << endl;
	//cout << "updated poly: " << printMonomials(originalPolynomial).c_str() << endl;
	//cout << "isTrivial" << isTrivial << endl;
}

void BoxOptimization::printNumberOfPoints() const
{
	cout << "Number of lattice points: " << N << endl;
}

bool BoxOptimization::isTrivial()
{
	return (N <= 1);
}

void BoxOptimization::enumerateProblem(const vec_ZZ &lowBound, const vec_ZZ &upBound, const monomialSum & poly)
{
	RationalNTL answer;
	if ( N != 1)
		THROW_LATTE_MSG(LattException::bug_NotImplementedHere, "enumerateProblem can only work with N=1 for right now");
	// so lowBound == upBound

	BTrieIterator<RationalNTL, int>* pItr =	new BTrieIterator<RationalNTL, int> ();
	pItr->setTrie(poly.myMonomials,	poly.varCount);
	pItr->begin();

	//evaluate f at one point.
	term<RationalNTL, int>* term;
	for (term = pItr->nextTerm(); term; term = pItr->nextTerm())
	{
		RationalNTL newCoef(term->coef);

		for (int i = 0; i < poly.varCount; ++i)
		{
			newCoef *= power(lowBound[i], term->exps[i]);
		}

		answer += newCoef;
	}

	U = to_RR(answer);
	L = to_RR(answer);
}


void BoxOptimization::findRange(int itr)
{


	//cout << "k: " << currentPower << " Starting range " << L << " f(x) " << U << endl;
	//cout << "k: " << currentPower << " Starting U: " << U << endl;

	RR oldU, oldL;
	for(int i = 0; i < itr; ++i)
	{
		oldU = U;
		oldL = L;
		findNewUpperbound();
		findNewLowerbound();

		//cout << "k: " << currentPower << " range " << L << " <=f(x)<= " << U << " max: " << maximumLowerBound() << " <=max(f)<= " << maximumUpperbound() << " gap:" << maximumUpperbound() - maximumLowerBound() << endl;


		if ( (oldU - U)  + (L - oldL) < 1)
			return; //difference between new and old range bounds is less than one, so just stop.
	}
}

// 0 <= f - L
// ==> f - L <= u
// f <= u + L := U_{i+1}
void BoxOptimization::findNewUpperbound()
{
	//currentPolynomial = (f + s)^k
	// f-L <= [\sum(f-L)^k]^{1/k} = [\sum(f+(-L))^k]^{1/k}
	//so s = -L;

	RR ans;
	ans = 0;

	RR s(L);
	s *= -1;

	ans = currentMap.eval(s);

	//cout << "(f+s)^" << currentPower << " ans=" << ans  << " with s=" << s << endl;
	RR newU;
	newU = L + pow(ans, to_RR(1)/to_RR(currentPower));
	if (newU < U)
		U = newU;
	else
	cout << newU << " newU >= U " << U << endl;
}

// 0 <= U - f
// ==> U - f  <= u
// L_{i+1} := U - u <= f
void BoxOptimization::findNewLowerbound()
{
	//currentPolynomial = (f + s)^k
	// U-f <= [\sum(U-f)^k]^{1/k} = [(-1)^k * \sum(f-U)^k]^{1/k}
	//so s = -U;

	RR ans;
	ans = 0;

	RR s(U);
	s *= -1;

	ans = currentMap.eval(s);

	if ( currentPower % 2)
		ans *= -1;

	//cout << "(s-f)^" << currentPower << " ans=" << ans << endl;
	RR newL;
	newL = U - pow(ans, inv(to_RR(currentPower)));
	if ( newL > L)
		L = newL;
	else
		cout << newL << " newL <= L " << L << endl;
}

RR BoxOptimization::maximumUpperbound()
{
	// 0 <= f - L <= u, where u:=pow(currentMap.eval(-L), to_RR(1)/to_RR(currentPower));
	return L + pow(currentMap.eval(-L), to_RR(1)/to_RR(currentPower));
}

RR BoxOptimization::maximumLowerBound()
{
	// 0 <= f - L <= u
	//then u/N^{1/k} <= max (f - L)
	return L + pow(currentMap.eval(-L), to_RR(1)/to_RR(currentPower))/pow(to_RR(N), inv(to_RR(currentPower)));
}


RR BoxOptimization::sampleLowerBound(monomialSum &poly, const vec_ZZ & point)
{
	RationalNTL ans;
	BTrieIterator<RationalNTL, int>* pItr =	new BTrieIterator<RationalNTL, int> ();
	pItr->setTrie(poly.myMonomials,	poly.varCount);
	pItr->begin();

	term<RationalNTL, int>* term;

	for (term = pItr->nextTerm(); term; term = pItr->nextTerm())
	{
		//cout << "term " << term->coef << endl;
		RationalNTL value;
		value = 1;
		for (int currentPower = 0; currentPower < poly.varCount; ++currentPower)
		{
			//cout << term->exps[currentPower] << ".";
			value *= power(point[currentPower], term->exps[currentPower]);
		}
		//cout << endl;
		value *= term->coef;
		ans += value;
	}

	delete pItr;
	//cout << ans << endl;
	return to_RR(ans);
}

/**
 * @param k: sets currentPolynomial = originalPolynomial^k if k > currentPower, else the currentPolynomial is not changed.
 * todo: maybe delete current poly as it is not needed. Also, are we deleting the linear forms 2x?
 */
void BoxOptimization::setPower(int k)
{
	assert(k >= 1);
	//currentPower = k;

	if ( k <= currentPower)
		return; //we only want to increase the power of the currentPolynomial.

	cout << "computing (f(x) + s)^" << k << "...\n" << flush;

	if( currentPower == 0)
	{
		currentPolynomial.varCount = originalPolynomial.varCount;
		currentPolynomial.termCount = 0;

		BTrieIterator<RationalNTL, int>* pItr =	new BTrieIterator<RationalNTL, int> ();
		pItr->setTrie(originalPolynomial.myMonomials,	originalPolynomial.varCount);
		pItr->begin();


		term<RationalNTL, int>* term;
		for (term = pItr->nextTerm(); term; term = pItr->nextTerm())
			insertMonomial(term->coef, term->exps, currentPolynomial);
		delete pItr;

		currentPower = 1;
	}//copy the origional polynomisl (f(x)+s) into currentPolynomial.



	//********************************************
	//next, take currentPolynomial to the kth power.

	BTrieIterator<RationalNTL, int>* it1 = new BTrieIterator<RationalNTL, int>();
	BTrieIterator<RationalNTL, int>* it2 = new BTrieIterator<RationalNTL, int>();
	it1->setTrie(originalPolynomial.myMonomials, originalPolynomial.varCount);
	it2->setTrie(originalPolynomial.myMonomials, originalPolynomial.varCount);

	for( int i  = currentPower; i < k; ++i)
	{
		it1->setTrie(originalPolynomial.myMonomials, originalPolynomial.varCount);
		it2->setTrie(currentPolynomial.myMonomials, currentPolynomial.varCount);

		monomialSum temp;
		temp.varCount = currentPolynomial.varCount;
		multiply(it1, it2, temp);
		destroyMonomials(currentPolynomial);
		currentPolynomial = temp;
		cout << "  finished computing (f(x) + s)^ " << i+1 << "\n";
	}
	currentPower = k;
	//cout << "power poly: " << printMonomials(currentPolynomial).c_str() << endl;
	//cout << "power poly deg: " << k << endl;


	//****************************************************
	//just for fun, lets count the number of terms we have
	it2->setTrie(currentPolynomial.myMonomials, currentPolynomial.varCount);
	it2->begin();
	term<RationalNTL, int> * t;
	long numMonomials = 0;
	for(t = it2->nextTerm(); t; t = it2->nextTerm())
		++numMonomials;


	delete it1;
	delete it2;
	cout << "done. " << numMonomials << " monomials computed\n";





}


/**
 * algoType:
 * 		lf      -- decompose to linear forms and find the s-polynomial
 * 		lfCache -- decompose to linear forms and build the cache table
 * 		naturalSummation -- do nothing here.
 */
void BoxOptimization::decomposePoly(AlgoType algoType)
{

	if (algoType != BoxOptimization::lf && algoType != BoxOptimization::lfCache)
		return;

	//****************************************************
	//next, decompose currentPolynomial into linear forms.

	cout << "decomposing (f(x) + s)^" << currentPower << " into powers of linear forms..." << flush;
	BTrieIterator<RationalNTL, int>* pitr =	new BTrieIterator<RationalNTL, int> ();
	pitr->setTrie(currentPolynomial.myMonomials, currentPolynomial.varCount);
	pitr->begin();

	//if(theTrie)
	//	delete theTrie;
	BurstTrie<PolynomialMap, ZZ> *theTrie;
	theTrie = new BurstTrie<PolynomialMap, ZZ>;


	BTrieIterator<RationalNTL, ZZ>* tempLFItr = new BTrieIterator<RationalNTL, ZZ> ();

	term<RationalNTL, int>* currentMonomial;
	ZZ numLinForms;
	numLinForms = 0;
	for (currentMonomial = pitr->nextTerm(); currentMonomial; currentMonomial = pitr->nextTerm())
	{
		//cout << "going to decompose " << currentMonomial->coef << "* ";
		//for(int i = 0; i < currentMonomial->length; ++i)
		//	cout << currentMonomial->exps[i] << "  ";
		//cout << "] deg=" << currentMonomial->degree << ", len=" << currentMonomial->length << endl;

		term<RationalNTL, int> originalMonomial;
		originalMonomial.coef = currentMonomial->coef;
		originalMonomial.degree = currentMonomial->degree;
		originalMonomial.exps = currentMonomial->exps + 1;
		originalMonomial.length = currentMonomial->length -1;

		//decompose the monomial into a power of linear forms.
		//s^k *monomial -> s^k *( c_1\ell_1^m_1 + ....)
		linFormSum tempLF;
		tempLF.varCount = originalMonomial.length;
		decompose(&originalMonomial, tempLF);

		//loop over the linear forms.
		//BTrieIterator<RationalNTL, ZZ>* tempLFItr = new BTrieIterator<RationalNTL, ZZ> ();
		tempLFItr->setTrie(tempLF.myForms, tempLF.varCount);
		tempLFItr->begin();

		term<RationalNTL, ZZ>* tempLFTerm;
		for(tempLFTerm =  tempLFItr->nextTerm(); tempLFTerm; tempLFTerm =  tempLFItr->nextTerm())
		{
			PolynomialMap newCoef;
			newCoef.terms[currentMonomial->exps[0]] = tempLFTerm->coef;

			//insert each linear form back in.
			theTrie->insertTerm(newCoef, tempLFTerm->exps, 0, tempLFTerm->length, tempLFTerm->degree);
			assert(tempLFTerm->length == originalMonomial.length);
		}

		destroyLinForms(tempLF);
	}
	delete tempLFItr;
	cout << "done.\n";

	//****************************************************
	//next, integrate the linear forms.

	cout << "Performing the weighted counting..." << flush;
	currentMap.terms.clear();
	BTrieIterator<PolynomialMap, ZZ>* mapitr = new BTrieIterator<PolynomialMap, ZZ> ();
	mapitr->setTrie(theTrie, lowerBound.length());
	mapitr->begin();


	while (mapitr->nextTerm())
		++numLinForms;
	mapitr->begin();

	term<PolynomialMap, ZZ>* mapitrTerm;
	RationalNTL one(1,1);
	int progress = 0;
	WeightedCountingBuffer wcb;
	for(mapitrTerm =  mapitr->nextTerm(); mapitrTerm; mapitrTerm =  mapitr->nextTerm())
	{
		//cout << mapitrTerm->coef << "*( ";
		//for(int j = 0; j < lowerBound.length(); ++j)
		//	cout << mapitrTerm->exps[j] << ", ";
		//cout << " )^" << mapitrTerm->degree;



		if (algoType == BoxOptimization::lfCache)
		{
			WeightedExponentialTable* t = computeWeightedCountingBox_singleForm(wcb, lowerBound.length(), mapitrTerm->exps, mapitrTerm->degree);

			//cout << " = " << temp << endl;
			t->linFormPow = mapitrTerm->degree;
			t->linForm.SetLength(lowerBound.length());
			for(int i = 0; i <lowerBound.length(); ++i)
				t->linForm[i] = mapitrTerm->exps[i];
			t->sPoly = mapitrTerm->coef;
			t->next = cacheWeights;
			cacheWeights = t;
		} //cache table
		if ( algoType == BoxOptimization::lf)
		{
			mpq_class temp = computeWeightedCountingBox_singleForm(wcb, lowerBound, upperBound, mapitrTerm->exps, mapitrTerm->degree, one);
			RationalNTL rTemp;
			rTemp = convert_mpq_to_RationalNTL(temp);
			//cout << rTemp << "*(";
			//mapitrTerm->coef.print(cout);
			mapitrTerm->coef.mult(rTemp);
			//mapitrTerm->coef.print(cout);


			currentMap += mapitrTerm->coef;
			//currentMap.print(cout);
		}

		++progress;
		if ( progress % 100 == 0)
		{
			cout << "linear forms processed " << progress << "/" <<  numLinForms << endl;
		}

	}
	cout << "done." << endl;
	delete theTrie;
	delete mapitr;
}


void BoxOptimization::printSpolynomial() const
{
	cout << "Final spoly "<< currentMap <<endl;
}

/**
 * algoType:
 * 		lf      -- do nothing, s-polynomial already computed
 * 		lfCache -- use the cache table to build the s-polynomial
 * 		naturalSummation -- build the s-polynomial
 */
void BoxOptimization::findSPolynomial(AlgoType algoType, const vec_ZZ &lowerBound, const vec_ZZ & upperBound)
{
	if (algoType == BoxOptimization::lfCache)
		findSPolynomial_lfCache(lowerBound, upperBound);
	if (algoType == BoxOptimization::naturalSummation)
		findSPolynomial_naturalSummation(lowerBound, upperBound);



}

void BoxOptimization::findSPolynomial_naturalSummation(const vec_ZZ &lowerBound, const vec_ZZ &upperBound)
{
	PolynomialMap spoly;
	int dim = lowerBound.length();

	if (summationFormulas == NULL)
		summationFormulas = new map<int, vector<RationalNTL> >;


	BTrieIterator<RationalNTL, int>* pitr =	new BTrieIterator<RationalNTL, int> ();
	pitr->setTrie(currentPolynomial.myMonomials, currentPolynomial.varCount);
	pitr->begin();
	term<RationalNTL, int>* term;
	for(term =  pitr->nextTerm(); term; term =  pitr->nextTerm())
	{
		RationalNTL prod(1,1);

		for(int i = 1; i<= dim; ++i)
			prod *= oneVarSum(*summationFormulas, lowerBound[i-1], upperBound[i-1], term->exps[i]);



		spoly.terms[term->exps[0]] += prod * term->coef;

	}//for each monomial,


	delete pitr;
	currentMap = spoly;
}

void BoxOptimization::findSPolynomial_lfCache(const vec_ZZ &lowerBound, const vec_ZZ & upperBound)
{
	PolynomialMap spoly;
	int dim = originalPolynomial.varCount - 1;
	int n; //dim of amb box
	int index;
	ZZ two;
	ZZ twon;
	ZZ vDotL;
	ZZ Namb; //number of int points in the projected smaller box
	ZZ Nother; //number of int points that was projected out.

	//cout << "start s poly "<< endl;
	for(WeightedExponentialTable * ptr = cacheWeights; ptr; ptr = ptr->next)
	{

		//cout << ptr->linForm << " pow " <<ptr->linFormPow << " num vertex " << ptr->weights.size() << " M+d " << ptr->weights[0].size() <<endl;
		//cout << "got here" <<endl;
		n = 0;
		Namb = 1;
		Nother = 1;

		for(int i = 0; i <dim; ++i)
		{
			//cout << "i " <<i << " dim" <<dim <<endl;
			if (ptr->linForm[i] != 0)
			{
				++n; //dim of amb box
				Namb *= (upperBound[i] - lowerBound[i] + 1);
			}
			else
			{
				Nother *= (upperBound[i] - lowerBound[i] + 1);
			}
		}

		two = to_ZZ(2);
		twon = power(two, n);//2^n

		//cout << "ptr" << ptr <<" lf " << ptr->linForm << "^" << ptr->linFormPow << " non-zero's "<< n << endl;

		//cout <<"found num N" <<endl;

		int wIndex = 0;
		mpq_class_lazy oneCone;
		for(ZZ i = to_ZZ(0); i < twon; ++i, ++wIndex)
		{


			vDotL = 0;
			index = 0;
			//if (ptr->linForm[0] == 2 && ptr->linForm[1] == 4 && ptr->linForm[2] == 1)
			//	cout << " vertex ";
			//use the binary representation of i to pick the current vertex along with its tangent cone.
			for(int j = 0; j < n; ++j)
			{
				while(ptr->linForm[index] == 0)
					++index;

				if ( bit(i, j) )
				{
					vDotL += ptr->linForm[index] * upperBound[index];
					//if (ptr->linForm[0] == 2 && ptr->linForm[1] == 4 && ptr->linForm[2] == 1)
					//	cout << ", " << upperBound[index];
				}
				else
				{
					vDotL += ptr->linForm[index] * lowerBound[index];
					//if (ptr->linForm[0] == 2 && ptr->linForm[1] == 4 && ptr->linForm[2] == 1)
					//	cout << ", " << lowerBound[index];
				}
				++index;
			}
			mpq_class_lazy vDotLmpq = convert_ZZ_to_mpq(vDotL);

			//check if vDotL is zero.

			//cout << "vdotl "<< vDotL <<endl;
			mpq_class_lazy oneVertex;
			mpq_class_lazy inner = 1;
			for(int j = 0; j <= n+ptr->linFormPow; ++j)
			{
				oneVertex += inner * ptr->weights[wIndex][j];
				inner = inner * vDotLmpq;
			}

			//if (ptr->linForm[0] == 2 && ptr->linForm[1] == 4 && ptr->linForm[2] == 1)
			//	cout << "oneVertex = " << oneVertex <<endl;
			oneCone += oneVertex;

		}//for each vertex of a box.

		oneCone *= convert_ZZ_to_mpq(Nother);

		//cout << ptr->sPoly << " * " << oneCone << endl;

		spoly.add(ptr->sPoly, convert_mpq_to_RationalNTL(oneCone));

	}//for each linear form

	cout << "Final spoly "<< spoly <<endl;
	currentMap = spoly;
}



WeightedBoxProducer::WeightedBoxProducer(const vec_ZZ & lowerb, const vec_ZZ & upperb):
		lowerBound(lowerb), upperBound(upperb)
{
}

WeightedBoxProducer::~WeightedBoxProducer() {}


void WeightedBoxProducer::Produce(ConeConsumer &consumer)
{
	unsigned int n = upperBound.length();
	ZZ two = to_ZZ(2);
	ZZ twon = power(two, n);//2^n

	vec_ZZ allOnes, theVertex;
	allOnes.SetLength(n);
	theVertex.SetLength(n);
	for(int i = 0; i < n; ++i)
		allOnes[i] = 1;

	mat_ZZ id; // n by n idenity matrix.
	ident(id, n);


	for(ZZ i = to_ZZ(0); i < twon; ++i)
	{

		listCone *oneCone = createListCone();
		oneCone->determinant = 1;
		oneCone->dual_determinant = 1;

		//use the binary representation of i to pick the current vertex along with its tangent cone.
		int det = 0;
		for(int j = 0; j < n; ++j)
		{
			if ( bit(i, j) )
			{
				theVertex[j] = upperBound[j];
				oneCone->rays = appendVectorToListVector(id[j]*(-1), oneCone->rays);
				//oneCone->facets = appendVectorToListVector(id[j], oneCone->facets);
				det++;
			}
			else
			{
				theVertex[j] = lowerBound[j];
				oneCone->rays = appendVectorToListVector(id[j], oneCone->rays);
				//oneCone->facets = appendVectorToListVector(id[j]*(-1), oneCone->facets);
			}
		}

		//determinant = -1^(n/2) * -1^(det)
		//                |           |
		//                |         (number of negative rays)
		//             (sorting columns to identity)
		if ( (det + n/2) % 2)
			oneCone->determinant = -1;
		if ( (n/2 + n-det) % 2)
			oneCone->dual_determinant = -1;

		oneCone->vertex = new Vertex(new rationalVector(theVertex, allOnes));

		//oneCone->facet_divisors = allOnes; //need this if computing facets too.

		oneCone->latticePoints = createListVector(theVertex); //cone is unimodular


		//finally, process this tangent cone.
		consumer.ConsumeCone(oneCone);
	}//for each vertex of a box.
}//WeightedBoxProducer::Produce








/**
 * assumes lowerbound[i] < upperBound[i]
 */
mpq_class computeWeightedCountingBox(const vec_ZZ &lowerBound, const vec_ZZ &upperBound, const linFormSum &originalLinearForms)
{
	mpq_class ans;

	ans = 0;

	BTrieIterator<RationalNTL, ZZ>* linearFormsItr = new BTrieIterator<RationalNTL, ZZ> ();
	linearFormsItr->setTrie(originalLinearForms.myForms, originalLinearForms.varCount);
	linearFormsItr->begin();


	term<RationalNTL, ZZ>* lform;

	//loop over the linear forms
	WeightedCountingBuffer wcb;
	for (lform = linearFormsItr->nextTerm(); lform; lform
			= linearFormsItr->nextTerm())
	{
		ans += computeWeightedCountingBox_singleForm(wcb, lowerBound, upperBound, lform->exps, lform->degree, lform->coef);
		//cout << "Running sum: " << ans << endl;
	}//for every term in the originalPolynomial

	delete linearFormsItr;


	return ans;
}

mpq_class computeWeightedCountingBox_singleForm(WeightedCountingBuffer & wcb, const vec_ZZ &lowerBound, const vec_ZZ &upperBound, const ZZ* linFormExps, const int degree, const RationalNTL & coef)
{
	mpq_class ans;
	ans = 0;

	//project the box into a smaller dimension if one of the linear form coeffs is zero.
	vec_ZZ linFormCoeffs, newLB, newUB;
	int n = lowerBound.length();
	linFormCoeffs.SetLength(n);
	newLB.SetLength(n);
	newUB.SetLength(n);

	ZZ zeroTerms;
	zeroTerms = 1;

	int j = 0;
	for(int i = 0; i < n; ++i)
	{
		if ( linFormExps[i] != 0)
		{
			linFormCoeffs[j] = linFormExps[i];
			newLB[j] = lowerBound[i];
			newUB[j] = upperBound[i];
			++j;
		}
		else
			zeroTerms *= (upperBound[i] - lowerBound[i] + 1);
	}
	linFormCoeffs.SetLength(j); //ignore the zeros at the end if any.
	newUB.SetLength(j);
	newLB.SetLength(j);

	//easy special cases.
	if ( coef == 0 ||  (j == 0 && degree != 0))
		return ans; // if coeff is zero or the linear form is the zero polynomial.
	if ( degree == 0)
	{
		ZZ t;
		t = 1;
		for(int i = 0; i < lowerBound.length(); ++i)
			t *= (upperBound[i] - lowerBound[i] + 1);
		ans = mpq_class(convert_ZZ_to_mpz(coef.getNumerator()), convert_ZZ_to_mpz(coef.getDenominator())) * convert_ZZ_to_mpz(t);
		return ans;
	}//if the linear form is a constant




	//construct a producer-consumer pair.
	WeightedBoxProducer wbp(newLB, newUB);
	Weighted_Exponential_Single_Cone_Parameters wParams;

	wParams.substitution = BarvinokParameters::PolynomialSubstitution;
	wParams.decomposition = BarvinokParameters::DualDecomposition;
	wParams.max_determinant = 1;
	wParams.Number_of_Variables = j;
	wParams.InitializeComputation();
	wParams.linForm = linFormCoeffs;
	wParams.linFormPow = degree;
	wParams.wcb = &wcb;
	wbp.Produce(wParams);

	ans = mpq_class(convert_ZZ_to_mpz(coef.getNumerator()), convert_ZZ_to_mpz(coef.getDenominator())) * wParams.result * convert_ZZ_to_mpz(zeroTerms);
/*
	cout << " sum of " << lowerBound << endl;
	cout << "        " << upperBound << endl;
	cout << "        " << coef << "*(";
	for(int i = 0;i < lowerBound.length(); ++i)
		cout << linFormExps[i] << "*x[" << i << "] + ";
	cout << ")^" << degree << endl;
	cout << "        " << " is " << ans << endl;
*/

	return ans;

}


//only the weights info will be saved.
WeightedExponentialTable* computeWeightedCountingBox_singleForm(WeightedCountingBuffer & wcb, const int n, const ZZ* linFormExps, const int degree)
{
	mpq_class ans;
	ans = 0;

	//project the box into a smaller dimension if one of the linear form coeffs is zero.
	vec_ZZ linFormCoeffs, newLB, newUB;
	linFormCoeffs.SetLength(n);
	newLB.SetLength(n);
	newUB.SetLength(n);

	int j = 0;
	for(int i = 0; i < n; ++i)
	{
		if ( linFormExps[i] != 0)
		{
			linFormCoeffs[j] = linFormExps[i];
			++j;
		}
	}
	linFormCoeffs.SetLength(j); //ignore the zeros at the end if any.
	newUB.SetLength(j);
	newLB.SetLength(j);


	//construct a producer-consumer pair.
	WeightedBoxProducer wbp(newLB, newUB);
	Weighted_Exponential_Single_Cone_Parameters_BranchBound wParams;

	wParams.substitution = BarvinokParameters::PolynomialSubstitution;
	wParams.decomposition = BarvinokParameters::DualDecomposition;
	wParams.max_determinant = 1;
	wParams.Number_of_Variables = j;
	wParams.InitializeComputation();
	wParams.linForm = linFormCoeffs;
	wParams.linFormPow = degree;
	wParams.wcb = &wcb;
	wbp.Produce(wParams);


	return wParams.table;
}

