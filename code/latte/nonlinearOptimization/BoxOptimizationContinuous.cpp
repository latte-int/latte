/*
 * BoxOptimizationContinuous.cpp
 *
 *  Created on: Mar 31, 2015
 *      Author: bedutra
 */

#include "BoxOptimizationContinuous.h"

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




BoxOptimizationContinuous::BoxOptimizationContinuous()
{
	U = 0;
	L = 0;
	V = 0;
	currentPower = 0;
}

BoxOptimizationContinuous::~BoxOptimizationContinuous()
{
	destroyMonomials(originalPolynomial);
	destroyMonomials(currentPolynomial);
}

BoxOptimizationContinuous & BoxOptimizationContinuous::operator=(const BoxOptimizationContinuous & rhs)
{
	U = rhs.U;
	L = rhs.L;
	V = rhs.V;

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


	return *this;
}


void BoxOptimizationContinuous::setPolynomial(const monomialSum & poly)
{

	originalPolynomial.varCount = poly.varCount +1;
	originalPolynomial.termCount = 0;


	BTrieIterator<RationalNTL, int>* pItr =	new BTrieIterator<RationalNTL, int> ();
	pItr->setTrie(poly.myMonomials,	poly.varCount);
	pItr->begin();


	term<RationalNTL, int>* term;
	int * exp;
	exp = new int[poly.varCount + 1];
	exp[0] = 0;

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


void BoxOptimizationContinuous::setBounds(const vec_RRorQQ &lowBound, const vec_RRorQQ &upBound)
{
	int ambDim = originalPolynomial.varCount - 1;
	bool notAPoint = false;
	//copy the bounds.
	lowerBound = lowBound;
	upperBound = upBound;


	V = 1;
	for(int i = 0; i < ambDim; ++i)
	{
		if ( lowBound[i] + zero < upBound[i])
		{
			V *= (upBound[i] - lowBound[i]);
			notAPoint = true;
		}
	}
	//cout << "Number of lattice points: " << N << endl;

	if ( notAPoint == false)
		V = 0; //really small box...this is just a point.


	BTrieIterator<RationalNTL, int>* pItr =	new BTrieIterator<RationalNTL, int> ();
	pItr->setTrie(originalPolynomial.myMonomials,	originalPolynomial.varCount);
	pItr->begin();


	term<RationalNTL, int>* term;
	//int * exp;
	//exp = new int[ambDim + 1];
	//exp[0] = 0;

	pItr->begin();
	RRorQQ u,l;
	u=0;
	l=0;
	for(term = pItr->nextTerm(); term; term = pItr->nextTerm())
	{
		//cout << "   monomial " << term->coef;
		//for (int i = 0; i < ambDim; ++i)
		//	cout << " * x[" << i << "]^" << term->exps[i+1];
		//cout << endl;
		if (term->exps[0] != 0)
			continue;

		RRorQQ uu,ll;
		uu = term->coef.to_RR();
		//uu = term->coef;
		//		to_RR(term->coef.getNumerator());
		//assert(term->coef.getDenominator() == 1);
		ll = uu;
		for(int i = 0; i < ambDim; ++i)
		{
			RRorQQ a,b;
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

			//cout << i << " up, lower bound" << lowBound[i] << ", " << upBound[i] << ", pow " <<term->exps[i+1] << endl;
			//cout << "uu=" << uu << ", ll=" << ll << endl;
			//cout << "a=" << a << ", b=" << b << endl;
			RRorQQ tempUU;
			tempUU = max(uu*a, max(uu*b, max(ll*a,ll*b)));
			RRorQQ tempLL;
			tempLL = min(uu*a, min(uu*b, min(ll*a,ll*b)));
			uu = tempUU;
			ll = tempLL;
			//cout << "nuu=" << uu << ", nll=" << ll << endl;

		}
		//cout << "adding uu, ll = " << uu << ", "<< ll << endl;
		u += uu;
		l += ll;

	}
	U = u;
	L = l;
	//cout << "Smarter bounds " << l << "<= f(x) <=" << u << endl;

	//delete [] exp;


	//Now find the lipschitz constant times M

	pItr->begin();
	lipschitz = 0;
	M = 0;
	RR maxB = upperBound[0];
	for(int i = 0; i < ambDim; ++i)
	{
		if ( M < upperBound[i] - lowerBound[i])
			M =upperBound[i] - lowerBound[i];
		maxB = max(maxB, max(abs(upperBound[i]), abs(lowerBound[i])));
	}


	lipschitz = 0;
	for(term = pItr->nextTerm(); term; term = pItr->nextTerm())
	{
		int deg = 0;
		if ( term->exps[0] != 0)
			continue;

		for (int i = 0; i < ambDim; ++i)
			deg += term->exps[i+1];
		lipschitz += abs(term->coef.to_RR()) * deg* power(maxB, (deg-1));
		//lipschitz += abs(term->coef) * power(M, (deg-1)) * deg;


	}
	//cout << "old lipschitz=" <<lipschitz << endl;

	vec_RRorQQ maxLen, del;
	maxLen.SetLength(ambDim);
	del.SetLength(ambDim);
	pItr->begin();

	for(int i = 0; i < ambDim; ++i)
		maxLen[i] = max(abs(upperBound[i]), abs(lowerBound[i]));

	for(term = pItr->nextTerm(); term; term = pItr->nextTerm())
	{
		for(int i = 0; i <ambDim; ++i)
		{
			RR l1bound;
			l1bound = 1;
			if ( term->exps[i+1] == 0)
				continue; //no x[i] var in this monomial.

			for(int j = 0; j < ambDim; ++j)
			{
				if (j == i)
					l1bound *= power(maxLen[j], term->exps[j+1] - 1 ) * term->exps[j+1];
				else
					l1bound *= power(maxLen[j], term->exps[j+1]);
			}

			l1bound *= abs(term->coef.to_RR());
			del[i] += l1bound;
		}//differentiate w.r.t x[i]
	}

	RR newLipschitz;
	for(int i = 0; i < ambDim; ++i)
		newLipschitz += del[i];

	//cout << "newLipschitz=" <<newLipschitz << endl;
	lipschitz = min(lipschitz, newLipschitz);

	//for(int i = 0; i < ambDim; ++i)
	//	cout << "x" << i << ": " << lowBound[i] << ", " << upBound[i] << endl;
	//cout << "\nM=" << M << ", Lip=" << lipschitz << ", V=" << V << endl;
	//cout << "L=" << L << ", U=" << U << endl;

	delete pItr;
}

void BoxOptimizationContinuous::printStats() const
{
	cout << "Volume of box: " << V
		 << "\nLipschitz    : " << lipschitz
		 << "\nL            : " << L
		 << "\nU            : " << U
		 << "\nM            : " << M << "\n";

}

void BoxOptimizationContinuous::setFmaxLowerbound(const double & fmaxLB)
{
	fmaxLowerBound = fmaxLB;
}


RR BoxOptimizationContinuous::maximumLowerBound()
{
	// 0 <= f - L <= u
	//then u/N^{1/k} <= max (f - L)
	//cout << "eval" << endl;
	//cout <<evalSpoly(-L) << endl;
	//cout << "top" << endl;
	//cout << pow(evalSpoly(-L), inv(to_RR(currentPower))) << endl;
	//cout << "bot" << endl;
	//cout << pow(V, inv(to_RR(currentPower))) << endl;


	return L + pow(evalSpoly(-L), inv(to_RR(currentPower)))/pow(V, inv(to_RR(currentPower)));
}


RR BoxOptimizationContinuous::maximumUpperbound()
{
	int d = originalPolynomial.varCount -1;

	RR p1, p2;
	p1 = d;
	p1 /= d+currentPower;
	p2 = currentPower;
	p2 /= d+currentPower;
	//p1 = d/(d+k)
	//p2 = k/(d+k)

	RR f1, f2, f3;

	if ( p1 * (U -L) / (M*lipschitz) < 1) // p1 *fmax/ML < 1
	{


		f1 = pow(evalSpoly(-L)/V, inv(to_RR(d+currentPower)));
		f2 = pow( (M*lipschitz * (d+currentPower)) * inv(to_RR(d)), p1);
		f3 = pow( to_RR(d+currentPower)/to_RR(currentPower), p2);

		//cout << "s-poly at " << -L << " is " << evalSpoly(-L) << endl;
		//cout << f1 << " * " << f2 << " *  " << f3 << " = " << f1*f2*f3 <<endl;
		cout << "maximumUpperbound()::fmax used\n";
		return L + f1*f2*f3;
	}
	else if (p1*(fmaxLowerBound - L)/ (M*lipschitz) < 1) //p1 * f(something)/ML < 1
	{
		f1 = pow(evalSpoly(-L)/V, inv(to_RR(currentPower)));
		f2 = pow(p1*(fmaxLowerBound - L)/ (M*lipschitz), to_RR(d)/currentPower) * p2;

		cout << "maximumUpperbound()::f(something) used\n";
		return L + f1*inv(f2);
	}
	else
	{
		f1 = pow(evalSpoly(-L)/V, inv(to_RR(currentPower)));
		cout << "maximumUpperbound():: > 1\n";
		return L +f1*inv(p2);
	}
}


/*
RR BoxOptimizationContinuous::sampleLowerBound(monomialSum &poly, const vec_RR & point)
{
	RR ans;
	BTrieIterator<RationalNTL, int>* pItr =	new BTrieIterator<RationalNTL, int> ();
	pItr->setTrie(poly.myMonomials,	poly.varCount);
	pItr->begin();

	term<RationalNTL, int>* term;

	for (term = pItr->nextTerm(); term; term = pItr->nextTerm())
	{
		//cout << "term " << term->coef << endl;
		RR value;
		value = 1;
		for (int currentPower = 0; currentPower < poly.varCount; ++currentPower)
		{
			//cout << term->exps[currentPower] << ".";
			value *= power(point[currentPower], term->exps[currentPower]);
		}
		//cout << endl;
		value *= to_RR(term->coef);
		ans += value;
	}

	delete pItr;
	//cout << ans << endl;
	return ans;
}
*/

/**
 * @param k: sets currentPolynomial = originalPolynomial^k if k > currentPower, else the currentPolynomial is not changed.
 * todo: maybe delete current poly as it is not needed. Also, are we deleting the linear forms 2x?
 */
void BoxOptimizationContinuous::setPower(int k)
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
	currentMap.SetLength(k+1);
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



void BoxOptimizationContinuous::printSpolynomial() const
{
	int n= currentMap.length();
	cout << "Final spoly ";
	for(int i = 0; i < n; ++i)
		if (currentMap[i] >= 0)
			cout <<  " + " << currentMap[i] << "*(s^" << i <<") ";
		else
			cout <<  "  " << currentMap[i] << "*(s^" << i <<") ";
	cout << endl;

}

/**
 * integrate the polynomial.
 */
void BoxOptimizationContinuous::findSPolynomial(AlgoType algoType, const vec_RRorQQ &lowerBound, const vec_RRorQQ & upperBound)
{

	BTrieIterator<RationalNTL, int>* itr = new BTrieIterator<RationalNTL, int>();
    itr->setTrie(currentPolynomial.myMonomials, currentPolynomial.varCount);
    itr->begin();
	term<RationalNTL, int> * t;

	clear(currentMap); //length unchanged.

	for(t = itr->nextTerm(); t; t = itr->nextTerm())
	{
		//cout << "integral of " << t->coef;
		//for (int i = 1; i <currentPolynomial.varCount; ++i)
		//	cout << "*x[" << i << "]^" << t->exps[i];
		//cout << " is ";

		RRorQQ value;
		value = 1;

		for(int i = 1; i < currentPolynomial.varCount; ++i)
			if  ( lowerBound[i-1] + zero < upperBound[i-1] )
				value *= (power(upperBound[i-1], (long)(t->exps[i] + 1)) - power(lowerBound[i-1], (long)(t->exps[i] + 1)))/(t->exps[i] + 1);
			else
				value *= power(upperBound[i-1], (long) t->exps[i]);
		value *= t->coef.to_RR();

		//cout << value << endl;
		currentMap[t->exps[0]] += value;

	}



	delete itr;

}

int BoxOptimizationContinuous::findRange(int nitr)
{

	RR newU, newL;
	bool improved = true;
	int improvedOnce = 0;
	for(int i = 0; i < nitr && improved; ++i)
	{
		improved = false;
		newU = maximumUpperbound();
		if ( newU < U)
		{
			//cout << "newU=" <<newU <<" < " << U << "=U" << endl;
			U = newU;
			improved = true;
			improvedOnce |= 0x1;
		}


		int d = originalPolynomial.varCount -1;
		RR p1, p2;
		p1 = d;
		p1 /= d+currentPower;
		p2 = currentPower;
		p2 /= d+currentPower;

		//p1 = d/(d+k)
		//p2 = k/(d+k)

		RR f1, f2, f3;

		int sign = 1;
		if (currentPower % 2 )
			sign = -1;

		if ( p1 * (U - fmaxLowerBound) / (M*lipschitz) < 1) // p1 *fmax/ML < 1
		{


			f1 = pow(evalSpoly(-U)*sign/V, inv(to_RR(d+currentPower)));
			f2 = pow( (M*lipschitz * (d+currentPower)) * inv(to_RR(d)), p1);
			f3 = pow( to_RR(d+currentPower)/to_RR(currentPower), p2);
			newL = U -  f1*f2*f3;
		}
		else
		{
			newL = L;
		}

		//f1 = pow(evalSpoly(-U)*sign/V, inv(to_RR(d+currentPower)));
		//f2 = pow( (M*lipschitz * (d+currentPower)) * inv(to_RR(d)), p1);
		//f3 = pow( to_RR(d+currentPower)/to_RR(currentPower), p2);
		//newL = U - f1*f2*f3;

		if ( L < newL)
		{
			//cout << "L=" << L <<" < " << newL << "=newL" << endl;
			L = newL;
			improved = true;
			improvedOnce |= 0x2;
			//cout << "L improved" <<endl;
		}
	}

	return improvedOnce;
}


//Evaluate the s-polynomial .
RR BoxOptimizationContinuous::evalSpoly(const RR & s) const
{
	//ZZ zs = to_ZZ(s);
	//RationalNTL qs(zs, 1);
	//RationalNTL sPow(qs);
	//RationalNTL ans;
	RR ans;
	RR sPow = s;
	int p = 1;
	if ( currentMap.length() == 0)
		return ans;
	ans = currentMap[0];


	for (p = 1; p < currentMap.length() ; ++p)
	{
		//cout << "p=" << p << endl;
		//cout << " adding " << currentMap[p] << " * " << sPow <<   "(P= " << p<< endl;
		ans += sPow*currentMap[p];
		sPow *= s;
	}

	return ans;

}




