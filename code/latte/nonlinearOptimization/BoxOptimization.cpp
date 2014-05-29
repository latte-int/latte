/*
 * BoxOptimizationParameters.cpp
 *
 *  Created on: Apr 23, 2014
 *      Author: bedutra
 */

#include "BoxOptimization.h"
#include "ramon.h"
#include "nonlinearOptimization/WeightedExponentialSubs.h"
#include <map>



	//std::map<int, RationalNTL> terms;
PolynomialMap & PolynomialMap::operator+=(const PolynomialMap &rhs)
{
	for (std::map<int,RationalNTL>::const_iterator it=rhs.terms.begin(); it!=rhs.terms.end(); ++it)
	{

	    terms[it->first] += it->second;
	}
	return *this;
}


bool  PolynomialMap::operator==(const int rhs)
{
	if ( terms.size() != 1)
		return false;

	std::map<int,RationalNTL>::const_iterator it = terms.find(0);
	if (it != terms.end() && it->second == rhs)
	    return true;
	return false;
}

void PolynomialMap::mult(const RationalNTL &rhs)
{
	for (std::map<int,RationalNTL>::iterator it=terms.begin(); it!=terms.end(); ++it)
	{

	    it->second *= rhs;
	}
}

void PolynomialMap::print(ostream &out)
{
	for (std::map<int,RationalNTL>::const_iterator it=terms.begin(); it!=terms.end(); ++it)
	{
		out << it->second << "* S^( " << it->first << " ) + ";
	}
	out << endl;
}


RR PolynomialMap::eval(const RR & s) const
{
	RR ans;
	for (std::map<int,RationalNTL>::const_iterator it=terms.begin(); it!=terms.end(); ++it)
	{
		ans += to_RR((it->second)) * power(s, it->first);
	}
	return ans;
}

std::ostream & operator<<(std::ostream& os, const PolynomialMap & rhs)
{
	os << "(";
	for (std::map<int,RationalNTL>::const_iterator it=rhs.terms.begin(); it!=rhs.terms.end(); ++it)
	{
	    os << " +" << it->second << "*(s^" << it->first << ")";
	}
	os << ")";

	return os;
}

BoxOptimization::BoxOptimization()
{
	U = 0;
	L = 0;
	N = 0;
	currentPower = 0;
	//theTrie = NULL;
}


void BoxOptimization::setPolynomial(const vec_ZZ &lowBound, const vec_ZZ &upBound, const monomialSum & poly)
{
	lowerBound = lowBound;
	upperBound = upBound;
	N = 1;

	vec_ZZ maxBound;
	U = 0;
	maxBound.SetLength(lowBound.length());
	for(int i = 0; i < lowerBound.length(); ++i)
	{
		maxBound[i] = max(abs(lowerBound[i]), abs(upperBound[i]));
		N *= (upperBound[i] - lowerBound[i] + 1);
	}

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
		for (int currentPower = 0; currentPower < poly.varCount; ++currentPower)
		{
			exp[currentPower + 1] = term->exps[currentPower];
			boundTerm *= power(maxBound[currentPower], term->exps[currentPower]);
		}
		insertMonomial(term->coef, exp, originalPolynomial);
		U += to_RR((term->coef) * (boundTerm * sign(term->coef)));
	}
	L = U;
	L *= -1;
	cout << L << " f(x) " << U << endl;

	for (int currentPower = 0; currentPower < poly.varCount; ++currentPower)
		exp[currentPower + 1] = 0;
	exp[0] = 1;
	insertMonomial(RationalNTL(1,1), exp, originalPolynomial);
	cout << "updated poly:" << originalPolynomial.varCount << ", " << originalPolynomial.termCount << endl;
	cout << "updated poly: " << printMonomials(originalPolynomial).c_str() << endl;
	delete [] exp;
}


void BoxOptimization::findRange(int itr)
{


	cout << " Starting L: " << L << endl;
	cout << " Starting U: " << U << endl;

	RR oldU, oldL;
	for(int i = 0; i < itr; ++i)
	{
		oldU = U;
		oldL = L;
		findNewUpperbound();
		findNewLowerbound();

		cout << i << " new U: " << U << endl;
		cout << i << " new L: " << L << endl;

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


	//BTrieIterator<PolynomialMap, ZZ>* xitr = new BTrieIterator<PolynomialMap, ZZ> ();
	//xitr->setTrie(theTrie, lowerBound.length());
	//xitr->begin();

	RR s(L);
	s *= -1;

	//term<PolynomialMap, ZZ>* xitrTerm;
	//RationalNTL one(1,1);
	//for(xitrTerm =  xitr->nextTerm(); xitrTerm; xitrTerm =  xitr->nextTerm())
	//{
		//cout << xitrTerm->coef << "*( ";
		//for(int j = 0; j < lowerBound.length(); ++j)
		//	cout << xitrTerm->exps[j] << ", ";
		//cout << " )^" << xitrTerm->degree << endl;

		//mpq_class temp = computeWeightedCountingBox_singleForm(lowerBound, upperBound, xitrTerm->exps, xitrTerm->degree, one);
		//RR rTemp;
		//rTemp = to_RR(convert_mpz_to_ZZ(temp.get_num()));
		//rTemp /= to_RR(convert_mpz_to_ZZ(temp.get_den()));
		//ans += xitrTerm->coef.eval(s) * rTemp;
//	}
	//delete xitr;
	ans = currentMap.eval(s);

	//cout << "(f+s)^" << currentPower << " ans=" << ans << endl;
	RR newU;
	newU = L + pow(ans, to_RR(1)/to_RR(currentPower));
	if (newU < U)
		U = newU;
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


	//BTrieIterator<PolynomialMap, ZZ>* xitr = new BTrieIterator<PolynomialMap, ZZ> ();
	//xitr->setTrie(theTrie, lowerBound.length());
	//xitr->begin();

	RR s(U);
	s *= -1;

	//term<PolynomialMap, ZZ>* xitrTerm;
	//RationalNTL one(1,1);
	//for(xitrTerm =  xitr->nextTerm(); xitrTerm; xitrTerm =  xitr->nextTerm())
	//{
		//cout << xitrTerm->coef << "*( ";
		//for(int j = 0; j < lowerBound.length(); ++j)
		//	cout << xitrTerm->exps[j] << ", ";
		//cout << " )^" << xitrTerm->degree << endl;

		//mpq_class temp = computeWeightedCountingBox_singleForm(lowerBound, upperBound, xitrTerm->exps, xitrTerm->degree, one);
		//RR rTemp;
		//rTemp = to_RR(convert_mpz_to_ZZ(temp.get_num()));
		//rTemp /= to_RR(convert_mpz_to_ZZ(temp.get_den()));
		//ans += xitrTerm->coef.eval(s) * rTemp;
	//}
	//delete xitr;

	ans = currentMap.eval(s);

	if ( currentPower % 2)
		ans *= -1;

	//cout << "(s-f)^" << currentPower << " ans=" << ans << endl;
	RR newL;
	newL = U - pow(ans, inv(to_RR(currentPower)));
	if ( newL > L)
		L = newL;
}


/**
 * @param k: sets currentPolynomial = originalPolynomial^k if k > currentPower, else the currentPolynomial is not changed.
 */
void BoxOptimization::setPower(int k)
{
	assert(k >= 1);
	//currentPower = k;


	cout << "computing (f(x) + s)^" << k << "..." << flush;
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

	if ( k <= currentPower)
		return; //we only want to increase the power of the currentPolynomial.

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
	}
	currentPower = k;
	//cout << "power poly: " << printMonomials(currentPolynomial).c_str() << endl;
	//cout << "power poly deg: " << k << endl;
	delete it1;
	delete it2;
	cout << "done. \n";



	//****************************************************
	//next, decompose currentPolynomial into linear forms.

	cout << "decomposing (f(x) + s)^" << k << " into powers of linear forms..." << flush;
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
			++numLinForms;
			assert(tempLFTerm->length == originalMonomial.length);
		}
	}
	delete tempLFItr;
	cout << "done. At most " << numLinForms << " decomposed.\n";

	//****************************************************
	//next, integrate the linear forms.

	cout << "Performing the weighted counting..." << flush;
	currentMap.terms.clear();
	BTrieIterator<PolynomialMap, ZZ>* mapitr = new BTrieIterator<PolynomialMap, ZZ> ();
	mapitr->setTrie(theTrie, lowerBound.length());
	mapitr->begin();


	term<PolynomialMap, ZZ>* mapitrTerm;
	RationalNTL one(1,1);
	for(mapitrTerm =  mapitr->nextTerm(); mapitrTerm; mapitrTerm =  mapitr->nextTerm())
	{
		//cout << xitrTerm->coef << "*( ";
		//for(int j = 0; j < lowerBound.length(); ++j)
		//	cout << xitrTerm->exps[j] << ", ";
		//cout << " )^" << xitrTerm->degree << endl;

		mpq_class temp = computeWeightedCountingBox_singleForm(lowerBound, upperBound, mapitrTerm->exps, mapitrTerm->degree, one);
		RationalNTL rTemp;
		rTemp = convert_mpq_to_RationalNTL(temp);
		//cout << rTemp << "*(";
		//mapitrTerm->coef.print(cout);
		mapitrTerm->coef.mult(rTemp);
		//mapitrTerm->coef.print(cout);


		currentMap += mapitrTerm->coef;
		//currentMap.print(cout);
		//cout << endl;
	}
	cout << "done." << endl;
	delete theTrie;
	delete mapitr;

}




RR BoxOptimization::maximumUpperbound() { return U; }

RR BoxOptimization::maximumLowerBound()
{
	if ( currentPower)
		return U/pow(to_RR(N), inv(to_RR(currentPower)));
	return L;
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
	for (lform = linearFormsItr->nextTerm(); lform; lform
			= linearFormsItr->nextTerm())
	{
		ans += computeWeightedCountingBox_singleForm(lowerBound, upperBound, lform->exps, lform->degree, lform->coef);
	}//for every term in the originalPolynomial

	delete linearFormsItr;


	return ans;
}

mpq_class computeWeightedCountingBox_singleForm(const vec_ZZ &lowerBound, const vec_ZZ &upperBound, const ZZ* linFormExps, const int degree, const RationalNTL & coef)
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

	//cout << "newUB" << newUB << endl;
	//cout << "newLB" << newLB << endl;


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
	wbp.Produce(wParams);

	ans = mpq_class(convert_ZZ_to_mpz(coef.getNumerator()), convert_ZZ_to_mpz(coef.getDenominator())) * wParams.result * convert_ZZ_to_mpz(zeroTerms);
	return ans;

}

