#ifndef MULTIPLY_H
#define MULTIPLY_H
#include "iterators.h"
#include "PolyTrie.h"
#include "PolyRep.h"

#define mult_DEBUG 0

// Multipies two monomial sums, storing the result in the third one
// Any values stored in result will be overwritten
// result is every term in the product of two monomial sums whose exponents are greater than min and lower than max.
// min, max point to int arrays of length result.varCount
template <class T>
void multiply(PolyIterator<RationalNTL, int>* it, PolyIterator<RationalNTL, int>* it2, monomialSum& result, int* min, int* max)
{
	result.myMonomials = new BurstTrie<T, int>();
	int* resExps = new int[result.varCount];
	
	term<T, int> *firstTerm, *secondTerm;
	
	it->begin();
	it2->begin();
	
	int i;
	while (firstTerm = it->nextTerm())
	{
		while (secondTerm = it2->nextTerm())
		{
			for (i = 0; i < result.varCount; i++)
			{
				resExps[i] = firstTerm->exps[i] + secondTerm->exps[i];
				if (resExps[i] < min[i] || resExps[i] > max[i]) { break; }
			}
			
			if (i == result.varCount)
			{
				result.myMonomials->insertTerm(firstTerm->coef * secondTerm->coef, resExps, 0, result.varCount, -1);
			}
		}
		it2->begin();
	}
}

template <class T>
void multiply(PolyIterator<RationalNTL, int>* it, PolyIterator<RationalNTL, int>* it2, _monomialSum& result, int* min, int* max)
{
	//result.myMonomials = new BurstTrie<T, int>();
	int* resExps = new int[result.varCount];
	
	it->begin();
	it2->begin();

	term<T, int> *firstTerm = it->nextTerm();
	term<T, int> *secondTerm = it2->nextTerm();
	
	int i;
	RationalNTL prod;
	while (firstTerm)
	{
		while (secondTerm)
		{
			for (i = 0; i < result.varCount; i++)
			{
				resExps[i] = firstTerm->exps[i] + secondTerm->exps[i];
				if (resExps[i] < min[i] || resExps[i] > max[i]) { break; }
			}
			
			if (i == result.varCount)
			{
				prod = firstTerm->coef * secondTerm->coef;
				_insertMonomial(prod, resExps, result);
			}
			secondTerm = it2->nextTerm();
		}
		it2->begin();
		firstTerm = it->nextTerm();
	}
}
#endif
