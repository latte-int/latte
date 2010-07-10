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
void multiply(PolyIterator<ZZ, int>* it, PolyIterator<ZZ, int>* it2, monomialSum& result, int* min, int* max)
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
void multiply(PolyIterator<ZZ, int>* it, PolyIterator<ZZ, int>* it2, _monomialSum& result, int* min, int* max)
{
	//result.myMonomials = new BurstTrie<T, int>();
	int* resExps = new int[result.varCount];
	
	it->begin();
	it2->begin();

	term<T, int> *firstTerm = it->nextTerm();
	term<T, int> *secondTerm = it2->nextTerm();
	
	int i;
	ZZ prod;
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


template <class T>
void _multiply(_monomialSum& first, _monomialSum& second, _monomialSum& result, int* min, int* max)
{
	//cout << "Old multiply" << endl;
	eBlock* firstExp = first.eHead; cBlock<T>* firstCoef = first.cHead;
	eBlock* secondExp; cBlock<T>* secondCoef;
	
	int* exponents = new int[result.varCount];
	bool valid;
	result.termCount = 0;
	for (int i = 0; i < first.termCount; i++)
	{
		if (i > 0 && i % BLOCK_SIZE == 0) //this block is done, get next one
		{
			firstExp = firstExp->next; firstCoef = firstCoef->next;
		}
		secondExp = second.eHead;
		secondCoef = second.cHead;
		for (int j = 0; j < second.termCount; j++)
		{
			if (j > 0 && j % BLOCK_SIZE == 0) //this block is done, get next one
			{
				secondExp = secondExp->next; secondCoef = secondCoef->next;
			}
			valid = true;
			for (int k = 0; k < result.varCount; k++)
			{
				exponents[k] = firstExp->data[(i % BLOCK_SIZE)*first.varCount + k] + secondExp->data[(j % BLOCK_SIZE)*second.varCount + k];
				//if (exponents[k] < min[k] || exponents[k] > max[k]) {valid = false; break; }
			}
			if (valid) //all exponents are within range
			{
				_insertMonomial<T>(firstCoef->data[i % BLOCK_SIZE] * secondCoef->data[j % BLOCK_SIZE], exponents, result);
			}
		}
		
	}
	delete [] exponents;
}
#endif
