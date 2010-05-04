#ifndef MULTIPLY_H
#define MULTIPLY_H
#include "PolyTrie.h"
#include "PolyRep.h"

#define mult_DEBUG 0

// Multipies two monomial sums, storing the result in the third one
// Any values stored in result will be overwritten
// result is every term in the product of two monomial sums whose exponents are greater than min and lower than max.
// min, max point to int arrays of length result.varCount
template <class T>
void multiply(monomialSum& first, monomialSum& second, monomialSum& result, int* min, int* max)
{
	result.myMonomials = new BurstTrie<T, int>();
	int* resExps = new int[result.varCount];
	
	BurstTerm<T, int> *firstTerm, *secondTerm;
	
	BTrieIterator<ZZ, int>* it = new BTrieIterator<ZZ, int>();
	it->setTrie(first.myMonomials, first.varCount);
	it->begin();
	
	BTrieIterator<ZZ, int>* it2 = new BTrieIterator<ZZ, int>();
	it2->setTrie(second.myMonomials, second.varCount);
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
	
	/*
	 
	
	//wlog, first has more terms than second
	if (second.termCount > first.termCount)
	{ return multiply<T>(second, first, result, min, max); }
	
	cout << first.termCount << " and " << second.termCount << " terms" << endl;
	
	result.myMonomials = new BurstTrie<T, int>();
	int* resExps1 = new int[result.varCount];
	
	BurstTerm<T, int> *fTerm1, *fTerm2;
	
	if (second.termCount == 1)
	{
		BTrieIterator<T, int>* it = new BTrieIterator<T, int>();
		it->setTrie(first.myMonomials, first.varCount);
		it->begin();
		fTerm1 = it->nextTerm();
		
		BTrieIterator<T, int>* it2 = new BTrieIterator<T, int>();
		it2->setTrie(second.myMonomials, second.varCount);
		it2->begin();
		fTerm2 = it2->nextTerm();
		
		int i;
		while (fTerm1)
		{
			for (i = 0; i < result.varCount; i++)
			{
				resExps1[i] = fTerm1->exps[i] + fTerm2->exps[i];
				if (resExps1[i] < min[i] || resExps1[i] > max[i]) { break; }
			}
			if (i == result.varCount)
			{ result.myMonomials->insertTerm(fTerm1->coef * fTerm2->coef, resExps1, 0, result.varCount, -1); }
			fTerm1 = it->nextTerm();
		}
		
	}
	else
	{
		fTerm1 = fTerm2 = NULL;
		
		int* resExps2 = new int[result.varCount];
		bool valid1, valid2;
		int state;
		
		BTrieIterator<T, int>* it = new BTrieIterator<T, int>();
		it->setTrie(second.myMonomials, second.varCount);
		it->begin();
		
		BTrieIterator<T, int>** iterators = new BTrieIterator<T, int>*[second.termCount];
		BurstTerm<T, int> ** iteratorTerms = new BurstTerm<T, int> *[second.termCount];
		int itIndex = 0;
		for (int i = 0; i < second.termCount; i++)
		{
			iterators[i] = new BTrieIterator<T, int>();
			iterators[i]->setTrie(first.myMonomials, first.varCount);
			iterators[i]->begin();
			
			iteratorTerms[i] = new BurstTerm<T, int>(result.varCount);
			BurstTerm<T, int>* temp = it->nextTerm();
			iteratorTerms[i]->coef = temp->coef;
			memcpy(iteratorTerms[i]->exps, temp->exps, sizeof(temp->exps));
		}
		delete it;
		
		while (itIndex < second.termCount)
		{
			
		}
	}
	
	
	
	
	curTerm = it->nextTerm(); 
	nextTerm = it->nextTerm();
	
	if (!nextTerm) //only one term
	{
		cout << "one term" << endl;
		BTrieIterator<ZZ, int>* it2 = new BTrieIterator<T, int>();
		it2->setTrie(first.myMonomials, first.varCount);
		it2->begin();
		nextTerm = it2->nextTerm();
		do
		{
			int i;
			for (i = 0; i < result.varCount; i++)
			{
				res1 = nextTerm->exps[i] + curTerm->exps[i];
				if (res1 < min[i] || res1 > max[i]) { break; }
				resExps[i] = res1;
				cout << nextTerm->exps[i] << ", " << curTerm->exps[i] << "; ";
			}
			cout << endl;
			if (i == result.varCount)
			{
				result.myMonomials->insertTerm(nextTerm->coef * curTerm->coef, resExps, 0, result.varCount, -1);
			}
			nextTerm = it2->nextTerm();
		}
		while (nextTerm);
		
		delete it2;
		delete it;
		delete resExps;
		return;
	}
	
	BTrieIterator<T, int>** iterators = new BTrieIterator<T, int>*[second.termCount];
	int itIndex = 0;
	for (int i = 0; i < second.termCount; i++)
	{
		iterators[i] = new BTrieIterator<T, int>();
		iterators[i]->setTrie(first.myMonomials, first.varCount);
		iterators[i]->begin();
	}
	
	while (itIndex < second.termCount)
	{
		if (!fTerm1)
		{ fTerm1 = iterators[itIndex]->nextTerm(); }
		
		if (fTerm1)
		{
			if (!fTerm2)
			{ fTerm2 = iterators[(itIndex + 1) % second.termCount]->nextTerm(); }
			
			int i;
			for (i = 0; i < result.varCount; i++)
			{
				res1 = fTerm1->exps[i] + curTerm->exps[i];
				res2 = fTerm2->exps[i] + nextTerm->exps[i];
				if (res1 < res2)
				{
					resExps[i] = res1;
					break;
				}
				else if (res2 < res1)
				{
					resExps[i] = res2;
					break;
				}
				else //equal
				{
					resExps[i] = res1;
				}
			}
			
			if (res1 < res2 || i == result.varCount) //don't advance down
			{
				for (; i < result.varCount; i++)
				{
					resExps[i] = fTerm1->exps[i] + curTerm->exps[i];
				}
				//insert fTerm1 * curTerm
				result.myMonomials->insertTerm(fTerm1->coef * curTerm->coef, resExps, 0, result.varCount, -1);
				fTerm1 = NULL;
			}
			else
			{
				for (; i < result.varCount; i++)
				{
					resExps[i] = fTerm1->exps[i] + curTerm->exps[i];
				}
				//insert fTerm1 * curTerm
				result.myMonomials->insertTerm(fTerm1->coef * curTerm->coef, resExps, 0, result.varCount, -1);
				for (int j = 0; j < result.varCount; j++)
				{
					resExps[j] = fTerm2->exps[j] + nextTerm->exps[j];
				}
				//insert fTerm2 * nextTerm
				result.myMonomials->insertTerm(fTerm2->coef * nextTerm->coef, resExps, 0, result.varCount, -1);
				fTerm1 = fTerm2 = NULL;
				itIndex++; //advancing down
				if (itIndex == second.termCount) { itIndex = 0; it->begin(); }
			}
		}
		else
		{
			fTerm1 = fTerm2;
			itIndex++;
			curTerm = nextTerm;
			nextTerm = it->nextTerm();
		}
	}*/
}


template <class T>
void _multiply(_monomialSum& first, _monomialSum& second, _monomialSum& result, int* min, int* max)
{
	cout << "Old multiply" << endl;
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
