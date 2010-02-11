#ifndef MULTIPLY_H
#define MULTIPLY_H
#include "PolyRep.h"

// Multipies two monomial sums, storing the result in the third one
// Any values stored in result will be overwritten
// result is every term in the product of two monomial sums whose exponents are greater than min and lower than max.
// min, max point to int arrays of length result.varCount
template <class T>
void multiply(monomialSum& first, monomialSum& second, monomialSum& result, int* min, int* max)
{
	if (first.termCount == 0 || second.termCount == 0) { cout << "Only one monomial sum given, aborting."; return; }
	
	eBlock* firstExp = first.eHead; cBlock<T>* firstCoef = first.cHead;
	eBlock* secondExp = second.eHead; cBlock<T>* secondCoef = second.cHead;
	
	int* exponents = new int[result.varCount];
	bool valid;
	result.termCount = 0;
	for (int i = 0; i < first.termCount; i++)
	{
		if (i > 0 && i % BLOCK_SIZE == 0) //this block is done, get next one
		{
			firstExp = firstExp->next; firstCoef = firstCoef->next;
		}
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
				if (exponents[k] < min[k] || exponents[k] > max[k]) {valid = false; break; }
			}
			if (valid) //all exponents are within range
			{
				insertMonomial<T>(firstCoef->data[i % BLOCK_SIZE] * secondCoef->data[j % BLOCK_SIZE], exponents, result);
			}
		}
		
	}
	delete [] exponents;
}

#endif
