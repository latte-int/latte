/*
Defines the BurstTrie class, used for storing monomials and powers of linear forms, as well as iterators over them
*/

#ifndef BURSTTRIE_H
#define BURSTTRIE_H

#define BURST_MAX 2
#include <NTL/ZZ.h>
#include <NTL/vec_ZZ.h>
#include <stdio.h>
#include <sstream>
#include <assert.h>
#include "rational.h"

#define BT_DEBUG 1

NTL_CLIENT

template <class T, class S> class BurstTrie;

struct monomialSum
{
	int termCount, varCount;
	BurstTrie<RationalNTL, int>* myMonomials;

	monomialSum(): termCount(0), varCount(0), myMonomials(NULL)
	{
	}//constructor

};

//linear forms: sort on degree first, then the form coefficients
struct linFormSum
{
	int termCount, varCount;
	BurstTrie<RationalNTL, ZZ>* myForms;

	linFormSum(): 	termCount(0), varCount(0), myForms(NULL)
	{
	}//constructor.

};

#include "burstTrie.hpp"

#endif
