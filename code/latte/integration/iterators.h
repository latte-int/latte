#ifndef ITERATORS_H
#define ITERATORS_H

#include "burstTrie.h"
#include "blockReps.h"
#include <NTL/vec_ZZ.h>
#include <NTL/ZZ.h>

NTL_CLIENT

//Generic struct used by base iterator
template <class T, class S>
struct term
{
	T coef;
	S* exps;
	int length;
	int degree;
};

template <class T, class S>
class PolyIterator
{
public:
	virtual void begin() = 0;
	virtual term<T, S>* nextTerm() = 0;
	virtual term<T, S>* getTerm() = 0;
};

#include "iterators.hpp"

#endif
