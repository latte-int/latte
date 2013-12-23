/*
 * GeneralMonomialSum.h
 *
 *  Created on: Jun 5, 2013
 *      Author: bedutra
 */

#ifndef GENERALMONOMIALSUM_H_
#define GENERALMONOMIALSUM_H_

#include "burstTrie.h"
#include <sstream>

using namespace std;

/**
 * A class to deal with polynomials without having to manually make BTrieIterator iterators.
 */


/**
 * Assumptions
 * 		S ? I don't know, we only use int for S.
 * 		T: needs to support ==(int), =, <<, *=.
 */
template <class T, class S>
class GeneralMonomialSum
{
public:
	int termCount, varCount;
	BurstTrie<T, S>* myMonomials;

	GeneralMonomialSum(): termCount(0), varCount(0), myMonomials(NULL)
	{
	}//constructor
	~GeneralMonomialSum();

	void setToConstant(T c);
	void insertMonomial(const T & coefficient, S * exponents);
	void destroyMonomials();
	void multiply(const GeneralMonomialSum<T, S>& gms, const int* min, const int* max);
	void add(const GeneralMonomialSum<T, S> & gms);
	void check(); //for error checking. to delete
	string printMonomials() const;
	string printMonomialsX() const;
	string printMonomials(const int *min, const int * max) const;
};




// #######################################################

template <class T, class S>
GeneralMonomialSum<T, S>::~GeneralMonomialSum()
{
	destroyMonomials();
}

template <class T, class S>
void GeneralMonomialSum<T, S>::destroyMonomials()
{
	if (myMonomials)
		delete myMonomials;
	myMonomials = NULL;
	termCount = varCount = 0;

}


template <class T, class S>
void GeneralMonomialSum<T, S>::multiply(const GeneralMonomialSum<T, S>& gms, const int* min, const int* max)
{
	assert(varCount == gms.varCount);
	if(termCount == 0)
		return; // this is the zero polynomial.
	else if ( gms.termCount == 0)
	{
		int vCount = varCount;
		destroyMonomials();
		varCount = vCount;
		return;
	}// gms is the zero polynomial so zero this out.

	BTrieIterator<T, S>*  it1 = new BTrieIterator<T, S>();
	BTrieIterator<T, S>*  it2 = new BTrieIterator<T, S>();
	BurstTrie<T, S>* oldMonomials = myMonomials;

	it1->setTrie(oldMonomials, varCount);
	it2->setTrie(gms.myMonomials, gms.varCount);

	myMonomials = new BurstTrie<T, S>();
	int* resExps = new int[varCount];

	term<T, S>   *firstTerm;
	term<T, S> *secondTerm;

	it1->begin();
	it2->begin();

	int i;
	termCount = 0;
	while (firstTerm = it1->nextTerm())
	{
		while (secondTerm = it2->nextTerm())
		{
			for (i = 0; i < varCount; i++)
			{
				resExps[i] = firstTerm->exps[i] + secondTerm->exps[i];
				if (resExps[i] < min[i] || resExps[i] > max[i]) { break; }
			}

			if (i == varCount)
			{
				T temp;
				temp = firstTerm->coef;
				temp *= secondTerm->coef;
				myMonomials->insertTerm(temp, resExps, 0, varCount, -1);
				termCount++;
			}
		}
		it2->begin();
	}
	delete [] resExps;
	delete oldMonomials;
	delete it1;
	delete it2;

}


template <class T, class S>
void GeneralMonomialSum<T, S>::check()
{
	assert(varCount > 0);
	if ( termCount == 0)
		return; //added a zero polynomial


	BTrieIterator<T, S>*  it1 = new BTrieIterator<T, S>();
	it1->setTrie(myMonomials, varCount);
	term<T, S> *term;

	it1->begin();

	while ( (term = it1->nextTerm()) )
	{
		for(int i = 0; i < varCount; ++i)
			if ( term->exps[i] < 0)
			{
				cout << "ERROR: polynomial with negative power..." << endl;
				exit(1);
			}
	}

	delete it1;
}

template <class T, class S>
void GeneralMonomialSum<T, S>::add(const GeneralMonomialSum<T, S> & gms)
{
	assert(varCount == gms.varCount);
	if ( gms.termCount == 0)
		return; //added a zero polynomial
	else if ( ! myMonomials)
	{
		myMonomials = new BurstTrie<T, S> ();
		termCount = 0;
	}//else, this is currently the zero polynomial.


	BTrieIterator<T, S>*  it1 = new BTrieIterator<T, S>();
	it1->setTrie(gms.myMonomials, gms.varCount);
	term<T, S> *term;

	it1->begin();

	while ( (term = it1->nextTerm()) )
	{
				myMonomials->insertTerm(term->coef, term->exps, 0, varCount, -1);
				termCount++;
	}

	delete it1;
}

template <class T, class S>
void GeneralMonomialSum<T, S>::insertMonomial(const T & coefficient, S * exponents)
{
	if ( coefficient == 0)
		return;


	if (termCount == 0) //need to construct the first burst trie (sorted on the first variable) and first container
	{
		myMonomials = new BurstTrie<T, S> ();
	}

	myMonomials->insertTerm(coefficient, exponents, 0, varCount, -1);

	termCount++;
}


template <class T, class S>
void GeneralMonomialSum<T, S>::setToConstant(T c)
{
	if (myMonomials)
		delete myMonomials;
	termCount = 0;

	T coeff;
	coeff = c;
	S *exp = new S[varCount];

	for(int i = 0; i < varCount; ++i)
		exp[i] = 0;

	insertMonomial(coeff, exp);

	delete [] exp;
}

//Prints a nested list representation of our sum of monomials
//monomial sum: c_{1}*(x_{1}^e_{1}...x_{varCount}^e_{varCount}) + ...
//nested lists: [[c_{1}, [e_{1}, e_{2}, ..., e_{varCount}]], .. ]
template <class T, class S>
string GeneralMonomialSum<T, S>::printMonomials() const
{
	stringstream output(stringstream::in | stringstream::out);

	if (myMonomials)
	{
		BTrieIterator<T, S>* it =
				new BTrieIterator<T, S> ();
		term<T, S>* temp;
		it->setTrie(myMonomials, varCount);
		it->begin();
		temp = it->nextTerm();
		while(temp)
		{
			if (output.str() != "")
			{
				output << ", ";
			}
			output << "[" << temp->coef << ", [";
			for (int j = 0; j < temp->length; j++)
			{
				output << temp->exps[j];
				if (j + 1 < temp->length)
				{
					output << ", ";
				}
			}
			output << "]]";
			temp = it->nextTerm();
		}
		delete it;
	}
	else
	{
		//print out the zero polynomial.
		output << "[ 0, [ 0";
		for(int i = 0; i < varCount-1; ++i)
			output << ", 0";
		output << "]]";
	}
	return "[" + output.str() + "]";
}

//Prints a nested list representation of our sum of monomials
//monomial sum: c_{1}*(x_{1}^e_{1}...x_{varCount}^e_{varCount}) + ...
//nested lists: [[c_{1}, [e_{1}, e_{2}, ..., e_{varCount}]], .. ]
template <class T, class S>
string GeneralMonomialSum<T, S>::printMonomialsX() const
{
	stringstream output(stringstream::in | stringstream::out);

	if (myMonomials)
	{
		BTrieIterator<T, S>* it =
				new BTrieIterator<T, S> ();
		term<T, S>* temp;
		it->setTrie(myMonomials, varCount);
		it->begin();
		temp = it->nextTerm();
		while(temp)
		{
			if (output.str() != "")
			{
				output << " + ";
			}
			output << "(" << temp->coef << ")*";
			for (int j = 0; j < temp->length; j++)
			{
				output << "x[" << j << "]^(" << temp->exps[j] << ")";
				if (j + 1 < temp->length)
				{
					output << "*";
				}
			}
			//output << "]]";
			temp = it->nextTerm();
		}
		delete it;
	}
	else
	{
		//print out the zero polynomial.
		output << "0";
	}
	return output.str();

}

template <class T, class S>
string GeneralMonomialSum<T, S>::printMonomials(const int *min, const int * max) const
{

	stringstream output(stringstream::in | stringstream::out);
	bool isZero = true;

	if ( myMonomials)
	{
		BTrieIterator<T, S>* it =
				new BTrieIterator<T, S> ();
		term<T, S>* temp;
		it->setTrie(myMonomials, varCount);
		it->begin();
		int i = 0;
		while(temp = it->nextTerm())
		{
			for(i = 0; i < varCount; ++i)
				if (temp->exps[i] < min[i] || max[i] < temp->exps[i])
					break;
			if(i != varCount)
				continue;//skip this while loop.

			isZero = false;
			if (output.str() != "")
			{
				output << ", ";
			}
			output << "[" << temp->coef << ", [";
			for (int j = 0; j < temp->length; j++)
			{
				output << temp->exps[j];
				if (j + 1 < temp->length)
				{
					output << ", ";
				}
			}
			output << "]]";
		}
		delete it;
	}

	if(isZero)
	{
		//print out the zero polynomial.
		output << "[ 0, [ 0";
		for(int i = 0; i < varCount-1; ++i)
			output << ", 0";
		output << "]]";
	}
	return "[" + output.str() + "]";
}



#endif /* GENERALMONOMIALSUM_H_ */
