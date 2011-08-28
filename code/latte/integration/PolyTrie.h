//will define multivariate polynomial representation, as well as input and output functions
#ifndef POLYTRIE_H
#define POLYTRIE_H

#include "iterators.h"
#include "consumers.h"
#include "burstTrie.h"
#include <NTL/ZZ.h>
#include <NTL/vec_ZZ.h>

NTL_CLIENT

void loadMonomials(monomialSum&, const string&);
//string parsing
void parseMonomials(MonomialConsumer<RationalNTL>*, const string&);
//data structure operations
void insertMonomial(const RationalNTL&, int*, monomialSum&);
string printMonomials(const monomialSum&);
void destroyMonomials(monomialSum&);

void loadLinForms(linFormSum&, const string);
//string parsing
void parseLinForms(FormSumConsumer<RationalNTL>*, const string&);
//data structure operations
void insertLinForm(const RationalNTL& coef, int degree, const vec_ZZ& coeffs, linFormSum&);
string printLinForms(const linFormSum&);
void destroyLinForms(linFormSum&);

void loadLinFormProducts(linFormProductSum &forms, const string line);
//string parsing
void parseLinFormProducts(FormProductLoadConsumer<RationalNTL>* consumer, const string& line);
//data structure operations
void destroyLinFormProducts(linFormProductSum &myProd);
string printLinFormProducts(const linFormProductSum &plf);


void decompose(term<RationalNTL, int>*, linFormSum&);
void decompose(BTrieIterator<RationalNTL, int>* it, linFormSum&);
#endif
