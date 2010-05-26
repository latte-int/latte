//will define multivariate polynomial representation, as well as input and output functions
#ifndef POLYREP_H
#define POLYREP_H

#include <NTL/vec_ZZ.h>
#include <NTL/ZZ.h>
#include "consumers.h"
#include "blockReps.h"
#include "iterators.h"

NTL_CLIENT

void _loadMonomials(_monomialSum&, const string&);
//string parsing
void _parseMonomials(_MonomialConsumer<ZZ>*, const string&);
//data structure operations
string _printMonomials(const _monomialSum&);
void _destroyMonomials(_monomialSum&);

void _loadLinForms(_linFormSum&, const string);
//string parsing
void _parseLinForms(_FormSumConsumer<ZZ>*, const string&);
//data structure operations
string _printLinForms(const _linFormSum&);
void _destroyLinForms(_linFormSum&);

void _decompose(_monomialSum&, _linFormSum&, int);
#endif
