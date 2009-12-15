//will define multivariate polynomial representation, as well as input and output functions
#ifndef POLYREP_H
#define POLYREP_H

#define BLOCK_SIZE 64
#include <NTL/vec_ZZ.h>
#include <NTL/ZZ.h>

NTL_CLIENT

struct eBlock //used to store exponents
{
	vec_ZZ data[BLOCK_SIZE]; //each vec_ZZ contains n elements for an n-variable polynomial, max BLOCK_SIZE of these
};

struct cBlock
{
	ZZ data[BLOCK_SIZE];
};

struct polynomial
{
	int termCount, varCount;
	eBlock* exponentBlocks;
	cBlock* coefficientBlocks;
};

void loadPolynomial(polynomial &, const string);
string printPolynomial(const polynomial &);
void destroyPolynomial(polynomial &);

#endif
