//will define multivariate polynomial representation, as well as input and output functions
#ifndef POLYREP_H
#define POLYREP_H

#define BLOCK_SIZE 64
#include <NTL/vec_ZZ.h>
#include <NTL/ZZ.h>

NTL_CLIENT

struct eBlock //used to store exponent vector
{
	eBlock* next;
	vec_ZZ* data; //each vec_ZZ contains n elements for an n-variable polynomial
};

struct lBlock //used to store linear forms (coefficient vector and degree)
{
	lBlock* next;
	vec_ZZ* data; //each vec_ZZ contains n elements for an n-variable polynomial
	ZZ* degree; //for linear forms only
};

struct cBlock
{
	cBlock* next;
	ZZ* data;
};

struct monomialSum
{
	int termCount, varCount;
	eBlock* eHead; //variable exponents
	cBlock* cHead; //monomial coefficients
};

struct linFormSum
{
	int termCount, varCount;
	lBlock* lHead; //linear forms
	cBlock* cHead; //linear form coefficients
};

void loadMonomials(monomialSum&, const string);
string printMonomials(const monomialSum&);
void destroyMonomials(monomialSum&);

void loadLinForms(linFormSum&, const string);
string printLinForms(const linFormSum&);
void destroyLinForms(linFormSum&);

void decompose(monomialSum&, linFormSum&, int);
ZZ Power(const ZZ&a, const ZZ& e);
#endif
