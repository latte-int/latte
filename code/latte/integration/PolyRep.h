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
	int* data; //represents 2D array of ints, storing an exponent vector for each monomial
	//in reality, it's 1D array of length dimension*BLOCK_SIZE
};

struct lBlock //used to store linear forms (coefficient vector and degree)
{
	lBlock* next;
	vec_ZZ* data; //represents 2D array of ints, storing a coefficient vector for each linear form
	//in reality, it's 1D array of length dimension*BLOCK_SIZE
	int degree[BLOCK_SIZE]; //total degree of each linear form
};

struct cBlock
{
	cBlock* next;
	ZZ* data; //this stores the monomial coefficient for each monomial
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
