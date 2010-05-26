#ifndef REPRESENTATIONS_H
#define REPRESENTATIONS_H

#include <NTL/vec_ZZ.h>
#include <NTL/ZZ.h>

NTL_CLIENT

#define BLOCK_SIZE 64

//block-based

//this defines the data structures used
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

template <class T> struct cBlock
{
	cBlock* next;
	T* data; //this stores the monomial coefficient for each monomial
};

struct _monomialSum
{
	int termCount, varCount;
	eBlock* eHead; //variable exponents
	cBlock<ZZ>* cHead; //monomial coefficients over ZZ
};

struct _linFormSum
{
	int termCount, varCount;
	lBlock* lHead; //linear forms
	cBlock<ZZ>* cHead; //linear form coefficients over ZZ
};

#endif
