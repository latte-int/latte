/*
 * mpqclasslazy.h
 *
 *  Created on: Nov 25, 2014
 *      Author: bedutra
 */

#ifndef MPQCLASSLAZY_H_
#define MPQCLASSLAZY_H_

#include "latte_gmp.h"
#include <vector>
#include <iostream>
#include <fstream>
using namespace std;

class MpqLazy
{
private:
	mpz_t numerator, denominator;
	int funCall;
	//MpqClassLazy & canonicalize();
public:

	// CONSTRUCTORS
	MpqLazy(); //initialize to 0
	~MpqLazy();
	MpqLazy(const mpq_class & r);
	MpqLazy(const MpqLazy & rhs);
	MpqLazy(const mpz_class &num, const mpz_class & denom);
	MpqLazy(const mpz_t &num, const mpz_t & denom);

	MpqLazy & tryCanonicalize();
	MpqLazy & canonicalize(); // reduces the fraction to lowest terms, and makes
						//	the denominator positive if it can be.

	//ADDITION
	MpqLazy operator+(const MpqLazy & rhs) const;

	MpqLazy operator-(const MpqLazy & rhs) const;
	MpqLazy & operator+=(const MpqLazy &rhs);
	MpqLazy & operator-=(const MpqLazy &rhs);

	//DIVISION
	MpqLazy operator/(const MpqLazy & rhs) const;
	friend MpqLazy operator/(int a, const MpqLazy & rhs);
	MpqLazy & operator/=(const MpqLazy & rhs) ;
	MpqLazy & operator/=(const mpz_class & rhs);

	//MULTIPLICATION
	MpqLazy operator*(const MpqLazy & rhs) const;
	MpqLazy operator*(const mpz_class & rhs) const;
	MpqLazy & operator*=(const MpqLazy & hrs);
	MpqLazy & operator*=(const mpz_class & hrs);
	MpqLazy operator-() const;


	//ASSIGNMENT
	MpqLazy & operator=(const MpqLazy & rhs);
	MpqLazy & operator=(const int rhs);

	//OUT
	friend ostream& operator <<(ostream &out, const MpqLazy & rhs);
	mpq_class to_mpq();
	const mpz_t & getNumerator() const;
	const mpz_t & getDenominator() const;

};


/**
 * MpqClassLazy is a wrapper for mpz_class fractions.
 * Operations are performed in a lazy way, meaning, the fraction is not automatically placed in simplest terms.
 */
class MpqClassLazy
{
private:
	mpz_class numerator, denominator;
	//MpqClassLazy & canonicalize();
public:

	// CONSTRUCTORS
	MpqClassLazy(); //initialize to 0
	~MpqClassLazy();
	MpqClassLazy(const MpqClassLazy & rhs);
	MpqClassLazy(const mpq_class & r);
	MpqClassLazy(const mpz_class &num, const mpz_class & denom);

	MpqClassLazy & tryCanonicalize();
	MpqClassLazy & canonicalize(); // reduces the fraction to lowest terms, and makes
						//	the denominator positive if it can be.

	//ADDITION
	MpqClassLazy operator+(const MpqClassLazy & rhs) const;
	MpqClassLazy operator-(const MpqClassLazy & rhs) const;
	MpqClassLazy & operator+=(const MpqClassLazy &rhs);
	MpqClassLazy & operator-=(const MpqClassLazy &rhs);

	//DIVISION
	MpqClassLazy operator/(const MpqClassLazy & rhs) const;
	friend MpqClassLazy operator/(int a, const MpqClassLazy & rhs);
	MpqClassLazy & operator/=(const MpqClassLazy & rhs) ;
	MpqClassLazy & operator/=(const mpz_class & rhs);

	//MULTIPLICATION
	MpqClassLazy operator*(const MpqClassLazy & rhs) const;
	MpqClassLazy operator*(const mpz_class & rhs) const;
	MpqClassLazy & operator*=(const MpqClassLazy & hrs);
	MpqClassLazy & operator*=(const mpz_class & hrs);
	MpqClassLazy operator-() const;


	//ASSIGNMENT
	MpqClassLazy & operator=(const MpqClassLazy & rhs);
	MpqClassLazy & operator=(const int rhs);

	//OUT
	friend ostream& operator <<(ostream &out, const MpqClassLazy & rhs);
	mpq_class to_mpq();
	const mpz_class & getNumerator() const;
	const mpz_class & getDenominator() const;

};

//typedef std::vector<MpqLazy> mpq_lazy_vector;
//typedef MpqLazy mpq_class_lazy;


//typedef std::vector<MpqClassLazy> mpq_lazy_vector;
//typedef MpqClassLazy mpq_class_lazy;

typedef std::vector<mpq_class> mpq_lazy_vector;
typedef mpq_class mpq_class_lazy;




#endif /* MPQCLASSLAZY_H_ */
