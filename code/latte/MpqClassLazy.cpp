/*
 * mpqclasslazy.cpp
 *
 *  Created on: Nov 25, 2014
 *      Author: bedutra
 */

#include "MpqClassLazy.h"
#include "latte_gmp.h"
#include <iostream>


MpqLazy::MpqLazy()
{
	mpz_init_set_si(numerator,0);
	mpz_init_set_si(denominator,1);
	funCall = 0;

}
MpqLazy::~MpqLazy()
{
	mpz_clear(numerator);
	mpz_clear(denominator);
	funCall = 0;
}

MpqLazy::MpqLazy(const mpq_class & r)
{
	mpz_init(numerator);
	mpz_init(denominator);

	mpz_set(numerator, r.get_num_mpz_t());
	mpz_set(denominator, r.get_den_mpz_t());
	funCall = 0;
}

MpqLazy::MpqLazy(const MpqLazy & rhs)
{
	mpz_init_set(numerator, rhs.numerator);
	mpz_init_set(denominator, rhs.denominator);
	funCall = rhs.funCall;
}

MpqLazy::MpqLazy(const mpz_class &num, const mpz_class & denom)
{
	mpz_init(numerator);
	mpz_init(denominator);

	mpz_set(numerator, num.get_mpz_t());
	mpz_set(denominator, denom.get_mpz_t());
	funCall = 0;
}

MpqLazy::MpqLazy(const mpz_t &num, const mpz_t & denom)
{
	mpz_init(numerator);
	mpz_init(denominator);

	mpz_set(numerator, num);
	mpz_set(denominator, denom);
	funCall = 0;
}


MpqLazy & MpqLazy::tryCanonicalize()
{
	++funCall;
	if (funCall > 3)
	{
		funCall = 0;
		return canonicalize();
	}
	return *this;
}


MpqLazy & MpqLazy::canonicalize()
{
	mpz_t gcd;
	mpz_init(gcd);

	mpz_gcd(gcd, numerator, denominator);

	if(mpz_sgn(gcd) > 0)
	{
		mpz_divexact(numerator, numerator, gcd);
		mpz_divexact(denominator, denominator, gcd);
	}

	if ( mpz_sgn(denominator) < 0)
	{
		mpz_neg(numerator, numerator);
		mpz_neg(denominator, denominator);
	}
	mpz_clear(gcd);

	return *this;
}

//ADDITION
MpqLazy MpqLazy::operator+(const MpqLazy & rhs) const
{
	//a/b + c/e = (ea + bc)/be
	MpqLazy ans;
	mpz_t op1,op2;
	mpz_init(op1);
	mpz_init(op2);

	mpz_mul(op1, numerator, rhs.denominator);
	mpz_mul(op2, denominator, rhs.numerator);
	mpz_add(op1, op1, op2);
	mpz_set(ans.numerator, op1);

	mpz_mul(op2, denominator, rhs.denominator);
	mpz_set(ans.denominator, op2);

	mpz_clear(op1);
	mpz_clear(op2);
	ans.funCall = funCall + rhs.funCall;
	return ans.tryCanonicalize();
}

MpqLazy MpqLazy::operator-(const MpqLazy & rhs) const
{
	//a/b - c/e = (ea - bc)/be
	MpqLazy ans;
	mpz_t op1,op2;
	mpz_init(op1);
	mpz_init(op2);

	mpz_mul(op1, numerator, rhs.denominator);
	mpz_mul(op2, denominator, rhs.numerator);
	mpz_sub(op1, op1, op2);
	mpz_set(ans.numerator, op1);

	mpz_mul(op2, denominator, rhs.denominator);
	mpz_set(ans.denominator, op2);

	mpz_clear(op1);
	mpz_clear(op2);
	ans.funCall = funCall +rhs.funCall;
	return ans;
}

MpqLazy & MpqLazy::operator+=(const MpqLazy &rhs)
{
	//a/b + c/e = (ea + bc)/be
	mpz_t op1;
	mpz_mul(numerator, numerator, rhs.denominator);
	mpz_init(op1);
	mpz_mul(op1, denominator, rhs.numerator);
	mpz_add(numerator, numerator, op1);
	mpz_mul(denominator, denominator, rhs.denominator);
	funCall += rhs.funCall;
	tryCanonicalize();
	return *this;
}


MpqLazy & MpqLazy::operator-=(const MpqLazy &rhs)
{
	//a/b - c/e = (ea - bc)/be
	mpz_t op1;
	mpz_mul(numerator, numerator, rhs.denominator);
	mpz_init(op1);
	mpz_mul(op1, denominator, rhs.numerator);
	mpz_sub(numerator, numerator, op1);
	mpz_mul(denominator, denominator, rhs.denominator);
	funCall += rhs.funCall;
	tryCanonicalize();
	return *this;
}

MpqLazy MpqLazy::operator/(const MpqLazy & rhs) const
{
	MpqLazy ans(*this);
	mpz_mul(ans.numerator, ans.numerator, rhs.denominator);
	mpz_mul(ans.denominator, ans.denominator, rhs.numerator);
	ans.funCall += rhs.funCall;
	return ans.tryCanonicalize();
}

MpqLazy operator/(int a, const MpqLazy & rhs)
{
	MpqLazy ans(rhs);
	mpz_mul_si(ans.numerator, ans.numerator, a);
	ans.funCall += rhs.funCall;
	return ans.tryCanonicalize();
}
MpqLazy & MpqLazy::operator/=(const MpqLazy & rhs)
{
	mpz_mul(numerator, numerator, rhs.denominator);
	mpz_mul(denominator, denominator, rhs.numerator);
	funCall += rhs.funCall;
	return tryCanonicalize();
}

MpqLazy & MpqLazy::operator/=(const mpz_class & rhs)
{
	mpz_mul(denominator, denominator, rhs.get_mpz_t());
	return tryCanonicalize();
}

MpqLazy MpqLazy::operator*(const MpqLazy & rhs) const
{
	MpqLazy ans(*this);
	mpz_mul(ans.numerator, numerator, rhs.numerator);
	mpz_mul(ans.denominator, denominator, rhs.denominator);
	ans.funCall += rhs.funCall;
	return ans.tryCanonicalize();
}
MpqLazy MpqLazy::operator*(const mpz_class & rhs) const
{
	MpqLazy ans(*this);
	mpz_mul(ans.numerator, ans.numerator, rhs.get_mpz_t());
	return ans.tryCanonicalize();
}

MpqLazy & MpqLazy::operator*=(const MpqLazy & rhs)
{
	mpz_mul(numerator, numerator, rhs.numerator);
	mpz_mul(denominator, denominator, rhs.denominator);
	funCall += rhs.funCall;
	return tryCanonicalize();
}

MpqLazy & MpqLazy::operator*=(const mpz_class & hrs)
{
	mpz_mul(numerator, numerator, hrs.get_mpz_t());
	return tryCanonicalize();
}

MpqLazy MpqLazy::operator-() const
{
	MpqLazy ans(*this);
	mpz_neg(ans.numerator, ans.numerator);
	return ans.tryCanonicalize();
}


MpqLazy & MpqLazy::operator=(const MpqLazy & rhs)
{
	mpz_set(numerator, rhs.numerator);
	mpz_set(denominator, rhs.denominator);
	funCall = rhs.funCall;
	return *this;
}
MpqLazy & MpqLazy::operator=(const int rhs)
{
	mpz_set_si(numerator, rhs);
	mpz_set_si(denominator, 1);
	funCall = 0;
	return *this;
}


ostream& operator <<(ostream &out, const MpqLazy & rhs)
{
	mpz_class n(rhs.numerator), d(rhs.denominator);
	out << n << "/" << d;
	return out;
}

mpq_class MpqLazy::to_mpq()
{
	mpq_class ans;
	ans = mpz_class(numerator);
	ans /= mpz_class(denominator);
	return ans;
}
const mpz_t & MpqLazy::getNumerator() const {return numerator; }
const mpz_t & MpqLazy::getDenominator() const {return denominator; }


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

MpqClassLazy::MpqClassLazy() {
	numerator = 0;
	denominator = 1;
}

MpqClassLazy::~MpqClassLazy() {
	// TODO Auto-generated destructor stub
}

MpqClassLazy::MpqClassLazy(const MpqClassLazy & rhs)
{
	numerator = rhs.numerator;
	denominator = rhs.denominator;
}

MpqClassLazy::MpqClassLazy(const mpq_class & r)
{
	numerator = r.get_num();
	denominator = r.get_den();
}


MpqClassLazy::MpqClassLazy(const mpz_class &num, const mpz_class & denom): numerator(num), denominator(denom)
{

}

MpqClassLazy & MpqClassLazy::tryCanonicalize()
{
	if (mpz_sizeinbase(numerator.get_mpz_t(), 2) > 50 || mpz_sizeinbase(denominator.get_mpz_t(), 2) > 50 )
	{
		cout << "try canonicalize is true! " << endl;
		cout << "canonicalize called on " << numerator << "/" << denominator << endl;
		return canonicalize();
	}
	return *this;
}

MpqClassLazy & MpqClassLazy::canonicalize()
{
	//cout << "canonicalize called on " << numerator << "/" << denominator << endl;
	mpz_t g;
	mpz_init(g);
	mpz_gcd(g, numerator.get_mpz_t(), denominator.get_mpz_t());
	mpz_class G(g);
	mpz_clear(g);
	numerator /= G;
	denominator /= G;
	if (denominator < 0)
	{
		denominator *= -1;
		numerator *= -1;
	}
	return *this;
}

MpqClassLazy MpqClassLazy::operator+(const MpqClassLazy & rhs) const
{
	//a/b + c/e = (ea + bc)/be
	MpqClassLazy ans(rhs.denominator * numerator + denominator*rhs.numerator, denominator*rhs.denominator);
	ans.tryCanonicalize();
	return ans;
}

MpqClassLazy MpqClassLazy::operator-(const MpqClassLazy & rhs) const
{
	//a/b - c/e = (ea - bc)/be
	MpqClassLazy ans(rhs.denominator * numerator - denominator*rhs.numerator, denominator*rhs.denominator);
	ans.tryCanonicalize();
	return ans;
}

MpqClassLazy & MpqClassLazy::operator+=(const MpqClassLazy &rhs)
{
	if (rhs.numerator == 0)
		return *this;
	//a/b + c/e = (ea + bc)/be
	numerator = rhs.denominator * numerator + denominator*rhs.numerator;
	denominator *= rhs.denominator;
	return tryCanonicalize();
}

MpqClassLazy & MpqClassLazy::operator-=(const MpqClassLazy &rhs)
{
	if (rhs.numerator == 0)
		return *this;
	//a/b - c/e = (ea - bc)/be
	numerator = rhs.denominator * numerator - denominator*rhs.numerator;
	denominator *= rhs.denominator;
	return tryCanonicalize();
}

MpqClassLazy MpqClassLazy::operator/(const MpqClassLazy & rhs) const
{
	if (numerator == 0)
		return *this;
	//a/b div c/e = ae/bc
	MpqClassLazy ans(numerator*rhs.denominator, denominator*rhs.numerator);
	ans.tryCanonicalize();
	return ans;
}

MpqClassLazy operator/(int a, const MpqClassLazy & rhs)
{
	// a div c/e = ae/c;
	MpqClassLazy ans(a*rhs.denominator, rhs.numerator);
	ans.tryCanonicalize();
	return ans;
}

MpqClassLazy & MpqClassLazy::operator/=(const MpqClassLazy & rhs)
{
	if (numerator == 0)
		return *this;
	numerator *= rhs.denominator;
	denominator *= rhs.numerator;
	return tryCanonicalize();
}

MpqClassLazy & MpqClassLazy::operator/=(const mpz_class & rhs)
{
	if (numerator == 0)
		return *this;
	denominator *= rhs;
	return tryCanonicalize();
}

MpqClassLazy MpqClassLazy::operator*(const MpqClassLazy & rhs) const
{
	if (numerator == 0)
		return *this;
	if ( rhs.numerator == 0)
		return rhs;
	//a/b * c/e = ac/be
	MpqClassLazy ans(numerator*rhs.numerator, denominator*rhs.denominator);
	ans.tryCanonicalize();
	return ans;
}

MpqClassLazy MpqClassLazy::operator*(const mpz_class & rhs) const
{
	if (numerator == 0)
		return *this;
	MpqClassLazy ans(numerator*rhs, denominator);
	ans.tryCanonicalize();
	return ans;
}

MpqClassLazy & MpqClassLazy::operator*=(const MpqClassLazy & rhs)
{
	if (numerator == 0)
		return *this;
	numerator *= rhs.numerator;
	denominator *= rhs.denominator;
	return tryCanonicalize();
}

MpqClassLazy & MpqClassLazy::operator*=(const mpz_class & rhs)
{
	if (numerator == 0)
		return *this;
	numerator *= rhs;
	if ( rhs == 0)
		denominator = 1;
	return tryCanonicalize();
}

MpqClassLazy MpqClassLazy::operator-() const
{
	MpqClassLazy a(*this);
	a.numerator *= -1;
	return a;
}


MpqClassLazy & MpqClassLazy::operator=(const MpqClassLazy & rhs)
{
	numerator = rhs.numerator;
	denominator = rhs.denominator;
	return tryCanonicalize();
}

MpqClassLazy & MpqClassLazy::operator=(const int rhs)
{
	numerator = rhs;
	denominator = 1;
	return *this;
}

ostream & operator<<(ostream & out, const MpqClassLazy & rhs)
{
	out << rhs.numerator << "/" << rhs.denominator;
	return out;
}


mpq_class MpqClassLazy::to_mpq()
{
	canonicalize();
	return mpq_class(numerator, denominator);
}

const mpz_class & MpqClassLazy::getNumerator() const
{
	return numerator;
}
const mpz_class & MpqClassLazy::getDenominator() const
{
	return denominator;
}

