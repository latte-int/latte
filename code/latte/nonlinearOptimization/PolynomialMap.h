/*
 * PolynomialMap.h
 *
 *  Created on: Feb 3, 2015
 *      Author: bedutra
 */

#ifndef POLYNOMIALMAP_H_
#define POLYNOMIALMAP_H_
#include <map>
#include <ostream>
#include "rational.h"


class PolynomialMap
{
public:
	std::map<int, RationalNTL> terms;
	PolynomialMap & operator+=(const PolynomialMap &rhs);
	bool  operator==(const int rhs);
	friend std::ostream & operator<<(std::ostream& os, const PolynomialMap & rhs);
	RR eval(const RR & s) const;
	void mult(const RationalNTL &rhs);
	void print(ostream &out);
	void add(const PolynomialMap & p, const RationalNTL c); //adds p*s
};

#endif /* POLYNOMIALMAP_H_ */
