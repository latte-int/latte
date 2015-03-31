/*
 * PolynomialMap.cpp
 *
 *  Created on: Feb 3, 2015
 *      Author: bedutra
 */

#include "PolynomialMap.h"


	//std::map<int, RationalNTL> terms;
PolynomialMap & PolynomialMap::operator+=(const PolynomialMap &rhs)
{
	for (std::map<int,RationalNTL>::const_iterator it=rhs.terms.begin(); it!=rhs.terms.end(); ++it)
	{

	    terms[it->first] += it->second;
	}
	return *this;
}

void PolynomialMap::add(const PolynomialMap & p, const RationalNTL c)
{
	for (std::map<int,RationalNTL>::const_iterator it=p.terms.begin(); it!=p.terms.end(); ++it)
	{

		    terms[it->first] += c * (it->second);
	}
}


bool  PolynomialMap::operator==(const int rhs)
{
	if ( terms.size() != 1)
		return false;

	std::map<int,RationalNTL>::const_iterator it = terms.find(0);
	if (it != terms.end() && it->second == rhs)
	    return true;
	return false;
}

void PolynomialMap::mult(const RationalNTL &rhs)
{
	for (std::map<int,RationalNTL>::iterator it=terms.begin(); it!=terms.end(); ++it)
	{

	    it->second *= rhs;
	}
}

void PolynomialMap::print(ostream &out)
{
	for (std::map<int,RationalNTL>::const_iterator it=terms.begin(); it!=terms.end(); ++it)
	{
		out << it->second << "* S^( " << it->first << " ) + ";
	}
	out << endl;
}


RR PolynomialMap::eval(const RR & s) const
{
	RR ans;
	for (std::map<int,RationalNTL>::const_iterator it=terms.begin(); it!=terms.end(); ++it)
	{
		ans += to_RR((it->second)) * power(s, it->first);
	}
	return ans;
}

std::ostream & operator<<(std::ostream& os, const PolynomialMap & rhs)
{
	os << "(";
	for (std::map<int,RationalNTL>::const_iterator it=rhs.terms.begin(); it!=rhs.terms.end(); ++it)
	{
		if ( it->second >= 0)
			os << " +" << it->second << "*(s^" << it->first << ")";
		else
			os << " " << it->second << "*(s^" << it->first << ")";
	}
	os << ")";

	return os;
}

