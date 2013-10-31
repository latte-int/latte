/*
 * PeriodicFunction.h
 *
 *  Created on: Jun 5, 2013
 *      Author: bedutra
 */

#ifndef PERIODICFUNCTION_H_
#define PERIODICFUNCTION_H_

#include "latte_ntl.h"
#include "rational.h"


class PeriodicFunction;

/**
 * A class to hold a function of periodic functions {x}:=x-floor(x).
 *
 * We make a binary tree of nodes and operators that connect them. This is best said with a graphic
 *
 *
 *                                                root: non-leaf, operation +
 *               left: non-leaf, operation power                           right: non-leaf: operation times
 *     (left: leaf. function {a/b * T}  right: leaf. num 10)          (left: leaf function {c/d*T}     right: leaf number 20
 *
 * This example tree holds the expression ({a/b*T}^10) + ({c/d*T}*20)
 *
 *
 *
 */
class PeriodicFunctionNode
{
private:
	bool isLeaf;      //if false, data has no meaning and is zero.
	bool isNumber;    //if true, data should be thought of as just a fraction.
					  //else, data represents the function {a/b*T}
	RationalNTL data;
	enum Operation { plus, minus, times, divide, power};
	Operation opt;    //this is only defined if isLeaf is false.

	PeriodicFunctionNode *left, *right;

	friend class PeriodicFunction;

public:
	PeriodicFunctionNode();
	~PeriodicFunctionNode();
	PeriodicFunctionNode(Operation operation, PeriodicFunctionNode * l, PeriodicFunctionNode * r);

	PeriodicFunctionNode(const PeriodicFunctionNode& p);
	PeriodicFunctionNode(const RationalNTL & d, bool isN);

	void print(int i) const; //for debugging.
	friend ostream& operator<<(ostream& out, const PeriodicFunctionNode & pfn);
	//~PeriodicFunctionNode();
};


class PeriodicFunction
{
private:
	PeriodicFunctionNode * head;

public:
	PeriodicFunction();
	PeriodicFunction(const RationalNTL & d, bool type);
	PeriodicFunction(const PeriodicFunction & p);
	~PeriodicFunction();

	void add(const PeriodicFunction & p);
	void subtract(const PeriodicFunction & p);
	void pow(int p);
	void div(const ZZ & d);
	void times(const RationalNTL & d);
	void times(const PeriodicFunction& p);
	void setToConstant(int c);
	void setToConstant(const RationalNTL & c);
	//void clear();

	PeriodicFunction & operator=(const PeriodicFunction & p);
	PeriodicFunction & operator=(const int c);
	bool operator==(const int) const;
	void operator+=(const PeriodicFunction &p);
	void operator*=(const PeriodicFunction & d);

	void print() const;
	friend ostream& operator<<(ostream& out, const PeriodicFunction & pf);
};


#endif /* PERIODICFUNCTION_H_ */
