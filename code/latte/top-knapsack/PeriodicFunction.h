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
#include "config.h"






class PeriodicFunction;
class PeriodicFunctionNode;

/**
 * PeriodicFunction.h: Classes to store periodic functions in the form {x}:=x-floor(x).
 * We do not evaluate expressions like {4/5 T}^10 + 4{4/5 T}^10 to 5{4/5 T}^10.
 * Instead, everything is saved in an operator tree using shared pointers.
 */

#undef HAVE_STD_SHARED_PTR
#undef HAVE_STD_TR1_SHARED_PTR

#if defined(HAVE_STD_SHARED_PTR)
#include <memory>
typedef std::shared_ptr<PeriodicFunctionNode> PeriodicFunctionNodePtr;
#elif defined(HAVE_STD_TR1_SHARED_PTR)
#include <tr1/memory>
typedef std::tr1::shared_ptr<PeriodicFunctionNode> PeriodicFunctionNodePtr;
#endif





/**
 * Internal node for the PeriodicFunction binary tree.
 * Shared pointers are used to save memory and construction time.
 */

class PeriodicFunctionNode
{
private:

	bool isNumber;    //!< if true, data should be thought of as just a fraction, else data represents the function {a/b*T}
	RationalNTL data; //!< a fraction or argument to a periodice function
	enum Operation { plus, minus, times, divide, power}; //!< set of poperators that this node can be.
	Operation opt;    //!< operator for this node. this is only defined if isLeaf() is false.

#if defined(HAVE_STD_SHARED_PTR) || defined(HAVE_STD_TR1_SHARED_PTR)
	PeriodicFunctionNodePtr left, right; //!< If not a leaf, this this node is an binary operator (+,-,*,/,^) and the left and right children are some expressions.
#else
	PeriodicFunctionNode *left, *right;
#endif
	friend class PeriodicFunction; //!< only the PeriodicFunction class should use this node class.


	PeriodicFunctionNode();
#if defined(HAVE_STD_SHARED_PTR) || defined(HAVE_STD_TR1_SHARED_PTR)
	PeriodicFunctionNode(Operation operation, PeriodicFunctionNodePtr  l, PeriodicFunctionNodePtr  r);
#else
	PeriodicFunctionNode(Operation operation, PeriodicFunctionNode * l, PeriodicFunctionNode * r);
#endif
	PeriodicFunctionNode(const PeriodicFunctionNode& p);
	PeriodicFunctionNode(const RationalNTL & d, bool isN);

	bool isLeaf() const; //!< true if this node is a leaf in the tree, representing a number or function.

	void print(int i) const; //!< for debugging.
	friend ostream& operator<<(ostream& out, const PeriodicFunctionNode & pfn);
public:
~PeriodicFunctionNode();
};

/**
 * PeriodicFunction is a class to hold an expression of periodic functions in the form {x}:=x-floor(x).
 *
 *
 * We make a binary tree of nodes and operators that connect them. This is best said with a graphic
 *
 *
 *                                             root(head): non-leaf, operation +
 *               left: non-leaf, operation power                           right: non-leaf: operation times
 *     (left: leaf. function {a/b * T},  right: leaf. num 10)          (left: leaf function {c/d*T},     right: leaf number 20
 *
 * This example tree holds the expression ({a/b*T}^10) + ({c/d*T}*20)
 */
class PeriodicFunction
{
private:
#if defined(HAVE_STD_SHARED_PTR) || defined(HAVE_STD_TR1_SHARED_PTR)
	PeriodicFunctionNodePtr  head; //!< pointer to the root of the tree.
#else
	PeriodicFunctionNode * head;
#endif
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


	PeriodicFunction & operator=(const PeriodicFunction & p);
	PeriodicFunction & operator=(const int c);
	bool operator==(const int) const;
	void operator+=(const PeriodicFunction &p);
	void operator*=(const PeriodicFunction & d);

	void print() const;
	friend ostream& operator<<(ostream& out, const PeriodicFunction & pf);
};


#endif /* PERIODICFUNCTION_H_ */
