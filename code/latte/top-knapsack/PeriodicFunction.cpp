/*
 * PeriodicFunction.cpp
 *
 *  Created on: Jun 5, 2013
 *      Author: bedutra
 */

#include "PeriodicFunction.h"



/**
 * Set the shared pointer to null.
 */
PeriodicFunction::~PeriodicFunction()
{
	head.reset();
}

/**
 * Set the shared pointer to null.
 */
PeriodicFunctionNode::~PeriodicFunctionNode()
{
	left.reset();
	right.reset();
}

/**
 * Set to constant function 0.
 */
PeriodicFunctionNode::PeriodicFunctionNode(): isNumber(true)//, left(nullptr), right(nullptr)
{
}

PeriodicFunctionNode::PeriodicFunctionNode(Operation operation, PeriodicFunctionNodePtr l, PeriodicFunctionNodePtr r):
		isNumber(false), opt(operation), left(l), right(r)
{
}

/**
 * Copy class values, but make a pointer copy of the left and right subtrees.
 */
PeriodicFunctionNode::PeriodicFunctionNode(const PeriodicFunctionNode& p):
		isNumber(p.isNumber), data(p.data), opt(p.opt)
{
	left = p.left;
	right = p.right;
}

PeriodicFunctionNode::PeriodicFunctionNode(const RationalNTL & d, bool isN):
		isNumber(isN), data(d)
{
}

bool PeriodicFunctionNode::isLeaf() const
{
	return (left == NULL && right == NULL);
}

/**
 * Make a pointer copy.
 */
PeriodicFunction::PeriodicFunction(const PeriodicFunction & p)
{
	head = p.head;
}

PeriodicFunction::PeriodicFunction(const RationalNTL & d, bool isN)
{
	head = PeriodicFunctionNodePtr(new PeriodicFunctionNode(d, isN));
}

/**
 * Set to the constant function 0
 */
PeriodicFunction::PeriodicFunction()
{
	head = PeriodicFunctionNodePtr(new PeriodicFunctionNode( RationalNTL(0,1), true));
}



void PeriodicFunction::setToConstant(int c)
{
	head = PeriodicFunctionNodePtr(new PeriodicFunctionNode(RationalNTL(c,1), true));
}


void PeriodicFunction::setToConstant(const RationalNTL & c)
{
	head = PeriodicFunctionNodePtr(new PeriodicFunctionNode(c, true));
}


/**
 * If this and p are constant functions, simply add the numbers, else update the tree.
 */
void PeriodicFunction::add(const PeriodicFunction &p)
{
	if (head->isLeaf() && head->isNumber && (p.head)->isLeaf() && (p.head)->isNumber )
	{
		head = PeriodicFunctionNodePtr(new PeriodicFunctionNode(head->data + (p.head)->data,true));
	}
	else
	{
		head = PeriodicFunctionNodePtr(new PeriodicFunctionNode(PeriodicFunctionNode::plus, head, p.head));
	}
}

/**
 * If this and p are constant functions, simply add the numbers, else update the tree.
 */
void PeriodicFunction::subtract(const PeriodicFunction & p)
{
	if (head->isLeaf() && head->isNumber && (p.head)->isLeaf() && (p.head)->isNumber )
	{
		head = PeriodicFunctionNodePtr(new PeriodicFunctionNode(head->data - (p.head)->data,true));
	}
	else
	{
		head = PeriodicFunctionNodePtr(new PeriodicFunctionNode(PeriodicFunctionNode::minus, head, p.head) );
	}
}

/**
 * Make a pointer copy.
 */
PeriodicFunction & PeriodicFunction::operator=(const PeriodicFunction & p)
{
	if( this == &(p))
		return *this;

	head = p.head;

	return *this;
}

PeriodicFunction & PeriodicFunction::operator=(const int c)
{
	setToConstant(c);
	return *this;
}

/**
 * Anything to the zero power is 1, even 0^0.
 * If this is a number, the perform the power.
 */
void PeriodicFunction::pow(int p)
{

	if ( p == 0)
	{
		setToConstant(1); //num^0 = 1, 0^0=1
	}
	else if ( head->isLeaf() && head->isNumber)
		head = PeriodicFunctionNodePtr(new PeriodicFunctionNode(RationalNTL::power(head->data, p),true));
	else
		head = PeriodicFunctionNodePtr(new PeriodicFunctionNode(PeriodicFunctionNode::power, head,  PeriodicFunctionNodePtr(new PeriodicFunctionNode(RationalNTL(p,1), true))));
}

/**
 * @param d is assumed to be non-zero.
 */
void PeriodicFunction::div(const ZZ & d)
{
	if(d == 1)
		return; //don't waste time dividing by 1.
	if ( head->isNumber && head->isLeaf() )
		head = PeriodicFunctionNodePtr(new PeriodicFunctionNode( (head->data)/d ,true));
	else
		head = PeriodicFunctionNodePtr(new PeriodicFunctionNode(PeriodicFunctionNode::divide, head, PeriodicFunctionNodePtr(new PeriodicFunctionNode(RationalNTL(d,1), true))));
}

/**
 * @param d we are not checking of d is zero or one.
 */
void PeriodicFunction::times(const RationalNTL & d)
{
	if (head->isNumber && head->isLeaf())
		head = PeriodicFunctionNodePtr(new PeriodicFunctionNode( (head->data)*d ,true));
	else
	{
		head = PeriodicFunctionNodePtr(new PeriodicFunctionNode(PeriodicFunctionNode::times,
				head, PeriodicFunctionNodePtr(new PeriodicFunctionNode(d, true))));
	}

}

void PeriodicFunction::times(const PeriodicFunction& p)
{
	if(head->isLeaf() && head->isNumber && (p.head)->isLeaf() && (p.head)->isNumber)
		head = PeriodicFunctionNodePtr(new PeriodicFunctionNode( (head->data) * ((p.head)->data) ,true));
	else
	{
		head = PeriodicFunctionNodePtr(new PeriodicFunctionNode(PeriodicFunctionNode::times,
				head, p.head));
	}
}

void PeriodicFunction::operator*=(const PeriodicFunction & d)
{
	times(d);
}

/**
 * @return true if this is a constant function equal to x.
 */
bool PeriodicFunction::operator==(const int x) const
{
	return (head->isLeaf() && head->isNumber && head->data == x);
}

void PeriodicFunction::operator+=(const PeriodicFunction &p)
{
	add(p);
}

/**
 * Mostly used for debugging.
 */
void PeriodicFunction::print() const
{
	if(head)
		head->print(0);
}

void PeriodicFunctionNode::print(int i) const
{
	string s;
	for(int j = 0; j < i; ++j)
		s+="  ";

	cout << s.c_str() << "node level " << i << endl;
	if(isLeaf())
	{
		if(isNumber)
			cout << s.c_str() <<" num " << data << endl;
		else
			cout <<s.c_str() << " fun " << data << endl;
	}
	else
	{
		cout << s.c_str() <<" node " << opt << endl;
		cout << s.c_str() <<" left:" << endl;
		if(left)
			left->print(i+1);
		cout << s.c_str() <<" right:" << endl;
		if(right)
			right->print(i+1);
	}
}

/**
 * The main print function
 */
ostream& operator<<(ostream& out, const PeriodicFunction & pf)
{
	if (pf.head)
		out << *(pf.head);
	return out;
}

ostream& operator<<(ostream& out, const PeriodicFunctionNode & pfn)
{
	if(pfn.isLeaf())
	{
		if(pfn.isNumber)
			out << "(" << pfn.data << ")";
		else
			out << "( Frac( T * (" << pfn.data << ") ) )";
	}
	else
	{
		out << "(" << *(pfn.left);

		switch(pfn.opt) //plus, minus, times, divide, power
		{
		case PeriodicFunctionNode::plus: out << " + "; break;
		case PeriodicFunctionNode::minus: out << " - "; break;
		case PeriodicFunctionNode::times: out << "*"; break;
		case PeriodicFunctionNode::divide: out << "/"; break;
		case PeriodicFunctionNode::power: out << "^"; break;
		}

		out << *(pfn.right) << ")";
	}

	return out;
}



