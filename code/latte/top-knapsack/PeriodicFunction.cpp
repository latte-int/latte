/*
 * PeriodicFunction.cpp
 *
 *  Created on: Jun 5, 2013
 *      Author: bedutra
 */

#include "PeriodicFunction.h"



PeriodicFunction::~PeriodicFunction()
{
	if (head)
		delete head;
}

PeriodicFunctionNode::~PeriodicFunctionNode()
{
	if(! isLeaf)
	{
		if(left)
			delete left;
		if(right)
			delete right;
	}
}

PeriodicFunctionNode::PeriodicFunctionNode(): isLeaf(true), isNumber(true), left(NULL), right(NULL)
{

}

PeriodicFunctionNode::PeriodicFunctionNode(Operation operation, PeriodicFunctionNode * l, PeriodicFunctionNode * r):
		isLeaf(false), isNumber(false), opt(operation), left(l), right(r)
{
}

/**
 * Recursivly call the copy constructor to make a deep copy of this tree.
 */
PeriodicFunctionNode::PeriodicFunctionNode(const PeriodicFunctionNode& p):
		isLeaf(p.isLeaf), isNumber(p.isNumber), data(p.data), opt(p.opt)
{
	if (p.left)
		left = new PeriodicFunctionNode(*(p.left));
	else
		left = NULL;

	if(p.right)
		right = new PeriodicFunctionNode(*(p.right));
	else
		right = NULL;

}

PeriodicFunctionNode::PeriodicFunctionNode(const RationalNTL & d, bool isN):
		isLeaf(true), isNumber(isN), data(d), left(NULL), right(NULL)
{

}


PeriodicFunction::PeriodicFunction(): head(NULL)
{}

/*
 * Removes ownership of any memory.

void PeriodicFunction::clear()
{
	head = NULL;
}
 */

void PeriodicFunction::setToConstant(int c)
{
	if (head)
		delete head;
	head = new PeriodicFunctionNode(RationalNTL(c,1), true);
}

/*
 * This object becomes the owner of the memory in p, and p is "cleared"
 */
void PeriodicFunction::add( PeriodicFunctionNode * p)
{
	if ( head ==  NULL)
	{
		head = p;
	}
	else
	{
		PeriodicFunctionNode * temp = new PeriodicFunctionNode(PeriodicFunctionNode::plus, head, p);
		head = temp;
	}
}


void PeriodicFunction::addProduct(const RationalNTL &coeff, const RationalNTL & function)
{
	//check if we can add these numbers.
	if(coeff.getNumerator() == 0 || function.getDenominator() == 1)
		return ; //coefficient or periodic function is zero so skip it.

	//temp is a plus note with two leaf children: one a number, the other a function
	PeriodicFunctionNode *temp;
	temp = new PeriodicFunctionNode(PeriodicFunctionNode::times,
			new PeriodicFunctionNode(coeff, true),
			new PeriodicFunctionNode(function, false));

	//if there is a head, we have to make a new plus node.
	if (head == NULL)
		head = temp;
	else
	{
		head = new PeriodicFunctionNode(PeriodicFunctionNode::plus, head, temp);
	}

}

PeriodicFunction & PeriodicFunction::operator=(const PeriodicFunction & p)
{
	if( this == &(p))
		return *this;

	if(head)
		delete head;
	if(p.head == NULL)
	{
		head = NULL;
		return *this;
	}
	head = new PeriodicFunctionNode(*(p.head));

	return *this;
}

void PeriodicFunction::pow(int p)
{

	if ( p == 0)
	{
		setToConstant(1); //num^0 = 1, 0^0=1
		return;
	}
	if ( !head)
		return; //0^p =0 for p!= 0
	if ( p == 1)
		return;


	head = new PeriodicFunctionNode(PeriodicFunctionNode::power,
			head,  new PeriodicFunctionNode(RationalNTL(p,1), true));
}

void PeriodicFunction::div(const ZZ & d)
{
	if(d == 1)
		return; //don't waste time dividing by 1.
	if ( !head)
		return; //0/d  = 0....not checking of d=0.
	head = new PeriodicFunctionNode(PeriodicFunctionNode::divide,
			head, new PeriodicFunctionNode(RationalNTL(d,1), true));
}

void PeriodicFunction::times(const RationalNTL & d)
{
	if (!head)
		return; //0*d=0
	if (head->isLeaf && head->isNumber)
		head->data *= d;
	else
	{
		head = new PeriodicFunctionNode(PeriodicFunctionNode::times,
				head, new PeriodicFunctionNode(d, true));
	}

}

bool PeriodicFunction::operator==(const int x) const
{
	if(!head && x == 0)
		return true;
	if (! head)
		return false;
	return (head->isLeaf && head->isNumber && head->data == x);
}

PeriodicFunction & PeriodicFunction::operator+=(const PeriodicFunction &p)
{
	if(p.head == NULL)
		return *this; //if adding the zero function

	PeriodicFunctionNode * copy = new PeriodicFunctionNode(*(p.head));

	if(head == NULL)
	{
		head = copy;
		return *this;
	}//if this was the zero function.

	PeriodicFunctionNode * temp = new PeriodicFunctionNode(PeriodicFunctionNode::plus, head, copy);

	head = temp;
	return *this;
}

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
	if(isLeaf)
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
ostream& operator<<(ostream& out, const PeriodicFunction & pf)
{

	if (pf.head != NULL)
		out << *(pf.head);
	else
		out << 0;
	return out;
}

ostream& operator<<(ostream& out, const PeriodicFunctionNode & pfn)
{

	if(pfn.isLeaf)
	{
		if(pfn.isNumber)
			out << "(" << pfn.data << ")";
		else
			out << "( FF(" << pfn.data << ", T) )";
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




