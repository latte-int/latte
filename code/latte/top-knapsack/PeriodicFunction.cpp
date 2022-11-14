/*
 * PeriodicFunction.cpp
 *
 *  Created on: Jun 5, 2013
 *      Author: bedutra
 */

#include "PeriodicFunction.h"


#if defined(HAVE_STD_SHARED_PTR) || defined(HAVE_STD_TR1_SHARED_PTR)
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
	return !(left || right);
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
 * @parm d is assumed to be non-zero.
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
 * @parm d we are not checking of d is zero or one.
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
			out << "( MOD( t * (" << pfn.data << "), 1 ) )";
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

#else

PeriodicFunction::~PeriodicFunction()
{
	if (head)
		delete head;
}

PeriodicFunctionNode::~PeriodicFunctionNode()
{
	if(! isLeaf())
	{
		if(left)
			delete left;
		if(right)
			delete right;
	}
}

PeriodicFunctionNode::PeriodicFunctionNode(): isNumber(true), left(NULL), right(NULL)
{

}

PeriodicFunctionNode::PeriodicFunctionNode(Operation operation, PeriodicFunctionNode * l, PeriodicFunctionNode * r):
		isNumber(false), opt(operation), left(l), right(r)
{
}

/**
 * Recursivly call the copy constructor to make a deep copy of this tree.
 */
PeriodicFunctionNode::PeriodicFunctionNode(const PeriodicFunctionNode& p):
		isNumber(p.isNumber), data(p.data), opt(p.opt)
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
		isNumber(isN), data(d), left(NULL), right(NULL)
{

}

bool PeriodicFunctionNode::isLeaf() const
{
	return (left == NULL && right == NULL);
}


PeriodicFunction::PeriodicFunction(const PeriodicFunction & p)
{
	if ( p.head)
		head = new PeriodicFunctionNode(*(p.head));
	else
		head = NULL;
}

PeriodicFunction::PeriodicFunction(const RationalNTL & d, bool type)
{
	head = new PeriodicFunctionNode(d, type);
}

PeriodicFunction::PeriodicFunction(): head(NULL)
{}



void PeriodicFunction::setToConstant(int c)
{
	if (head)
		delete head;
	head = new PeriodicFunctionNode(RationalNTL(c,1), true);
}


void PeriodicFunction::setToConstant(const RationalNTL & c)
{
	if (head)
		delete head;
	head = new PeriodicFunctionNode(c, true);
}


void PeriodicFunction::add(const PeriodicFunction &p)
{
	if ( p.head == NULL)
		return;
	else if (head == NULL)
	{
		*this = p;
	}
	else if ( head->isNumber && p.head->isNumber)
	{
		head->data += p.head->data;
	}
	else
		head = new PeriodicFunctionNode(PeriodicFunctionNode::plus, head, new PeriodicFunctionNode(*(p.head)) );
}

void PeriodicFunction::subtract(const PeriodicFunction & p)
{
	if ( p.head == NULL)
		return;
	else if (head == NULL)
	{
		*this = p;
	}
	else if ( head->isNumber && p.head->isNumber)
		head->data -= p.head->data;
	else
		head = new PeriodicFunctionNode(PeriodicFunctionNode::minus, head, new PeriodicFunctionNode(*(p.head)) );
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

PeriodicFunction & PeriodicFunction::operator=(const int c)
{

	if(head)
		delete head;

	head = new PeriodicFunctionNode(RationalNTL(c,1), true);

	return *this;
}

void PeriodicFunction::pow(int p)
{

	if ( p == 0)
	{
		setToConstant(1); //num^0 = 1, 0^0=1
	}
	else if ( !head)
		return; //0^p =0 for p!= 0
	else if ( p == 1)
		return;
	else if ( head->isNumber)
		head->data.power(p);
	else
	head = new PeriodicFunctionNode(PeriodicFunctionNode::power,
			head,  new PeriodicFunctionNode(RationalNTL(p,1), true));
}

void PeriodicFunction::div(const ZZ & d)
{
	if(d == 1)
		return; //don't waste time dividing by 1.
	if ( !head)
		return; //0/d  = 0....not checking of d=0.
	if (head->isNumber)
		head->data.div(d);
	else
		head = new PeriodicFunctionNode(PeriodicFunctionNode::divide,
			head, new PeriodicFunctionNode(RationalNTL(d,1), true));
}

void PeriodicFunction::times(const RationalNTL & d)
{
	if (!head)
		return; //0*d=0
	if (head->isNumber)
		head->data *= d;
	else
	{
		head = new PeriodicFunctionNode(PeriodicFunctionNode::times,
				head, new PeriodicFunctionNode(d, true));
	}

}

void PeriodicFunction::times(const PeriodicFunction& p)
{
	if(!head)
		return; //0*p=0
	if(head->isLeaf() && head->isNumber && p.head->isLeaf() && p.head->isNumber)
		head->data *= p.head->data;
	else
	{
		head = new PeriodicFunctionNode(PeriodicFunctionNode::times,
				head, new PeriodicFunctionNode(*(p.head)));
	}
}

void PeriodicFunction::operator*=(const PeriodicFunction & d)
{
	times(d);
}

bool PeriodicFunction::operator==(const int x) const
{
	if(!head && x == 0)
		return true;
	if (! head)
		return false;
	return (head->isLeaf() && head->isNumber && head->data == x);
}

void PeriodicFunction::operator+=(const PeriodicFunction &p)
{
	add(p);
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

#endif

