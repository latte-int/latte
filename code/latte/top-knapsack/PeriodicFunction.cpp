/*
 * PeriodicFunction.cpp
 *
 *  Created on: Jun 5, 2013
 *      Author: bedutra
 */

#include "PeriodicFunction.h"



PeriodicFunction::~PeriodicFunction()
{
	head.reset();	
}

PeriodicFunctionNode::~PeriodicFunctionNode()
{
	left.reset();
	right.reset();
}

PeriodicFunctionNode::PeriodicFunctionNode(): isNumber(true)//, left(nullptr), right(nullptr)
{

}

PeriodicFunctionNode::PeriodicFunctionNode(Operation operation, PeriodicFunctionNodePtr l, PeriodicFunctionNodePtr r):
		isNumber(false), opt(operation), left(l), right(r)
{
}

/**
 * Recursivly call the copy constructor to make a deep copy of this tree.
 */
PeriodicFunctionNode::PeriodicFunctionNode(const PeriodicFunctionNode& p):
		isNumber(p.isNumber), data(p.data), opt(p.opt)
{
	//if (p.left)
	//	left = new PeriodicFunctionNode(*(p.left));
	//else
	//	left = NULL;

	//if(p.right)
	//	right = new PeriodicFunctionNode(*(p.right));
	//else
	//	right = NULL;

	left = p.left;
	right = p.right;
}

PeriodicFunctionNode::PeriodicFunctionNode(const RationalNTL & d, bool isN):
		isNumber(isN), data(d) //, left(nullptr), right(nullptr)
{

}

bool PeriodicFunctionNode::isLeaf() const
{
	return (left == NULL && right == NULL);
}


PeriodicFunction::PeriodicFunction(const PeriodicFunction & p)
{
	//head = new PeriodicFunctionNode(*(p.head));
	head = p.head;
}

PeriodicFunction::PeriodicFunction(const RationalNTL & d, bool isN)
{
	head = PeriodicFunctionNodePtr(new PeriodicFunctionNode(d, isN));
}

PeriodicFunction::PeriodicFunction()
{
	head = PeriodicFunctionNodePtr(new PeriodicFunctionNode( RationalNTL(0,1), true));
}



void PeriodicFunction::setToConstant(int c)
{
	//if (head)
	//	delete head; //this should always be true.
	head = PeriodicFunctionNodePtr(new PeriodicFunctionNode(RationalNTL(c,1), true));
}


void PeriodicFunction::setToConstant(const RationalNTL & c)
{
	//if (head)
	//	delete head;
	head = PeriodicFunctionNodePtr(new PeriodicFunctionNode(c, true));
}


void PeriodicFunction::add(const PeriodicFunction &p)
{

	//static int help = 1;
	//help++;
	if (head->isLeaf() && head->isNumber && (p.head)->isLeaf() && (p.head)->isNumber )
	{
		/*
		RationalNTL temp, temp2;
		temp = head->data;
		temp2 = (p.head)->data;
		if (help == 2658412)
		{
			ZZ a,b,c,d;
			cout << "h" << help << endl;
			a = temp.getNumerator();
			b = temp.getDenominator();

			c = temp2.getNumerator();
			d = temp2.getDenominator();


			cout << "temp: " << temp << endl;
			cout << "      " << a << "/" << b << endl;
			cout << "temp2: " << temp2 << endl;
			cout << "       " << c << "/" << d << endl;

			//a/b + c/d = ad/bd + cb/db
			ZZ num, denom, gcd;
			num = a*d + c*b;
			denom = b*d;

			cout << "num=" << num << endl;
			cout << "denom=" << denom << endl;
			cout << "going to find gcd...." << endl;
			gcd = GCD(num, denom);
			cout << "what?" << endl;
			cout << "gcd = " << gcd << endl;
			num /= gcd;
			denom /= gcd;
			RationalNTL ans(num, denom);
			cout << "ans = " << ans << endl;
			cout << "going to add it in temp" << endl;
			temp = temp + temp2;
			cout << "end of if block" << endl;
		}
		else
			temp = temp + temp2;
		*/
		head = PeriodicFunctionNodePtr(new PeriodicFunctionNode(head->data + (p.head)->data,true));
		//head = PeriodicFunctionNodePtr(new PeriodicFunctionNode(temp,true));
	}
	else
	{
		head = PeriodicFunctionNodePtr(new PeriodicFunctionNode(PeriodicFunctionNode::plus, head, p.head));
	}
}

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
PeriodicFunction & PeriodicFunction::operator=(const PeriodicFunction & p)
{
	if( this == &(p))
		return *this;

	//if(head)
	//	delete head;
	head = p.head;

	return *this;
}

PeriodicFunction & PeriodicFunction::operator=(const int c)
{
	setToConstant(c);
	return *this;
}

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

void PeriodicFunction::div(const ZZ & d)
{
	if(d == 1)
		return; //don't waste time dividing by 1.
	if ( head->isNumber && head->isLeaf() )
		head = PeriodicFunctionNodePtr(new PeriodicFunctionNode( (head->data)/d ,true));
	else
		head = PeriodicFunctionNodePtr(new PeriodicFunctionNode(PeriodicFunctionNode::divide, head, PeriodicFunctionNodePtr(new PeriodicFunctionNode(RationalNTL(d,1), true))));
}

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

bool PeriodicFunction::operator==(const int x) const
{
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




