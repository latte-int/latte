#include <NTL/vec_ZZ.h>
#include <NTL/ZZ.h>
#include "PolyTrie.h"
#include "multiply.h"

ZZ Power_ZZ(ZZ a, int b) //power function computes a^b
{
	if (b == 0)
		return to_ZZ(1);
	int bi[20];
	int digit = 0;
	while (b > 0)
	{
		digit++;
		bi[digit - 1] = b % 2;
		b = b / 2;
	};
	ZZ t = a;
	for (int i = digit - 2; i >= 0; i--)
	{
		t *= t;
		if (bi[i] == 1)
			t *= a;
	};
	return t;
}

ZZ AChooseB(int a, int b)
{
	ZZ t = to_ZZ(1);
	if (b > a)
		return to_ZZ(0);
	if (2 * b > a)
		b = a - b;
	for (int i = 1; i <= b; i++)
	{
		t = t * (a - i + 1) / i;
	};
	return t;
}
//This function is called when a vertex is irregular. The function computes the residue at the irregular vertex
void computeResidue(int d, int M, const vec_ZZ &innerProDiff, const ZZ &p,
		ZZ &a, ZZ &b)
{
	if (p == 0)
	{
		a = 0;
		b = 1;
		return;
	}; //vertex vanishes, return 0;
	int k, i, j;
	int counter[1000];//counter counts number of appearances of each index[i]
	vec_ZZ index;//collecting different terms in the innerProDiff passed in
	bool found;
	ZZ de, nu, g; //c was ZZ
	RationalNTL c; //coefficient? --Brandon
	int e[1];
	int mindeg[1];
	int maxdeg[1];
	k = 1;
	index.SetLength(d);
	index[k - 1] = 0;
	counter[k - 1] = 0;
	for (i = 0; i <= d; i++)
	{
		found = 0;
		for (j = 0; j < k; j++)
		{
			if (innerProDiff[i] == index[j])
			{
				counter[j]++;
				found = 1;
				break;
			}
		};
		if (!found)
		{
			k++;
			index[k - 1] = innerProDiff[i];
			counter[k - 1] = 1;
		};
	};
	counter[0]--; //excluding the vertex itself
	//so far I've been doing book keeping stuff, i.e. how many different differences of each kind are there? index stores the differences and counter stores the multiplicity. 0th entry means 0 itself. It is very important because it's the multiplicity of a vertex itself.
	//actual calculations, I want the appropriate coefficient in a product of one polynomial and some power series. (everything is truncated).
	nu = 1;
	de = 1;
	//for (i=1;i<=counter[0];i++) nu*=i;
	for (j = 1; j <= k - 1; j++)
		de = de * Power_ZZ(index[j], counter[j] + counter[0]);
	monomialSum m1;
	monomialSum sub; //sub is the substitution for m1, which alternatively stores the product for each other
	m1.varCount = 1;
	m1.termCount = 0;
	sub.varCount = 1;
	sub.termCount = 0;
	for (i = 0; i <= counter[0]; i++)
	{
		c = AChooseB(M + d, i) * Power_ZZ(p, M + d - i);
		e[0] = i;
		insertMonomial(c, e, m1);
	};
	for (i = 1; i < k; i++)
	{
		monomialSum m2;
		m2.varCount = 1;
		m2.termCount = 0;
		for (j = 0; j <= counter[0]; j++)
		{
			c = AChooseB(counter[i] + j - 1, j) * Power_ZZ(index[i], counter[0]
					- j);
			if (j % 2 == 1)
				c.mult(to_ZZ(-1), to_ZZ(1));
			e[0] = j;
			insertMonomial(c, e, m2);
		};
		mindeg[0] = 0;
		maxdeg[0] = counter[0];
		BTrieIterator<RationalNTL, int>* it = new BTrieIterator<RationalNTL, int> ();
		BTrieIterator<RationalNTL, int>* it2 = new BTrieIterator<RationalNTL, int> ();
		if (i % 2 == 1)
		{
			it->setTrie(m1.myMonomials, m1.varCount);
			it2->setTrie(m2.myMonomials, m2.varCount);
			multiply<RationalNTL> (it, it2, sub, mindeg, maxdeg);
		}//cout<<"times "<<printMonomials(m2)<<" gives "<<printMonomials(sub)<<endl;}
		else
		{
			it->setTrie(sub.myMonomials, sub.varCount);
			it2->setTrie(m2.myMonomials, m2.varCount);
			multiply<RationalNTL> (it, it2, m1, mindeg, maxdeg);
		}//cout<<"times "<<printMonomials(m2)<<"gives "<<printMonomials(m1)<<endl;};
	};
	//ZZ findCoeff = to_ZZ(0);
	RationalNTL findCoeff;
	/*    This part does the same thing as the uncommented part following this, but using the original block structure for multiplication
	 if (k % 2) //m1
	 {
	 //search m1 for the first term whose exponent vector is equal to [counter[0]]
	 }
	 else //sub
	 {
	 //search sub for the first term whose exponent vector is equal to [counter[0]]
	 }
	 eBlock* myExps; cBlock<ZZ>* myCoeffs;
	 if (k % 2==1) 					//choose which one to pick result from
	 {myExps = m1.eHead; myCoeffs = m1.cHead;
	 for (i=0;i<m1.termCount;i++)
	 {
	 if (i>0 && i % BLOCK_SIZE ==0)
	 {
	 myExps = myExps->next; myCoeffs=myCoeffs->next;
	 };
	 if (myExps->data[i % BLOCK_SIZE]== counter[0]) {findCoeff=myCoeffs->data[i % BLOCK_SIZE];break;};
	 };
	 }
	 else
	 {myExps = sub.eHead; myCoeffs = sub.cHead;
	 for (i=0;i<sub.termCount;i++)
	 {
	 if (i>0 && i % BLOCK_SIZE ==0)
	 {
	 myExps = myExps->next; myCoeffs=myCoeffs->next;
	 };
	 if (myExps->data[i % BLOCK_SIZE]== counter[0]) {findCoeff=myCoeffs->data[i % BLOCK_SIZE];break;};
	 };
	 };*/

	//The following part is trying to find a monomial that has the degree we want and returns its coefficient to findCoeff
	BurstTerm<RationalNTL, int>* temp;
	BurstTrie<RationalNTL, int>* myTrie;
	BTrieIterator<RationalNTL, int>* it = new BTrieIterator<RationalNTL, int> ();
	if (k % 2 == 1)
	{
		temp = new BurstTerm<RationalNTL, int> (m1.varCount);
		myTrie = m1.myMonomials;
		it->setTrie(myTrie, m1.varCount);
	} else
	{
		temp = new BurstTerm<RationalNTL, int> (sub.varCount);
		myTrie = sub.myMonomials;
		it->setTrie(myTrie, sub.varCount);
	}
	it->begin();
	term<RationalNTL, int>* storedTerm;
	while (storedTerm = it->nextTerm())
	{
		if (storedTerm->exps[0] == counter[0])
		{
			findCoeff = storedTerm->coef;
			break;
		};
	};

	a = nu * findCoeff.getNumerator();
	b = de * findCoeff.getDenominator();
	g = GCD(a, b);
	if (g != 0)
	{
		a = a / g;
		b = b / g;
	};
	return;
}
;
