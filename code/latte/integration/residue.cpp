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



/**
 *This function is called when a vertex is irregular. The function computes the residue at the irregular vertex
 * d = dimension
 * m = power of the linear form.
 * innerProDiff[j] =  <l, s_i - s_j>
 * p = <l, s_i>
 * a, b: output parameters: answer: a/b.
 *
 * Formula: Res_{z = 0} \defrac{(z + l_k)^{M + d}}{z^{m_k} \prod_{k \neq i} ((z + l_{i})^{m_i}}
 *
 * Implementation: we expand the polynomial (z + l_k)^{M + d} and we expand the
 *                 series of 1 / (z + l_{i})^{m_i} up to degree m_k and multiply everything together.
 *                 Then, we return the coefficient of the m_k term.
 *
 * Background Math:
 * (x + y) ^ r , for any (real or complex r) = \sum _{k = 0}^{\inf} (r, k) x^{r - k}y^{k},
 * where (r, k) = r (r-1)(r-2)...(r-k+1)/k! = (r)_k/k!.
 * See http://en.wikipedia.org/wiki/Binomial_series
 *
 * So, let r = -s, and using the fact that (-s, k) = -s(-s -1)(-s -2)..(-s-k + 1)/k!= (s +k -1, k)*(-1)^k
 * 1/(x + b)^s = \sum _ {k = 0}^{\inf} (s +k -1, k)*(-1)^k * x^k * b^{-s - j}
 *             = \dfrac{1}{b^{m + s}} ( \sum _ {k = 0}^{\inf} (s +k -1, k)*(-1)^k * x^k * b^{m - j} )
 * 				such that m >= j. This is done because we cannot work with rational-coeff. series....yet.
 */
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
	RationalNTL c; //coefficient.
	int e[1]; //this is an array of size one because this is the exponent "vector" using the BurstTrie.
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

	//Brandon's notes: index[k] keeps the unique terms <l, s_i - s_j> for some fixed i.
	//				 : and counter[k] keeps track of how many times index[k] appears in the denominator.
	//				 : counter[0] = number of terms <l, s_i - s_j> that equals zero.

	//actual calculations, I want the appropriate coefficient in a product of one polynomial and some power series. (everything is truncated).
	nu = 1;
	de = 1;
	//for (i=1;i<=counter[0];i++) nu*=i;
	for (j = 1; j <= k - 1; j++)
		de = de * Power_ZZ(index[j], counter[j] + counter[0]);
	//cout << "de = ";
	//for (j = 1; j <= k - 1; j++)
	//	cout << index[j] << "^(" <<  counter[j] << '+' << counter[0] << ")  ";
	//cout << endl;

	monomialSum m1;
	monomialSum sub; //sub is the substitution for m1, which alternatively stores the product for each other
	m1.varCount = 1;
	m1.termCount = 0;
	sub.varCount = 1;
	sub.termCount = 0;
	//cout << "(x + " << p << ")^(" << M << '+' << d << ") = ";
	for (i = 0; i <= counter[0]; i++)
	{
		c = AChooseB(M + d, i) * Power_ZZ(p, M + d - i);
		e[0] = i;
		insertMonomial(c, e, m1); //now, m1 = <l, s_i>^(M+d) but in expanded form.
		//cout << c << "x^" << i << "  ";
	}
	//cout << endl;
	for (i = 1; i < k; i++)
	{
		monomialSum m2;
		m2.varCount = 1;
		m2.termCount = 0;
		//cout << "1/" << Power_ZZ(index[i], counter[i] + counter[0]) << "*(x + " << index[i] << ")^-" << counter[i] << "= ";
		for (j = 0; j <= counter[0]; j++)
		{
			c = AChooseB(counter[i] + j - 1, j) * Power_ZZ(index[i], counter[0]
					- j); //counter[0] - j - (counter[j] + counter[0] (which is in de)) = -counter[j] - j. :)
			if (j % 2 == 1)
				c.mult(to_ZZ(-1), to_ZZ(1));
			e[0] = j;
			insertMonomial(c, e, m2); //now, m2 = ???
			//cout << c/RationalNTL(Power_ZZ(index[i], counter[i] + counter[0]), to_ZZ(1)) << "x^" << j << "  ";
		};
		//cout << endl;
		mindeg[0] = 0;
		maxdeg[0] = counter[0];
		BTrieIterator<RationalNTL, int>* it = new BTrieIterator<RationalNTL, int> (); //MEMORY LEAK !?!??!?!
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
	//cout << "final poly = ";
	while (storedTerm = it->nextTerm())
	{
		//cout << storedTerm->coef/RationalNTL(de, to_ZZ(1)) << "*x^" << storedTerm->exps[0] << "  ";
		if (storedTerm->exps[0] == counter[0])
		{
			findCoeff = storedTerm->coef;
			//cout << endl;
			break;
		};
	};

	//cout << " all times " << nu << '/' << de << endl;
	a = nu * findCoeff.getNumerator();
	b = de * findCoeff.getDenominator();
	g = GCD(a, b);
	if (g != 0)
	{
		a = a / g;
		b = b / g;
	};
	return;
}//computeResidue


/**
 *This function is called when a vertex is irregular. The function computes the residue at the irregular vertex
 * d = dimension
 * m = power of the linear form.
 * innerProDiff[j] =  <l, s_i - s_j>
 * p = <l, s_i>
 * answer: output parameters.
 *
 * Formula: Res_{z = 0} \defrac{(z + l_k)^{M + d}}{z^{m_k} \prod_{k \neq i} ((z + l_{i})^{m_i}}
 *
 * Implementation: we expand the polynomial (z + l_k)^{M + d} and we expand the
 *                 series of 1 / (z + l_{i})^{m_i} up to degree m_k and multiply everything together.
 *                 Then, we return the coefficient of the m_k term.
 *
 * Background Math:
 * (x + y) ^ r , for any (real or complex r) = \sum _{k = 0}^{\inf} (r, k) x^{r - k}y^{k},
 * where (r, k) = r (r-1)(r-2)...(r-k+1)/k! = (r)_k/k!.
 * See http://en.wikipedia.org/wiki/Binomial_series
 *
 * So, let r = -s, and using the fact that (-s, k) = -s(-s -1)(-s -2)..(-s-k + 1)/k!= (s +k -1, k)*(-1)^k
 * 1/(x + b)^s = \sum _ {k = 0}^{\inf} (s +k -1, k)*(-1)^k * x^k * b^{-s - j}
 * */
void computeResidue(int d, int M, const vec_RationalNTL &innerProDiff, const RationalNTL &p,
		RationalNTL & finalAnswer)
{

	/**
	 * TODO: fix the memory leaks....counter, it*'s,
	 * print each polynomial series expansion (mul by de if needed) and then compare.
	 */

	finalAnswer = 0;
	if (p == 0)
	{
		return;
	}; //vertex vanishes, return 0;
	int k, i, j;
	int counter[1000];//counter counts number of appearances of each index[i]
	vec_RationalNTL index;//collecting different terms in the innerProDiff passed in
	bool found;
	RationalNTL de, nu, g; //c was ZZ
	RationalNTL c; //coefficient.
	c.setCanonicalizeFraction(false);
	int e[1]; //this is an array of size one because this is the exponent "vector" using the BurstTrie.
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

	//Brandon's notes: index[k] keeps the unique terms <l, s_i - s_j> for some fixed i.
	//				 : and counter[k] keeps track of how many times index[k] appears in the denominator.
	//				 : counter[0] = number of terms <l, s_i - s_j> that equals zero.

	//actual calculations, I want the appropriate coefficient in a product of one polynomial and some power series. (everything is truncated).
	nu = 1;
	de = 1;
	//for (i=1;i<=counter[0];i++) nu*=i;
	//for (j = 1; j <= k - 1; j++)
	///	de.mult(RationalNTL::power(index[j], counter[j] + counter[0]));
	//cout << "de = ";
	//for (j = 1; j <= k - 1; j++)
	//	cout << index[j] << "^(" <<  counter[j] << '+' << counter[0] << ")  ";
	//cout << endl;

	monomialSum m1;
	monomialSum sub; //sub is the substitution for m1, which alternatively stores the product for each other
	m1.varCount = 1;
	m1.termCount = 0;
	sub.varCount = 1;
	sub.termCount = 0;
	//cout << "(x + " << p << ")^(" << M << '+' << d << ") = ";
	for (i = 0; i <= counter[0]; i++)
	{
		//c = AChooseB(M + d, i) * Power_ZZ(p, M + d - i);
		c = p;
		c.power(M + d - i);
		c.mult(AChooseB(M + d, i));
		e[0] = i;
		insertMonomial(c, e, m1); //now, m1 = <l, s_i>^(M+d) but in expanded form.
		//cout << c << "x^" << i << "  ";
	}
	//cout << endl;
	for (i = 1; i < k; i++)
	{
		monomialSum m2;
		m2.varCount = 1;
		m2.termCount = 0;
		//cout << "1/" << Power_ZZ(index[i], counter[i] + counter[0]) << "*(x + " << index[i] << ")^-" << counter[i] << "= ";
		for (j = 0; j <= counter[0]; j++)
		{
			//c = AChooseB(counter[i] + j - 1, j) * Power_ZZ(index[i], counter[0]
			//		- j); //counter[0] - j - (counter[j] + counter[0] (which is in de)) = -counter[j] - j. :)

			c = index[i];
			c.power(((-1)*counter[i]) - j);
			c.mult(AChooseB(counter[i] + j - 1, j));
			if (j % 2 == 1)
				c.mult(to_ZZ(-1), to_ZZ(1));
			e[0] = j;
			insertMonomial(c, e, m2); //now, m2 = ???
			//cout << c/RationalNTL(Power_ZZ(index[i], counter[i] + counter[0]), to_ZZ(1)) << "x^" << j << "  ";
		};
		//cout << endl;
		mindeg[0] = 0;
		maxdeg[0] = counter[0];
		BTrieIterator<RationalNTL, int>* it = new BTrieIterator<RationalNTL, int> (); //MEMORY LEAK !?!??!?!
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
	//cout << "final poly = ";
	while (storedTerm = it->nextTerm())
	{
		//cout << storedTerm->coef/RationalNTL(de, to_ZZ(1)) << "*x^" << storedTerm->exps[0] << "  ";
		if (storedTerm->exps[0] == counter[0])
		{
			findCoeff = storedTerm->coef;
			//cout << endl;
			break;
		};
	};
	finalAnswer = findCoeff;
	finalAnswer.canonicalize();
	//finalAnswer = nu/de;
	//finalAnswer.mult(findCoeff);

	//cout << " all times " << nu << '/' << de << endl;
	//a = nu * findCoeff.getNumerator();
	//b = de * findCoeff.getDenominator();
	//g = GCD(a, b);
	//if (g != 0)
	//{
	//	a = a / g;
	//	b = b / g;
	//};
	return;



	/*
	finalAnswer = RationalNTL(); // = 0
	if (p == 0)
	{
		return;
	}; //vertex vanishes, return 0;

	int k, i, j;
	int * counter = new int[d];//counter counts number of appearances of each index[i]
	vec_RationalNTL index;//collects unique terms in the innerProDiff passed in
	bool found;
	ZZ de, nu;
	RationalNTL c; //coefficient.
	int e[1]; //this is an array of size one because this is the exponent "vector" using the BurstTrie.
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
	//count[0] = number of additional terms that subtract to zero, so count[0] +1 = number of zero subtractions.

	//int zeroCont = 0;
	//for(int cc = 0; cc <=d; ++cc)
	//	if ( innerProDiff[cc] == 0)
	//		++zeroCont;
	//cout << "zeroCount= " << zeroCont << endl;
	//cout << "counter[0]=" << counter[0] << endl;
    //So far:
	//index[k] keeps the unique terms <l, s_i - s_j> for some fixed i, that is passed in as innerProDiff
	//and counter[k] keeps track of how many times index[k] appears in the denominator.
	//counter[0] = number of terms <l, s_i - s_j> that equals zero.

	//actual calculations, I want the appropriate coefficient in a product of one polynomial and some power series. (everything is truncated).
	//nu = 1;
	//de = 1;
	//for (i=1;i<=counter[0];i++) nu*=i;
	//for (j = 1; j <= k - 1; j++)
	//	de = de * Power_ZZ(index[j], counter[j] + counter[0]);

	//alternatively save the product of the truncated polynomial series of 1 / (z + l_{i})^{m_i} and the current "answer polynomial" is m1 or sub.
	monomialSum m1;
	monomialSum sub; //sub is the substitution for m1, which alternatively stores the product.
	m1.varCount = 1;
	m1.termCount = 0;
	sub.varCount = 1;
	sub.termCount = 0;
	//cout << "(x + " << p << ")^(" << M << '+' << d << ") = ";
	for (i = 0; i <= counter[0]; i++)
	{
		//c = AChooseB(M + d, i) * Power_ZZ(p, M + d - i);
		c = p;
		c.power(M + d - i).mult(AChooseB(M + d, i));

		e[0] = i;
		insertMonomial(c, e, m1); //now, m1 = <l, s_i>^(M+d) but in expanded form.
		//cout << c << "x^" << i << "  ";
	}
	//cout << endl;
	for (i = 1; i < k; i++)
	{
		monomialSum m2;
		m2.varCount = 1;
		m2.termCount = 0;
		cout << "(x + " << index[i] << ")^-" << counter[i] << "= ";
		for (j = 0; j <= counter[0]; j++)
		{

			//c = AChooseB(counter[i] + j - 1, j) * Power_ZZ(index[i], counter[0]
			//		- j); //counter[0] - j - (counter[j] + counter[0] (which is in de)) = -counter[j] - j. :)

			c = index[i];
			c.power(-1 * counter[i] - j).mult(AChooseB(counter[i] + j - 1, j));
			//cout << index[i] << "^" << (-1 * counter[0] - j) << "=> " << RationalNTL(index[i]).power(-1 * counter[0] - j)
			//	 << " times " << AChooseB(counter[i] + j - 1, j) << " = " << c << endl;


			if (j % 2 == 1)
				c.mult(to_ZZ(-1));
			e[0] = j;
			insertMonomial(c, e, m2); //now, m2 = ???
			cout << c << "x^" << j << "  ";
		};
		cout << endl;
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
	//cout << "final poly = ";
	while (storedTerm = it->nextTerm())
	{
		//cout << storedTerm->coef << "*x^" << storedTerm->exps[0] << "  ";
		if (storedTerm->exps[0] == counter[0])
		{
			finalAnswer = storedTerm->coef;
			//cout << endl;
			return;
		};
	};

	cout << "computeResidue()(Rational)::Exponent " << counter[0] << "was not found" << endl;
	exit(1);
	//a = nu * findCoeff.getNumerator();
	//b = de * findCoeff.getDenominator();
	//g = GCD(a, b);
	//if (g != 0)
	//{
	//	a = a / g;
	//	b = b / g;
	//};
	//delete [] counter;
	//return;

	 */
}//computeResidue (RationalNTL)


