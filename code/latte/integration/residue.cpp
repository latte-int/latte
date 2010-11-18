#include <NTL/vec_ZZ.h>
#include <NTL/ZZ.h>
#include "PolyTrie.h"
#include "multiply.h"
#include "valuation/Perturbation.h"

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
 * @parm d: dimension
 * @parm m: power of the linear form.
 * @parm innerProDiff[j] =  <l, s_i - s_j>
 * @parm p = <l, s_i>
 * @parm a, b: output parameters: answer: a/b.
 *
 * Formula: Res_{z = 0} \defrac{(z + l_k)^{M + d}}{z^{m_k} \prod_{k \neq i} ((z + l_{i})^{m_i}}
 *
 * Implementation Goal: we expand the polynomial (z + l_k)^{M + d} and we expand the
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
 * 				such that m >= j. This is done because we cannot work with rational-coeff. series...we
 * 				don't want b^(negative number), so we increase the powers of b in the series and divide by
 * 				b^{m + s} at the end.
 */
void computeResidue(int d, int M, const vec_ZZ &innerProDiff, const ZZ &p,
		ZZ &a, ZZ &b)
{
	//cout << "Compute residue called" << endl;
	//cout << "d=" << d << "M=" << M << " innerprod=" << innerProDiff << "p= " << p << endl;
	if (p == 0)
	{
		a = 0;
		b = 1;
		return;
	}; //vertex vanishes, return 0;
	int k, i, j;
	int counter[1000];//counter counts number of appearances of each index[i]
					//again, put this on the stack. Don't want the time requesting memory from the heap because this function is called many times.
	vec_ZZ index;//collecting different terms in the innerProDiff passed in
	bool found;
	ZZ de, nu, g;
	RationalNTL c; //coefficient
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
	//so far we've been doing book keeping stuff: index stores the UNIQUE differences and counter stores the multiplicity


	//Brandon's notes: index[k] keeps the unique terms <l, s_i - s_j> for some fixed i.
	//				 : and counter[k] keeps track of how many times index[k] appears in the denominator.
	//				 : counter[0] = number of terms <l, s_i - s_j> that equals zero, which has to at least one, otherwise computeResidue would not have been called.
	//				 : 		Thus, we take counter[0]--, and so counter[0] now is equal to the number of additional terms that vanish.
	//				 : 		So now, counter[0] is our upper bound for how far we need to take the series expansion.
	//				 :		To find the residue, we need to find the coeff. of the counter[0]-degree term because we are dividing by z^{counter[0]+1}

	//actual calculations, I want the appropriate coefficient in a product of one polynomial and some power series. (everything is truncated).
	nu = 1;
	de = 1;
	//for (i=1;i<=counter[0];i++) nu*=i;
	for (j = 1; j <= k - 1; j++)
		de = de * Power_ZZ(index[j], counter[j] + counter[0]);
	monomialSum m1;
	monomialSum sub; //sub is the substitution for m1, which alternatively stores the product for each other
					 // we alternatively do sum:= m1 * m2  and then m1:=sub * m2.
	m1.varCount = 1;
	m1.termCount = 0;
	sub.varCount = 1;
	sub.termCount = 0;
	sub.myMonomials = NULL;
	for (i = 0; i <= counter[0]; i++)
	{
		c = AChooseB(M + d, i) * Power_ZZ(p, M + d - i);
		e[0] = i;
		insertMonomial(c, e, m1);
	}//now, m1 = <l, s_i>^(M+d) but in expanded form.
	for (i = 1; i < k; i++)
	{
		monomialSum m2;
		m2.varCount = 1;
		m2.termCount = 0;
		for (j = 0; j <= counter[0]; j++)
		{
			c = AChooseB(counter[i] + j - 1, j) * Power_ZZ(index[i], counter[0]
					- j);
			//index[i]^{counter[0] - j - (counter[j] + counter[0] (which is in de))} = -counter[j] - j. :)

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
			if ( sub.myMonomials != NULL)
				destroyMonomials(sub);
			sub.varCount = 1;
			multiply<RationalNTL> (it, it2, sub, mindeg, maxdeg);
		}//cout<<"times "<<printMonomials(m2)<<" gives "<<printMonomials(sub)<<endl;}
		else
		{
			it->setTrie(sub.myMonomials, sub.varCount);
			it2->setTrie(m2.myMonomials, m2.varCount);
			destroyMonomials(m1);
			m1.varCount = 1;
			multiply<RationalNTL> (it, it2, m1, mindeg, maxdeg);
		}//cout<<"times "<<printMonomials(m2)<<"gives "<<printMonomials(m1)<<endl;};
		delete it;
		delete it2;
		destroyMonomials(m2);
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
	//BurstTerm<RationalNTL, int>* temp;
	BurstTrie<RationalNTL, int>* myTrie;
	BTrieIterator<RationalNTL, int>* it = new BTrieIterator<RationalNTL, int> ();
	if (k % 2 == 1)
	{
		//temp = new BurstTerm<RationalNTL, int> (m1.varCount);
		myTrie = m1.myMonomials;
		it->setTrie(myTrie, m1.varCount);
	} else
	{
		//temp = new BurstTerm<RationalNTL, int> (sub.varCount);
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
		}//again, counter[0] +1 = degree of z, which we are dividing by.
	}//while

	a = nu * findCoeff.getNumerator();
	b = de * findCoeff.getDenominator();
	g = GCD(a, b);
	if (g != 0)
	{
		a = a / g;
		b = b / g;
	};

	//cout << "compute residue: d=" << d << ", M=" << M << ", p=" << p << endl;
	//cout << "  innerProdDiff" << innerProDiff << endl;
	//cout << "  a/b = " << RationalNTL(a, b) << endl;

	delete it;
	destroyMonomials(m1);
	destroyMonomials(sub);
	return;
}//computeResidue


/**
 * This function is very similar to the computeResidue() function, but
 * instead of finding the series expansion of 1/(a+e)^m about e= 0,
 * we need to take the expansion 1/(a+b*e)^m about e=0.
 *
 * Assumes the numerator does not vanish.
 * Assumes the first entry of lDotR is zero and the first entry of leDotRPower is the order of the pole.
 * We find Residue about e=0 of( (lDotV +eDotV*e)^d+m / (e^leDotRPower[0]+1 * (lDotR[i] +  eDotR[i]*e)^leDotRPower[i]))
 *                                               \_>note the +1.
 * This function is a friend function to LinearLawrenceIntegration
 * @parm d: dimension
 * @parm m: power of the linear form.
 * @parm coneTerm: <v,l+e>/ e^m1*(a1 +b1e)^m2*...(ak +bke)^mk, k<d.(note that the power of <v,e+l> is not saved here.
 * @parm numerator, denominator:output parameters
 *
 * Formula: Res_{z = 0} \defrac{(z + l_k)^{M + d}}{z^{m_k} \prod_{k \neq i} ((z + l_{i})^{m_i}}
 *
 * Implementation: Use regular binomial theorem to expand <v,l+e>^(m+d) and delete terms of larger than degree m1.
 *                 Expand 1/(a+be)^mi using the negative binomial theorem and delete powers larger than degree m1.
 *                 Multiply these finite powers, keeping the degree <= m1 terms.
 *                 Then find the constant coeff. term in this final polynmial product.
 * For the negative binomial theorem: http://en.wikipedia.org/wiki/Binomial_series
 *
 * Like in computeResideu, when expanding 1/(a + be)^m up to degree m1,
 * we factor a 1/(a^(m+m1)) so that the finite-degree polynomial will have
 * integer values.See the notes for computeResidue()
 */
void computeResidueLawrence(const int d, const int M, const LinearLawrenceIntegration & coneTerm, ZZ &numerator, ZZ &denominator)
{

	//cout << "computeResidueLawrence" << endl;
	//cout << "  d=" << d << ", M=" << M << " ";
	//coneTerm.printTerm();

	int k, i, j;

	//int counter[1000];//counter counts number of appearances of each index[i]
					//again, put this on the stack. Don't want the time requesting memory from the heap because this function is called many times.
	//vec_ZZ index;//collecting different terms in the innerProDiff passed in
	//bool found;

	ZZ de, nu, g;
	RationalNTL c; //coefficient.
	int e[1]; //this is an array of size one because this is the exponent "vector" using the BurstTrie.
	int mindeg[1];
	int maxdeg[1];

	int truncateDegree;
	truncateDegree = coneTerm.rayDotProducts[0].power; //want to find the coef. of the truncateDegree-degree term in the final polynomial.
			//we assume the power of the (0+c1*e) term is in the first array index.


	nu = 1;
	de = coneTerm.rayDotProducts[0].epsilon; //factor the coeff. of e^m out.

	monomialSum products;
	monomialSum tempProducts;
	products.varCount = 1;
	products.termCount = 0;
	tempProducts.varCount = 1;
	tempProducts.termCount = 0;
	tempProducts.myMonomials = NULL;
	//cout << "(" << coneTerm.numeratorDotProduct.constant << "+ " <<  coneTerm.numeratorDotProduct.epsilon << ") ^" << M + d << "==" << endl;
	for (i = 0; i <= truncateDegree; i++)
	{
	    //MATH: (a + be)^M+d = sum_k=0^m+d (m+d choose k) (be)^k * (a)^{m+d -k}
		//                                                  \         \_>numeratorDotProduct.constant
		//                                                   \_>coneTerm.numeratorDotProduct.epsilon
		c = AChooseB(M + d, i) * Power_ZZ(coneTerm.numeratorDotProduct.constant, M + d - i) * Power_ZZ(coneTerm.numeratorDotProduct.epsilon, i);
		e[0] = i;
		//cout << c << "e^" << e[0] << " + ";
		insertMonomial(c, e, products);
	}//now, m1 = <l +e, vertex>^(M+d) but in expanded form and truncated
	//cout << endl;
	//start i at 1 because first index (0) assumed to be order of pole....(0+ce)
	for (i = 1; i < d; i++)
	{
		//cout << "going to do case i=" << i << endl;
		if (coneTerm.rayDotProducts[i].power <= 0 )
			continue; //really, at this point, the power should not be zero. It could be negative if this term is a repeat or positive.
		if (coneTerm.rayDotProducts[i].epsilon == 0)
		{
			//cout << "factored anoter term: " << coneTerm.rayDotProducts[i].constant << "^" << coneTerm.rayDotProducts[i].power << endl;
			de *= Power_ZZ(coneTerm.rayDotProducts[i].constant, coneTerm.rayDotProducts[i].power);
			continue;
		}//factor the constant out: (a + 0*e)^power.
		//factor out a so that the polynomial is integer.
		//cout << "factoring out: " << coneTerm.rayDotProducts[j].constant << "^" << coneTerm.rayDotProducts[j].power + truncateDegree << endl;
		de = de * Power_ZZ(coneTerm.rayDotProducts[i].constant, coneTerm.rayDotProducts[i].power + truncateDegree);


		monomialSum m2;
		m2.varCount = 1;
		m2.termCount = 0;
		//cout << "series: ";
		for (j = 0; j <= truncateDegree; j++)
		{
			//MATH: (a + be)^{-s} truncated to degree m in terms of e:
			//      = b^{-s -m} sum_j=0^k=m (-1)^j * (s+j-1 choose j) * (be)^j * a^{m -k}
			//              \                                             \        \_>coneTerm.rayDotProducts[i].constant
			//               \                                             \_>coneTerm.rayDotProducts[i].epsilon
			//                \_>factored out into de agove.
			c = AChooseB(coneTerm.rayDotProducts[i].power + j - 1, j) * Power_ZZ(coneTerm.rayDotProducts[i].constant, truncateDegree	- j) * Power_ZZ(coneTerm.rayDotProducts[i].epsilon, j);

			if (j % 2 == 1)//this sign comes from (-s, k) = -s(-s -1)(-s -2)..(-s-k + 1)/k!= (s +k -1, k)*(-1)^k
				c.mult(to_ZZ(-1));
			e[0] = j;
			//cout << c << "e^" << j << " + ";
			insertMonomial(c, e, m2);
		};
		//cout << endl;

		//I took out the whole polynomial multiplication flip-flop code found in computeResidue().

		mindeg[0] = 0;
		maxdeg[0] = truncateDegree;
		BTrieIterator<RationalNTL, int>* it = new BTrieIterator<RationalNTL, int> ();
		BTrieIterator<RationalNTL, int>* it2 = new BTrieIterator<RationalNTL, int> ();
		it->setTrie(products.myMonomials, products.varCount);
		it2->setTrie(m2.myMonomials, m2.varCount);
		multiply<RationalNTL> (it, it2, tempProducts, mindeg, maxdeg);//tempProducts = m2 * products.
		destroyMonomials(products);
		products.myMonomials = tempProducts.myMonomials;
		products.termCount = tempProducts.termCount;
		products.varCount =  tempProducts.varCount;

		delete it;
		delete it2;
		destroyMonomials(m2);
	};
	//cout << "end of denominator processing" << endl;
	//ZZ findCoeff = to_ZZ(0);
	RationalNTL findCoeff;

	///find a monomial that has degree truncateDegree and get its coefficient
	BTrieIterator<RationalNTL, int>* it = new BTrieIterator<RationalNTL, int> ();
	//myTrie = products.myMonomials;
	it->setTrie(products.myMonomials, products.varCount);
	it->begin();
	term<RationalNTL, int>* storedTerm;
	while (storedTerm = it->nextTerm())
	{
		if (storedTerm->exps[0] == truncateDegree)
		{
			findCoeff = storedTerm->coef;
			break;
		}//found coeff. of highest term.
	}//while



	numerator = nu * findCoeff.getNumerator();
	denominator = de * findCoeff.getDenominator();


	//cout << "*compute residue: d=" << d << ", M=" << M << " " ;
	//coneTerm.printTerm();
	//cout << "*  num/den = " << RationalNTL(numerator, denominator) << endl;


	delete it;
	destroyMonomials(products);
	return;

}//computeResidueLawrence




