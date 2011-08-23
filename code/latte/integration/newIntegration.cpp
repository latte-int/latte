#include "newIntegration.h"
#include "iterators.h"
#include <NTL/vec_ZZ.h>
#include <NTL/ZZ.h>
#include <NTL/mat_ZZ.h>
#include <iostream>
#include "print.h"

using namespace std;

NTL_CLIENT

//this function deletes space from a given string
void delSpace(string &line)
{
	for (int i = 0; i < line.length(); i++)
	{
		while ((i < line.length()) && (line.at(i) == 32))
		{
			line.erase(i, 1);
		}
	}
}//delSpace


//this function converts a given string into a simlexZZ mySimplex
//for example, string [[0,0],[1,1],[7,8]] is converted to a two-dimensional vector of ZZs.
void convertToSimplex(simplexZZ &mySimplex, string line)
{
	delSpace(line);
	int index, i, t, j, c;
	string temp, subtemp;
	t = 2;
	mySimplex.d = 1;
	t = line.find("[", t) + 1;
	t = line.find("[", t) + 1;
	temp = line.substr(t, line.find("]", t) - t);
	for (i = 0; i < temp.length(); i++)
		mySimplex.d += (temp.at(i) == ',');
	c = 0;
	for (i = 0; i < line.length(); i++)
		c += (line.at(i) == ']');
	if (c - 2 != mySimplex.d)
	{
		cout << "The d-simplex should have d+1 vertices. Please check." << endl;
		exit(1);
	};
	(mySimplex.s).SetLength(mySimplex.d + 1);
	index = 1;
	for (i = 0; i <= mySimplex.d; i++)
	{
		temp = line.substr(index, line.find("]", index) - index + 1);
		c = 0;
		for (j = 0; j < temp.length(); j++)
			c += (temp.at(j) == ',');
		if (c != mySimplex.d - 1)
		{
			cout << "Each vertex should have d coordinates. Please check."
					<< endl;
			exit(1);
		};
		(mySimplex.s[i]).SetLength(mySimplex.d);
		t = 1;
		for (j = 0; j < mySimplex.d - 1; j++)
		{
			subtemp = temp.substr(t, temp.find(",", t) - t);
			t = temp.find(",", t) + 1;
			mySimplex.s[i][j] = to_ZZ(subtemp.c_str());
		};
		subtemp = temp.substr(t, temp.find(",", t) - t + 1);
		t = temp.find(",", t);
		mySimplex.s[i][mySimplex.d - 1] = to_ZZ(subtemp.c_str());
		index = line.find("]", index) + 2;
	};
	mat_ZZ matt;
	matt.SetDims(mySimplex.d, mySimplex.d);
	for (i = 1; i <= mySimplex.d; i++)
		matt[i - 1] = mySimplex.s[i] - mySimplex.s[0];
	mySimplex.v = determinant(matt);
	if (mySimplex.v < 0)
		mySimplex.v = -mySimplex.v;
}
;

/**
 * Integrate a simplex over a linear form.
 *
 * @parm a, b: ouput parameters, we return a/b += integration answer.
 * @parm l: a linear form.
 * @parm mySimplex: integer simplex
 * @parm m: the power the linear form is raised to
 * @parm coe: the coefficient of a linear form
 * @parm de: is the extra factor in the formulae that we we multiply the result by
 *
 * ASSUMES the polytope has dimension less than 1000.
 *
 * Paper Citation: @ARTICLE
 * {2008arXiv0809.2083B,
 *  author = {{Baldoni}, V. and {Berline}, N. and {De Loera}, J. and {K{\"o}ppe}, M. and
 *	{Vergne}, M.},
 *    title = "{How to Integrate a Polynomial over a Simplex}",
 *  journal = {ArXiv e-prints},
 *archivePrefix = "arXiv",
 *   eprint = {0809.2083},
 *primaryClass = "math.MG",
 * keywords = {Mathematics - Metric Geometry, Computer Science - Computational Complexity, Computer Science - Symbolic Computation},
 *     year = 2008,
 *    month = sep,
 *   adsurl = {http://adsabs.harvard.edu/abs/2008arXiv0809.2083B},
 *  adsnote = {Provided by the SAO/NASA Astrophysics Data System}
 * }
 *
 *BACKGROUND MATH 1: \int_\Delta l^m \d m' = d!\vol(\Delta, \d m')\frac{m!}{(m+d)!}
 *                            \Big(\sum_{i=1}^{d+1}
 *                                  \frac{ <l, s_i >^{M+d}}
 *                                       {\prod_{j\neq i} <l, s_i- s_j >}
 *                            \Big),
 * Where \Delta is a regular integer-vertex simplex, d is the dimension, l is a linear form, s_i are the verties.
 *
 *BACKGROUND MATH 2:
 *  \int_{\Delta} l^m  \d m' = d!\vol(\Delta, \d m') \frac{m!}{(m+d)!}
 *                             \sum_{k\in K} \Res_{z=0}
 *                                  \frac{(z + <l, s_k>)^{m+d}}
 *                                       {\z^{m_k} {\prod_{i\in K, i \neq k} {(z + <l, s_k - s_i> )}^{m_i}} }
 * as above and if the simplex is not regular, and where $K\subseteq\{1,\dots,d+1\}$ is an index set of the different poles
 * $t= 1/\langle \l ,s_k\rangle$, and for $k\in K$ let $m_k$ denote the order of the pole, i.e.,
 *   $m_k = sizeof the set { i \in \{1,\dots,d+1\} : <l ,s_i> = <l ,s_k> }.
 *
 * Implementation:
 * Instead of finding d!\vol(\Delta, \d m') directly, we just find the volume of the parallelepiped of the simplex.
 * The data structure also assumes the m! is part of the linear form's coefficient.
 */
void update(ZZ &a, ZZ &b, vec_ZZ l, simplexZZ mySimplex, int m, RationalNTL coe, ZZ de)
{

	ZZ sum, lcm, total, g, tem;
	int i, j;
	vec_ZZ inner_Pro; //inner_Pro[i] = <l, s_i>
	vec_ZZ sum_Nu, sum_De; // (sum_Nu/sum_De)[i] = <l, s_i>^d/ (\prod_{j \neq i} <l, s_i - s_j>)
	inner_Pro.SetLength(mySimplex.d + 1);
	sum_Nu.SetLength(mySimplex.d + 1);
	sum_De.SetLength(mySimplex.d + 1);
	total = 0;
	lcm = 1;
	bool repeat[1000]; //put this on the stack, do not waste the time requesting memory from the heap because this function is called many, many, many times.
					 //Why is this bool (vs int): if there are no repeats in the <l, s_i> terms, the simplex is regular on l and we compute the integral as in the first case of the theory.
					 // Otherwise we will have to compute the residue. It is in the residue-function where we worry about the multiplicity of things.
	for (i = 0; i <= mySimplex.d; i++)
	{
		sum = 0;
		for (j = 0; j < mySimplex.d; j++)
			sum = sum + l[j] * mySimplex.s[i][j];

		inner_Pro[i] = sum; // inner_Pro_i= <l, s_i>
		repeat[i] = 0;
		for (j = 0; j < i; j++)
			if (inner_Pro[j] == inner_Pro[i])
			{
				repeat[i] = 1;
				break;
			};//record repetitions
	};//stores inner product for use

	for (i = 0; i <= mySimplex.d; i++)
		if (!repeat[i])
		{
			sum_Nu[i] = 1;
			//cout << "update: l= " << l;
			//cout << "update: v=";

			//for(int index = 0; index <= mySimplex.d; ++index) cout << mySimplex.s[index] << ", ";
			//cout << endl;
			//cout << "update::l^dim+m=" << inner_Pro[i] << "^" << mySimplex.d << "+" << m;
			for (j = 0; j < m + mySimplex.d; j++)
				sum_Nu[i] = sum_Nu[i] * inner_Pro[i]; // sum_Nu_i = inner_pro ^ (m + d)
			//cout << "=" << sum_Nu[i] << endl;
			sum_De[i] = 1;
			for (j = 0; j <= mySimplex.d; j++)
				if (i != j)
					sum_De[i] = sum_De[i] * (inner_Pro[i] - inner_Pro[j]);  //sum_de_i	 = \prod _{i != j} <l, s_i - s_j>
			if ((sum_Nu[i] < 0) && (sum_De[i] < 0))
			{
				sum_Nu[i] = -sum_Nu[i];
				sum_De[i] = -sum_De[i];
			};
			if (sum_De[i] == 0)
			{
				vec_ZZ ProDiff;
				ProDiff.SetLength(mySimplex.d + 1);
				for (j = 0; j <= mySimplex.d; j++)
					ProDiff[j] = inner_Pro[i] - inner_Pro[j];

				computeResidue(mySimplex.d, m, ProDiff, inner_Pro[i],
						sum_Nu[i], sum_De[i]);
			}
			if (sum_De[i] != 0)
			{
				lcm = lcm * sum_De[i] / (GCD(lcm, sum_De[i]));
			};
			//cout << "update:i= " << i << "num/dem= " << RationalNTL(sum_Nu[i], sum_De[i]) << endl;
			//cout << "update:i= " << i << "num/dem= " << sum_Nu[i]<< " / "<< sum_De[i] << endl;
			//cout << "update:i= " << i << "num/dem= " << sum_Nu[i]/GCD(sum_Nu[i],sum_De[i])<< " / "<< sum_De[i]/GCD(sum_Nu[i],sum_De[i]) << endl;
			//cout << RationalNTL(coe.getNumerator()*mySimplex.v*sum_Nu[i],de*coe.getDenominator()*sum_De[i]) << endl;

		};

	for (i = 0; i <= mySimplex.d; i++)
		if ((!repeat[i]) && (sum_De[i] != 0))
		{
			total += sum_Nu[i] * (lcm / sum_De[i]);
		}


	lcm = lcm * de * coe.getDenominator();
	total = total * mySimplex.v * coe.getNumerator();

	//cout << "update: total/lcm = " << RationalNTL(total,lcm) << endl;
	if (a == 0)
	{
		a = total;
		b = lcm;
	} //return a/b = total/lcm * coe
	else if ((lcm != 0) && (b != 0))
	{
		// a/b = a/b + total/lcm
		tem = b * lcm / GCD(b, lcm); //find LCM of b and lcm.
		a = a * tem / b + total * tem / lcm;
		b = tem;

	} //return a/b := a/b + total/lcm.
	g = GCD(a, b);
	if (g != 0)
	{
		a = a / g;
		b = b / g;
	}
}//update



//This function computes a given fraction a/b, the integral of the linear form forms, over the simplex mySimplex
void integrateLinFormSum(ZZ& numerator, ZZ& denominator,
		PolyIterator<RationalNTL, ZZ>* it, const simplexZZ &mySimplex)
{
	ZZ v, de, counter, tem; //, coe;
	RationalNTL coe;
	int i, j, index, k, m;
	vec_ZZ l;
	//if (forms.varCount!=mySimplex.d) {cout<<"The dimensions of the polynomial and simplex don't match. Please check!"<<forms.varCount<<"<>"<<mySimplex.d<<endl;exit(1);};
	l.SetLength(mySimplex.d);
	numerator = 0;
	denominator = 0;
	it->begin();
	term<RationalNTL, ZZ>* temp;
	while (temp = it->nextTerm())
	{
		coe = temp->coef;
		m = temp->degree; //obtain coefficient, power
		l.SetLength(temp->length); //obtain exponent vector
		for (j = 0; j < temp->length; j++)
		{
			l[j] = temp->exps[j];
		}
		de = 1;
		for (i = 1; i <= mySimplex.d + m; i++)
		{
			de = de * i;
		} //de is (d+m)!. Note this is different from the factor in the paper because in our storage of a linear form, any coefficient is automatically adjusted by m!
		update(numerator, denominator, l, mySimplex, m, coe, de);//We are ready to compute the integral of one linear form over the simplex
	}
	delete temp;
	if (denominator < 0)
	{
		denominator *= to_ZZ(-1);
		numerator *= to_ZZ(-1);
	}
}//integrateLinFormSum

/* computes the integral of a product of powers of linear forms over a  simplex
 * @input it: iterator to the linear forms.
 * @input mySimplex: a simplex
 * @input productCount: the number of products in the linear form.
 * @return RationalNTL: the integral over the simplex.
 *
 * Math: integral(over simplex) <l_1, x>^m_1 ... <l_d, x>^m_D =
 *
 *               abs(det(matrix formed by the rays))* M!/(d+|M|)!
 *  --------------------------------------------------------------------------
 *  product(over j)(1 - <l_1, s_j>t_1 - <l_2, s_j>t_2 - ... - <l_D, s_j>t_D )
 *                                 \          \_these are numbers (tVector)
 *                                  \_the t are symbolic.
 *
 *  where we want the coefficient of t_1^m_1 ... t_D^m_D in the polynomial
 *  expansion of the RHS,
 *
 * and where s_j is a vertex,
 *       l_i is a linear form
 *       D is any number of products
 *       M is the power vector
 *       M! = m_1! m_2! ... m_D!
 *       |M| = m_1 + m_2 + ... + m_D.
 *
 * Note that we cannot divide by zero, so don't worry.
 *
 * To find such an expansion, we use a Tayler expansion on each
 * 1/(1 +a_1t_1 + ... +a_Dt_D ) up to degree M for each such product. We then
 * multiply all these polynomials together (ignoring degrees larger than M)
 * and then we find the coefficient of t^M. (Again, M is a vector of powers)
 *
 * Taylor Expansion about zero (x is a vector)
 * f(x) = 1/(1-a_1x_1 + ... + a_d x_d)
 *      = sum_{n \in \Z^d_{>= 0}}  x^n / n! * d^{n|} f(x)     (0,0,...0)
 *                                         ------------------
 *                                        dx_1^{n_1}...dx_d^{n_d}
 *      = sum_{n \in \Z^d_{>= 0}}  x^n / n! * (-1)(-2)...(-1 * |n|) (1 - 0)^{-|n|-1} a_1^{n_1} ... a_d^{n_d}
 *      = sum_{n \in \Z^d_{>= 0}}  x^n (-1)^|n|  (|n| choose n_1, ..., n_d) a_1^{n_1} ... a_d^{n_d}
 *      where (|n| choose n_1, ..., n_d) is a multinomial coefficient.
 *
 * See the paper: "How to Integrate a Polynomial over a Simplex" by  V. BALDONI, N. BERLINE, J. A. DE LOERA, M. VERGNE.
 */
RationalNTL integrateLinFormProducts(PolyIterator<RationalNTL, ZZ>* it, const simplexZZ &mySimplex, const int productCount)
{
	//cout << "integrateLinFormProducts called" << endl;
	int * M; //the power vector.
	ZZ lenM; // |M|
	RationalNTL coef; //temp. coefficient.
	RationalNTL answer;
	ZZ monomialCount; //number of power vectors less than M component wise.

	monomialCount = 1;
	coef = 1;

	M = new int[productCount];
	it->begin();
	term<RationalNTL, ZZ>* temp;
	int i = 0;
	while (temp = it->nextTerm())
	{
		M[i] = temp->degree; //save the power
		++i;
		lenM += temp->degree; //add the power

		coef *= temp->coef; // M1! M2! ... MD! * (coefficents ^ powers)

		monomialCount *= (temp->degree+1); //monomialCount = number of monomials (m1, ..., md) <= (deg1, ..., degD).
	}
	if ( i != productCount)
		THROW_LATTE_MSG(LattException::ue_BadPolynomialLinFormInput, "count of terms differ");

	ZZ factorialDim = to_ZZ(1); // = (|M|+d)!
	for (ZZ j = to_ZZ(2); j <=  lenM + mySimplex.d; ++j)
	{
		factorialDim *= j;
	}

	// 1/factorialDim = 1 / (|M| + d)!

	answer = coef * mySimplex.v;
	answer.div(factorialDim);

	//now, answer = [vol(simplex)*d!] M! / (|M|+d)!.


	//ok, now we just need to find the coeff of M in the polynomial expansion.


	vec_ZZ tVector; //the coefficent vector of ( 1- a_1t_1 - ... - a_Dt_D) (we don't save the leading 1)
	int* counter; //current power n
	tVector.SetLength(productCount);
	counter = new int[productCount];
	int* minDegree = new int[productCount];

	monomialSum polynomialProduct;
	polynomialProduct.termCount = 0;
	polynomialProduct.varCount = productCount;

	for(int i = 0; i < productCount; ++i)
		minDegree[i] = 0;

	//insert 1.
	insertMonomial(RationalNTL(1,1), minDegree, polynomialProduct);



	//for every 1/(1 + stuff) term in the denominator.
	for(int i = 0; i < mySimplex.d +1; ++i)
	{
		it->begin();
		int j = 0;
		while (temp = it->nextTerm())
		{
			tVector[j] = 0;
			for (int len = 0; len < temp->length; ++len)
			{
				tVector[j] += mySimplex.s[i][len] * temp->exps[len];
			}//dot l  and vertex i
			tVector[j] *= -1;

			++j;
		}//while. build the t-vector.

		//now tVector is in the form [-<l_1, s_i>, ..., -<l_D, s_i>]
		//expand tVecotr into a polynomial series.


		monomialSum thePolynomial;
		thePolynomial.termCount = 0;
		thePolynomial.varCount = productCount;

		for(j = 0; j < productCount; ++j)
			counter[j] = 0;

		insertMonomial(RationalNTL(1,1), counter, thePolynomial);

		//for every monomial less than or equal to M.
		for(ZZ currentMonomialCount = to_ZZ(0); currentMonomialCount < monomialCount-1 /*-1 because we already added 1x^0*/; ++currentMonomialCount)
		{
			counter[0] += 1;
			for (int myIndex = 0; counter[myIndex] > M[myIndex]; myIndex++)
			{
				counter[myIndex] = 0;
				counter[myIndex + 1] += 1;
			}//add one and do a carry if you have to.

			//insert the polynomial in counter.

			RationalNTL c(1,1);
			int lenC = 0; // = |n|
			for(int k = 0; k < productCount; ++k)
				lenC += counter[k];

			int lenC_copy = lenC;

			for(int k = 0; k < productCount; ++k)
			{
				c *= AChooseB(lenC, counter[k]);
				c *= power(tVector[k], counter[k]);
				lenC -= counter[k];
			}// find (|n| choose n_1, ..., n_d) and a^n

			if ( lenC_copy % 2 == 1)
				c.changeSign(); // (-1)^|n|

			//cout << "going to insert polynomail term" << endl;
			//cout << " c=" << c << " exp= ";
			//for(int k = 0; k < productCount; ++k)
			//	cout << counter[k] << ", " ;
			//cout << endl;
			insertMonomial(c, counter, thePolynomial);
			//cout << "end going to insert polynomail term" << endl;

		}//make the polynomial.

		//cout << "\n\ndoing poly multiplication" << endl;
		BTrieIterator<RationalNTL, int>* it3 = new BTrieIterator<RationalNTL, int> ();
		BTrieIterator<RationalNTL, int>* it2 = new BTrieIterator<RationalNTL, int> ();
		it3->setTrie(thePolynomial.myMonomials, thePolynomial.varCount);
		it2->setTrie(polynomialProduct.myMonomials, polynomialProduct.varCount);

		monomialSum tempProducts;
		tempProducts.varCount = productCount;
		multiply<RationalNTL> (it3, it2, tempProducts, minDegree, M);
		destroyMonomials(thePolynomial);
		destroyMonomials(polynomialProduct);
		delete it3;
		delete it2;

		polynomialProduct.myMonomials = tempProducts.myMonomials;
		polynomialProduct.termCount = tempProducts.termCount;
		polynomialProduct.varCount =  tempProducts.varCount;
		//cout << "end poly multiplication" << endl;
	}//for every simplex vertex

	//now find the coeff of M.

	BTrieIterator<RationalNTL, int>* finalPoly = new BTrieIterator<RationalNTL, int> ();
	finalPoly->setTrie(polynomialProduct.myMonomials, polynomialProduct.varCount);
	finalPoly->begin();
	term<RationalNTL, int>* storedTerm;
	bool found;
	while (storedTerm = finalPoly->nextTerm())
	{
		found = true;
		for(int i = 0; i < productCount; ++i)
			if ( storedTerm->exps[i] != M[i])
			{
				found = false;
				break;
			}

		if (found == true )
		{
			answer.mult(storedTerm->coef);
			break;
		}//found coeff. of highest term.
	}//while
	assert(found == true);

	return answer;
}//integrateLinFormProducts

void integrateMonomialSum(ZZ &a, ZZ &b, monomialSum &monomials,
		const simplexZZ &mySimplex)//integrate a polynomial stored as a Burst Trie
{
	linFormSum forms;
	forms.termCount = 0;
	forms.varCount = monomials.varCount;
	BTrieIterator<RationalNTL, int>* it = new BTrieIterator<RationalNTL, int> ();
	it->setTrie(monomials.myMonomials, monomials.varCount);
	decompose(it, forms); //decomposition
	delete it;
	BTrieIterator<RationalNTL, ZZ>* it2 = new BTrieIterator<RationalNTL, ZZ> ();
	it2->setTrie(forms.myForms, forms.varCount);
	integrateLinFormSum(a, b, it2, mySimplex);
}


void _integrateMonomialSum(ZZ &a, ZZ &b, _monomialSum &monomials,
		const simplexZZ &mySimplex)
{
	_linFormSum forms;
	forms.termCount = 0;
	forms.varCount = monomials.varCount;
	for (int i = 0; i < monomials.termCount; i++)
		_decompose(monomials, forms, i);
	LBlockIterator<RationalNTL>* it_ = new LBlockIterator<RationalNTL> ();
	it_->setLists(forms.lHead, forms.cHead, forms.varCount, forms.termCount);
	integrateLinFormSum(a, b, it_, mySimplex);
}

