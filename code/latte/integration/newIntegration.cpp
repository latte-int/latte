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

			for(int index = 0; index <= mySimplex.d; ++index) cout << mySimplex.s[index] << ", ";
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

void updateLawrence(ZZ &a, ZZ &b, vec_ZZ l, listCone *cone, int m, RationalNTL coe, ZZ de, int dim)
{
	//cout << "m=" << m << " dim=" << dim << endl;
	ZZ sum, lcm, total, g, tem, det;
	int i, j;
	ZZ temp_a, temp_b;
	temp_a = 0;
	temp_b = 1;

	//vec_ZZ inner_Pro; //inner_Pro[i] = <l, s_i>
	//vec_ZZ sum_Nu, sum_De; // (sum_Nu/sum_De)[i] = <l, s_i>^d/ (\prod_{j \neq i} <l, s_i - s_j>)
	//inner_Pro.SetLength(mySimplex.d + 1);
	//sum_Nu.SetLength(mySimplex.d + 1);
	//sum_De.SetLength(mySimplex.d + 1);
	total = 0;
	lcm = 1;
	//bool repeat[1000]; //put this on the stack, do not waste the time requesting memory from the heap because this function is called many, many, many times.
					 //Why is this bool (vs int): if there are no repeats in the <l, s_i> terms, the simplex is regular on l and we compute the integral as in the first case of the theory.
					 // Otherwise we will have to compute the residue. It is in the residue-function where we worry about the multiplicity of things.

	mat_ZZ mat;
	mat.SetDims(dim, dim); //(n, m) makes mat have dimension n x m.


	//printListCone(cone, dim);
	//find <l, v>^(dim+m).
	scaleRationalVectorToInteger(cone->vertex->vertex, dim, temp_b);
	assert(temp_b == 1);
	cout << "updateLawrence: l = " << l << endl;
	cout << "updateLawrence: numerators: " << cone->vertex->vertex->numerators() << endl;
	cout << "updateLawrence:l^dim+m= ";
	temp_a = l * cone->vertex->vertex->numerators();
	cout << temp_a << "^" << dim << "+ "<< m;
	temp_a = power(temp_a, dim + m);
	cout << "= " << temp_a << endl;
	//temp_b = power(temp_b, dim + m);

	int col = 0;
	for(listVector *ray = cone->rays; ray; ray = ray->rest, col++){

		temp_b *= -1*(ray->first * l); //find <ray, l>
			//why times -1 you ask? The paper says <l, vertex - other vtertex> which is a ray directed TO the vertex, not AWAY from the vertex.
		for (int row = 0; row < dim; row++)
		{
			mat[row][col] = ray->first[row];
		}
	}
	//cout << "hey, I made it past the matrix building: newIntegration::251" << endl;
	if(temp_a < 0 && temp_b < 0){
		temp_a *= -1;
		temp_b *= -1;
	};
	if(temp_b == 0 && temp_a != 0){
		vec_ZZ ProDiff;
		ProDiff.SetLength(dim+1);
		listVector * temp = cone->rays;
		for (i = 0; i < dim; i++)
		{
			//cout << "going to take inner prod " << i << endl;
			//cout << temp->first << endl;
			ProDiff[i] = -1*temp->first * l;
			cout << "ProDiff:" << temp->first << "*" << l << "= " << ProDiff[i] << endl;
			temp = temp->rest;
		}
		ProDiff[dim] = 0;
		//cout << "done with inner prods.!" << endl;
		//computeResidue(dim, m, ProDiff, l * scaleRationalVectorToInteger(
		//	cone->vertex->vertex, dim, temp_b), temp_a, temp_b);

		//TODO: fix computeResidue. questions: how to modify prodDiff and other imputs to trick it into computing the right residue
		//			or we have to do a cute-pate-edit solution which is messay :(
		computeResidue(dim, m, ProDiff, l * cone->vertex->vertex->numerators(), temp_a, temp_b);

	}
	else if ( temp_b == 0 && temp_a == 0)
	{
		temp_b = 1;
	}



	determinant(det, mat);
	cout << "de=" << de << ", det=" << det << ", coe=" << coe << ", tempa= " << temp_a << ", tempb=" << temp_b << endl;
	temp_a *= abs(det) * coe.getNumerator();// we should add to a/b. ???
	temp_b *= de * coe.getDenominator();

	cout << "updateLawrence:total/lcm(L) = " << RationalNTL(temp_a, temp_b) << endl;

	//total/lcm =
	//total = total * mySimplex.v * coe.getNumerator();
	//lcm = lcm * de * coe.getDenominator();


	g = GCD(temp_a, temp_b);
	if (g != 0)
	{
		temp_a = temp_a / g;
		temp_b = temp_b / g;
	}
	//add a/b + temp_a/temp_b

	cout << "updateLawrence:" << a << "/" << b << " + " << temp_a << "/" << temp_b << "==";
	if (a == 0)
	{
		a = temp_a;
		b = temp_b;
	} //return a/b = temp_a/temp_b
	else //if ((lcm != 0) && (b != 0))
	{

		ZZ lcmB;
		lcmB = b * temp_b / GCD(b, temp_b); //find LCM of b and temp_b
		a = a * lcmB / b + temp_a * lcmB / temp_b;
		b = lcmB;
	} // a/b = a/b + temp_a/temp_b.
	cout << a << "/" << b << endl;

} // updateLawrence

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

//This function computes a given fraction a/b, the integral of the linear form forms, over the cone
void integrateLinFormSumLawrence(ZZ& numerator, ZZ& denominator, PolyIterator<RationalNTL, ZZ>* it, listCone *cone, int dim)
{
	ZZ v, de, counter, tem; //, coe;
	RationalNTL coe;
	int i, j, index, k, m;
	vec_ZZ l;
	//if (forms.varCount!=mySimplex.d) {cout<<"The dimensions of the polynomial and simplex don't match. Please check!"<<forms.varCount<<"<>"<<mySimplex.d<<endl;exit(1);};
	l.SetLength(dim);
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
		for (i = 1; i <= dim + m; i++)
		{
			de = de * i;
		} //de is (d+m)!. Note this is different from the factor in the paper because in our storage of a linear form, any coefficient is automatically adjusted by m!
		updateLawrence(numerator, denominator, l, cone, m, coe, de, dim);//We are ready to compute the integral of one linear form over the simplex
		cout << "integrateLinFormSumLawrence:: partial sum:" << numerator << "/" << denominator << endl;
	}
	delete temp;
	if (denominator < 0)
	{
		denominator *= to_ZZ(-1);
		numerator *= to_ZZ(-1);
	}
}//integrateLinFormSum

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

void integrateMonomialSumLawrence(ZZ &a, ZZ &b, monomialSum &monomials,
		listCone *cone, int dim)//integrate a polynomial stored as a Burst Trie
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
	integrateLinFormSumLawrence(a, b, it2, cone, dim);
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

