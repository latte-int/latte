#include "newIntegration.h"
#include "iterators.h"
#include <NTL/vec_ZZ.h>
#include <NTL/ZZ.h>
#include <NTL/mat_ZZ.h>
#include <iostream>

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
		};
	};
}
;

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
}//convert to simplex (ratio


//The purpose of this function is to modify a given fraction a/b. 
//Input: l is the exponents of a linear form, mySimplex is the Simplex we are integrating over with d+1 vertices
//m is the power the linear form is raised to
//coe is the coefficient of a linear form
//de is the extra factor in the formulae that we we multiply the result by
void update(RationalNTL & addToAnswer, const vec_RationalNTL &l,
		const simplexRationalNTL &mySimplex, int m, const RationalNTL &coe,
		const ZZ &de)
{
	ZZ lcm, total, g, tem;
	RationalNTL sum;
	sum.setCanonicalizeFraction(false);
	int i, j;
	vec_RationalNTL inner_Pro; //inner_Pro[i] = <l, s_i>
	vec_RationalNTL sum_Nu, sum_De;
	vec_RationalNTL rationalTerms; // rationalTerms[i] = <l, s_i>^d/ (\prod_{j \neq i} <l, s_i - s_j>)
	inner_Pro.SetLength(mySimplex.d + 1);
	rationalTerms.SetLength(mySimplex.d + 1);
	sum_De.SetLength(mySimplex.d + 1);
	sum_Nu.SetLength(mySimplex.d + 1);
	//sum_Nu.SetLength(mySimplex.d + 1);
	//sum_De.SetLength(mySimplex.d + 1);
	total = 0;
	lcm = 1;
	bool repeat[1000];//*repeat = new bool[mySimplex.d + 1];
	for (i = 0; i <= mySimplex.d; i++)
	{
		//sum = 0;

		//for (j = 0; j < mySimplex.d; j++)
		//	sum = sum + l[j] * mySimplex.s[i][j];
		//inner_Pro[i] = sum; // inner_Pro_i= <l, s_i> where <.,.> is the std. real inner product.
		inner_Pro[i] = vec_RationalNTL::innerProduct(l, mySimplex.s[i]);
		repeat[i] = 0;
		for (j = 0; j < i; j++)
			if (inner_Pro[j] == inner_Pro[i])
			{
				repeat[i] = 1;
				break;
			};//record repetitions
	};//stores inner product for use

	//cout << "l-vector = ";
	//for(int kk = 0; kk < l.length(); ++kk)
	//	cout << l[kk] << ", ";
	//cout << endl;
	//mySimplex.print(cout);
	//for (i = 0; i <= mySimplex.d; i++)
	//{
		//cout << "inner_Pro[" << i << "] = " << inner_Pro[i] << ", repeat" << repeat[i] << endl;
	//}

	for (i = 0; i <= mySimplex.d; i++)
		if (!repeat[i])
		{
			sum_Nu[i] = 1;
			sum_Nu[i].setCanonicalizeFraction(false);
			for (j = 0; j < m + mySimplex.d; j++)
				sum_Nu[i].mult(inner_Pro[i]); // sum_Nu_i = inner_pro ^ (m + d)
			sum_De[i] = 1;
			sum_De[i].setCanonicalizeFraction(false);
			for (j = 0; j <= mySimplex.d; j++)
				if (i != j)
					sum_De[i].mult(inner_Pro[i] - inner_Pro[j]); //sum_de_i	 = \prod _{i != j} <l, s_i - s_j>
			//if ((sum_Nu[i].getNumerator() < 0) && (sum_De[i].getNumerator() < 0))
			//{
			//	sum_Nu[i] = sum_Nu[i] * to_ZZ(-1);
			//	sum_De[i] = sum_De[i] * to_ZZ(-1);
			//};
			if (sum_De[i] == 0)
			{
				vec_RationalNTL ProDiff;
				ProDiff.SetLength(mySimplex.d + 1);
				for (j = 0; j <= mySimplex.d; j++)
					ProDiff[j] = inner_Pro[i] - inner_Pro[j];
				computeResidue(mySimplex.d, m, ProDiff, inner_Pro[i],
						rationalTerms[i]);
				//cout << " computeResidue = " << i << "; " << rationalTerms[i] << endl;
			}
			if (sum_De[i] != 0)
			{
				rationalTerms[i] = sum_Nu[i]/sum_De[i]; //lcm = lcm * sum_De[i] / (GCD(lcm, sum_De[i]));
			};
		};
	//for (i = 0; i <= mySimplex.d; i++)
	//{
	//	if (!(sum_De[i] == 0))
			//cout << "sum_Nu[" << i << "]/ sum_De[i] = " << RationalNTL(sum_Nu[i]/sum_De[i]) << endl;
	//}
	sum = 0;
	//cout << endl;
	for (i = 0; i <= mySimplex.d; i++)
		if ((!repeat[i]) )// && (sum_De[i] != 0))
		{
			sum.add(rationalTerms[i]); //total += sum_Nu[i] * (lcm / sum_De[i]);
			//cout << "rationalTerms[" << i << "] = " << rationalTerms[i] << endl;
			//cout << rationalTerms[i] * coe * RationalNTL(to_ZZ(1), de) * mySimplex.v << " + ";
		}
	//cout << endl;
	//cout << "coe = " << coe << endl;
	sum *= coe;
	sum *= mySimplex.v;
	sum.div(de);
	sum.canonicalize();

	addToAnswer.add(sum);

	//cout << "coe" << coe << endl;
	//cout << "volume =" << mySimplex.v << endl;
	//cout << "de =" << de << endl;
	//cout << "sum = " << sum << endl;
	//cout << "for simplex = " << endl;
	//mySimplex.print(cout);
	//cout << endl;
	/*
	lcm = lcm * de * coe.getDenominator();
	total = total * mySimplex.v * coe.getNumerator(); //was just * coe.
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
	*/
	//delete[] repeat;
}//update (rational).


//The purpose of this function is to modify a given fraction a/b.
//Input: l is the exponents of a linear form, mySimplex is the Simplex we are integrating over with d+1 vertices
//m is the power the linear form is raised to
//coe is the coefficient of a linear form
//de is the extra factor in the formulae that we we multiply the result by

void update(ZZ &a, ZZ &b, const vec_ZZ &l, const simplexZZ &mySimplex, int m,
		const RationalNTL &coe, const ZZ &de)
{
	ZZ sum, lcm, total, g, tem;
	int i, j;
	vec_ZZ inner_Pro, sum_Nu, sum_De;
	inner_Pro.SetLength(mySimplex.d + 1);
	sum_Nu.SetLength(mySimplex.d + 1);
	sum_De.SetLength(mySimplex.d + 1);
	total = 0;
	lcm = 1;
	bool repeat[1000];
	for (i = 0; i <= mySimplex.d; i++)
	{
		sum = 0;
		for (j = 0; j < mySimplex.d; j++)
			sum = sum + l[j] * mySimplex.s[i][j];
		inner_Pro[i] = sum; // inner_Pro_i= <l, s_i> where <.,.> is the std. real inner product.
		repeat[i] = 0;
		for (j = 0; j < i; j++)
			if (inner_Pro[j] == inner_Pro[i])
			{
				repeat[i] = 1;
				break;
			};//record repetitions
	};//stores inner product for use

	//cout << "l-vector = ";
	//for(int kk = 0; kk < l.length(); ++kk)
	//	cout << l[kk] << ", ";
	//cout << endl;
	//mySimplex.print(cout);

	//for (i = 0; i <= mySimplex.d; i++)
	//{
		//cout << "inner_Pro[" << i << "] = " << inner_Pro[i] << ", repeat" << repeat[i] << endl;
	//}

	for (i = 0; i <= mySimplex.d; i++)
		if (!repeat[i])
		{
			sum_Nu[i] = 1;
			for (j = 0; j < m + mySimplex.d; j++)
				sum_Nu[i] = sum_Nu[i] * inner_Pro[i]; // sum_Nu_i = inner_pro ^ (m + d)
			sum_De[i] = 1;
			for (j = 0; j <= mySimplex.d; j++)
				if (i != j)
					sum_De[i] = sum_De[i] * (inner_Pro[i] - inner_Pro[j]); //sum_de_i	 = \prod _{i != j} <l, s_i - s_j>
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
				//cout << " computeResidue < " << i << "; " << sum_Nu[i] << '/' << sum_De[i] << endl;
			}
			if (sum_De[i] != 0)
			{
				lcm = lcm * sum_De[i] / (GCD(lcm, sum_De[i]));
			};
		};

	//for (i = 0; i <= mySimplex.d; i++)
	//{
	//	if  ( sum_De[i] != 0)
			//cout << "sum_Nu[" << i << "]/ sum_De[i] = " << RationalNTL(sum_Nu[i], sum_De[i]) << endl;
	//}

	//cout << endl;
	for (i = 0; i <= mySimplex.d; i++)
		if ((!repeat[i]) && (sum_De[i] != 0))
		{
			total += sum_Nu[i] * (lcm / sum_De[i]);


			//cout << RationalNTL(sum_Nu[i], sum_De[i]) * coe * RationalNTL(to_ZZ(1), de) * mySimplex.v << " + ";
		}
	//cout << endl;

	lcm = lcm * de * coe.getDenominator();
	total = total * mySimplex.v * coe.getNumerator(); //was just * coe.
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
	//cout << "coe < " << coe << endl;
	//cout << "coe" << coe << endl;
	//cout << "volume =" << mySimplex.v << endl;
	//cout << "de =" << de << endl;
	//cout << "sum < " << total/GCD(total, lcm) << "/" << lcm/GCD(total, lcm) << " / " << newIntegrationFile::debugDilationFactor << "^" << mySimplex.d
	//	 << "= " << RationalNTL(total, lcm).div(power(newIntegrationFile::debugDilationFactor, mySimplex.d)) << endl;
	//cout << "for simplex = " << endl;
	//mySimplex.print(cout);
	//cout << endl;
}//update (integer)


//This function computes a given fraction a/b, the integral of the linear form forms, over the simplex mySimplex
void integrateLinFormSum(RationalNTL & answer,
		PolyIterator<RationalNTL, ZZ>* it, const simplexRationalNTL &mySimplex)
{
	ZZ v, de, counter, tem; //, coe;
	RationalNTL coe;
	int i, j, index, k, m;
	vec_RationalNTL l;
	//if (forms.varCount!=mySimplex.d) {cout<<"The dimensions of the polynomial and simplex don't match. Please check!"<<forms.varCount<<"<>"<<mySimplex.d<<endl;exit(1);};
	l.SetLength(mySimplex.d);
	answer = 0;

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
		}; //de is (d+m)!. Note this is different from the factor in the paper because in our 								storage of a linear form, any coefficient is automatically adjusted by m!
		update(answer, l, mySimplex, m, coe, de);//We are ready to compute the integral of one linear form over the simplex
	};
	delete temp;
}//integrateLinFormSum (rational)

//This function computes a given fraction a/b, the integral of the linear form forms, over the simplex mySimplex
void integrateLinFormSum(ZZ& numerator, ZZ& denominator, PolyIterator<
		RationalNTL, ZZ>* it, const simplexZZ &mySimplex)
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
		}; //de is (d+m)!. Note this is different from the factor in the paper because in our 								storage of a linear form, any coefficient is automatically adjusted by m!
		update(numerator, denominator, l, mySimplex, m, coe, de);//We are ready to compute the integral of one linear form over the simplex
	};
	delete temp;
	if (denominator < 0)
	{
		denominator *= to_ZZ(-1);
		numerator *= to_ZZ(-1);
	};
}//integrateLinFormSum (integer)


void integrateMonomialSum(ZZ &a, ZZ &b, monomialSum &monomials,
		const simplexZZ &mySimplex)//integrate a polynomial stored as a Burst Trie
{
	linFormSum forms;
	forms.termCount = 0;
	forms.varCount = monomials.varCount;
	BTrieIterator<RationalNTL, int>* it =
			new BTrieIterator<RationalNTL, int> ();
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

