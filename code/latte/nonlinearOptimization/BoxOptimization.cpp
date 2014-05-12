/*
 * BoxOptimizationParameters.cpp
 *
 *  Created on: Apr 23, 2014
 *      Author: bedutra
 */

#include "BoxOptimization.h"
#include "ramon.h"
#include "nonlinearOptimization/WeightedExponentialSubs.h"

WeightedBoxProducer::WeightedBoxProducer(const vec_ZZ & lowerb, const vec_ZZ & upperb):
		lowerBound(lowerb), upperBound(upperb)
{
}

WeightedBoxProducer::~WeightedBoxProducer() {}


void WeightedBoxProducer::Produce(ConeConsumer &consumer)
{
	unsigned int n = upperBound.length();
	ZZ two = to_ZZ(2);
	ZZ twon = power(two, n);//2^n

	vec_ZZ allOnes, theVertex;
	allOnes.SetLength(n);
	theVertex.SetLength(n);
	for(int i = 0; i < n; ++i)
		allOnes[i] = 1;

	mat_ZZ id; // n by n idenity matrix.
	ident(id, n);


	for(ZZ i = to_ZZ(0); i < twon; ++i)
	{

		listCone *oneCone = createListCone();
		oneCone->determinant = 1;
		oneCone->dual_determinant = 1;

		//use the binary representation of i to pick the current vertex along with its tangent cone.
		int det = 0;
		for(int j = 0; j < n; ++j)
		{
			if ( bit(i, j) )
			{
				theVertex[j] = upperBound[j];
				oneCone->rays = appendVectorToListVector(id[j]*(-1), oneCone->rays);
				//oneCone->facets = appendVectorToListVector(id[j], oneCone->facets);
				det++;
			}
			else
			{
				theVertex[j] = lowerBound[j];
				oneCone->rays = appendVectorToListVector(id[j], oneCone->rays);
				//oneCone->facets = appendVectorToListVector(id[j]*(-1), oneCone->facets);
			}
		}

		//determinant = -1^(n/2) * -1^(det)
		//                |           |
		//                |         (number of negative rays)
		//             (sorting columns to identity)
		if ( (det + n/2) % 2)
			oneCone->determinant = -1;
		if ( (n/2 + n-det) % 2)
			oneCone->dual_determinant = -1;

		oneCone->vertex = new Vertex(new rationalVector(theVertex, allOnes));

		//oneCone->facet_divisors = allOnes; //need this if computing facets too.

		oneCone->latticePoints = createListVector(theVertex); //cone is unimodular

		//finally, process this tangent cone.
		consumer.ConsumeCone(oneCone);
	}//for each vertex of a box.
}//WeightedBoxProducer::Produce








/**
 *
 *
 * assumes lowerbound[i] < upperBound[i]
 */
mpq_class computeWeightedCountingBox(const vec_ZZ &lowerBound, const vec_ZZ &upperBound, const linFormSum &originalLinearForms)
{
	mpq_class ans;

	ans = 0;

	BTrieIterator<RationalNTL, ZZ>* linearFormsItr = new BTrieIterator<RationalNTL, ZZ> ();
	linearFormsItr->setTrie(originalLinearForms.myForms, originalLinearForms.varCount);
	linearFormsItr->begin();


	term<RationalNTL, ZZ>* lform;

	//loop over the linear forms
	for (lform = linearFormsItr->nextTerm(); lform; lform
			= linearFormsItr->nextTerm())
	{
		//project the box into a smaller dimension if one of the linear form coeffs is zero.
		vec_ZZ linFormCoeffs, newLB, newUB;
		linFormCoeffs.SetLength(lform->length);
		newLB.SetLength(lform->length);
		newUB.SetLength(lform->length);

		ZZ zeroTerms;
		zeroTerms = 1;

		int j = 0;
		for(int i = 0; i < lform->length; ++i)
		{
			if ( lform->exps[i] != 0)
			{
				linFormCoeffs[j] = lform->exps[i];
				newLB[j] = lowerBound[i];
				newUB[j] = upperBound[i];
				++j;
			}
			else
				zeroTerms *= (upperBound[i] - lowerBound[i] + 1);
		}
		linFormCoeffs.SetLength(j); //ignore the zeros at the end if any.
		newUB.SetLength(j);
		newLB.SetLength(j);

		//easy special cases.
		if ( lform->coef == 0 ||  (j == 0 && lform->degree != 0))
			continue; // if coeff is zero or the linear form is the zero polynomial.
		if ( lform->degree == 0)
		{
			ZZ t;
			t = 1;
			for(int i = 0; i < lowerBound.length(); ++i)
				t *= (upperBound[i] - lowerBound[i] + 1);
			ans += mpq_class(convert_ZZ_to_mpz(lform->coef.getNumerator()), convert_ZZ_to_mpz(lform->coef.getDenominator())) * convert_ZZ_to_mpz(t);
			continue;
		}//if the linear form is a constant

		//cout << "newUB" << newUB << endl;
		//cout << "newLB" << newLB << endl;


		//construct a producer-consumer pair.
		WeightedBoxProducer wbp(newLB, newUB);
		Weighted_Exponential_Single_Cone_Parameters wParams;

		wParams.substitution = BarvinokParameters::PolynomialSubstitution;
		wParams.decomposition = BarvinokParameters::DualDecomposition;
		wParams.max_determinant = 1;
		wParams.Number_of_Variables = j;
		wParams.InitializeComputation();
		wParams.linForm = linFormCoeffs;
		wParams.linFormPow = lform->degree;
		wbp.Produce(wParams);

		ans += mpq_class(convert_ZZ_to_mpz(lform->coef.getNumerator()), convert_ZZ_to_mpz(lform->coef.getDenominator())) * wParams.result * convert_ZZ_to_mpz(zeroTerms);

	}//for every term in the originalPolynomial

	delete linearFormsItr;


	return ans;
}

