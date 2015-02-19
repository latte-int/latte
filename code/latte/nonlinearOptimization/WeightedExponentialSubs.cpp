/*
 * WeightedExponentialSubs.cpp
 *
 *  Created on: Apr 24, 2014
 *      Author: bedutra
 */


#include <iostream>
#include "WeightedExponentialSubs.h"
#include "ExponentialSubst.h"
#include "todd/todd-expansion.h"
#include "genFunction/piped.h"
#include "ramon.h"
#include "print.h"
#include "LattException.h"

using namespace std;

int Weighted_Exponential_Single_Cone_Parameters::ConsumeCone(listCone *cone)
{
	//BarvinokParameters p;
	mpq_class r;
	r = computeExponentialResidue_Single_boxOpt(*wcb, generic_vector, cone, linForm.length(), this, linForm, linFormPow);
	//r = computeExponentialResidue_Single_boxOpt_clean(generic_vector, cone, linForm.length(), this, linForm, linFormPow);
	//cout << "valuation on one cone:" << r << endl;
	result += r;
	freeCone(cone);
	return 1;

}

void Weighted_Exponential_Single_Cone_Parameters::InitializeComputation()
{
  Generic_Vector_Single_Cone_Parameters::InitializeComputation();
  result = 0;
}

//*****************************
int Weighted_Exponential_Single_Cone_Parameters_BranchBound::ConsumeCone(listCone *cone)
{

    computeExponentialResidueWeights_boxOpt(*wcb, generic_vector, cone, linForm.length(), linForm,linFormPow);
    table->weights[index] = wcb->weights;
	++index;
	freeCone(cone);
	return 1;
}

void Weighted_Exponential_Single_Cone_Parameters_BranchBound::InitializeComputation()
{
  Generic_Vector_Single_Cone_Parameters::InitializeComputation();
  index = 0;
  table = new WeightedExponentialTable();
  table->weights.resize(pow(2, Number_of_Variables));
}





mpq_class computeWeightedExponentialResidue(listCone *cone, BarvinokParameters *params, linFormSum &originalLinearForms)
{
	mpq_class ans;
	ans = 0;
	WeightedCountingBuffer wcb;
	
	BTrieIterator<RationalNTL, ZZ>* linearFormsItr = new BTrieIterator<RationalNTL, ZZ> ();
	linearFormsItr->setTrie(originalLinearForms.myForms, originalLinearForms.varCount);
	linearFormsItr->begin();


	term<RationalNTL, ZZ>* lform;
	
	//loop over the original polynomial, and insert the updated monomial into the transformedPolynomial
	for (lform = linearFormsItr->nextTerm(); lform; lform
			= linearFormsItr->nextTerm())
	{

		vec_ZZ linFormCoeffs;
		linFormCoeffs.SetLength(params->Number_of_Variables);
		
		for(int i = 0; i < lform->length; ++i)
			linFormCoeffs[i] = lform->exps[i];

		while( true)
		{
			try
			{
				vec_ZZ generic_vector = guess_generic_vector(params->Number_of_Variables);
				ans += convert_RationalNTL_to_mpq(lform->coef) * computeWeightedExponentialResidue_singleForm(wcb, cone, params, linFormCoeffs, generic_vector, lform->degree);
				break;
			}
			catch (NotGenericException)
			{
				THROW_LATTE_MSG( LattException::bug_NotImplementedHere, "linear form is not generic. Perturbation is not implemented yet" );
			}
		}//while divided by zero, try again.

	
	}//for every term in the originalPolynomial

	delete linearFormsItr;
	
	
	return ans;
}


mpq_class computeWeightedExponentialResidue_singleForm(WeightedCountingBuffer & wcb, listCone *cone, BarvinokParameters *params, const vec_ZZ &linFormCoeffs, const  vec_ZZ &generic_vector, int M)
{
	mpq_class ans;
	ans = 0;

	for(listCone * c = cone; c; c = c->rest)
	{
		ans += computeExponentialResidue_Single_boxOpt(wcb, generic_vector, c, params->Number_of_Variables, params, linFormCoeffs, M);
	}
	
	return ans;
}
