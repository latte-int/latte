/*
 * TopEhrhart.cpp
 *
 *  Created on: Jan 2, 2012
 *      Author: Koeppe, Brandon
 */

#include "TopEhrhart.h"

TopEhrhart::TopEhrhart(Polyhedron * polyhedron, BarvinokParameters & para, int numTopCoeff, bool real, string saveTopEhrhartPoly ) :
	poly(polyhedron), parameters(para), numTopCoefficients(numTopCoeff), realDilations(real), saveTopEhrhartPolynomial(saveTopEhrhartPoly)
{
	assert(poly != NULL);
	assert(poly->cones != NULL);
	assert(poly->cones->rays != NULL);
	assert(poly->homogenized == false);
	assert(poly->dualized == false);

	if (numTopCoefficients != -1 && numTopCoefficients <= 0)
		THROW_LATTE_MSG(LattException::ue_BadCommandLineOption, "unexpedted numTopCoefficients given");
}

TopEhrhart::~TopEhrhart()
{
	// TODO Auto-generated destructor stub
}

/**
 * Calls a maple script to compute the weighted Ehrhart polynomial.
 * @parm polynomial: weight each integer by a polynomial.
 * The polynomial is decomposed to linear forms.
 */
void TopEhrhart::computeTopEhrhartPolynomial(const monomialSum & polynomial)
{
	linFormSum linearForms;

	//make an iterator for the transformed polynomial and decompose it into linear forms.
	BTrieIterator<RationalNTL, int>* polynomialItr = new BTrieIterator<RationalNTL, int> ();
	linearForms.termCount = 0;
	linearForms.varCount = polynomial.varCount;

	//decompose to linear forms. Constant polynomials are decomposed to [0,0,00,0]^0

	polynomialItr->setTrie(polynomial.myMonomials, polynomial.varCount);
	decompose(polynomialItr, linearForms);

	computeTopEhrhartPolynomial(linearForms);

	destroyLinForms(linearForms);
	delete polynomialItr;
}

/**
 * Calls a maple script to compute the weighted Ehrhart polynomial.
 * @parm linForm: weight each integer by sum of powers of linear forms.
 */
void TopEhrhart::computeTopEhrhartPolynomial(const linFormSum & linForm)
{

	ofstream maple("compute-top-ehrhart.mpl");

	maple << "read(\"" << relocated_pathname(MAPLE_SCRIPT_DIR) << "/"
			<< "Conebyconeapproximations_08_11_2010.mpl" << "\"):\n\n";

	maple << "\n seed:=randomize():" << endl;

	//print the polytope's vertex-rays in maple format out.
	maple << "\n simpleCones := [";
	for ( listCone *cone = poly->cones; cone; cone = cone->rest)
	{
		// add [[vertex], [[ray1], [ray2], ..., [ray d]] ]

		//add the vertex
		maple << "\n[";
		vec_ZZ num = cone->vertex->vertex->numerators();
		vec_ZZ den = cone->vertex->vertex->denominators();
		maple << "[";
		for(int i = 0; i < poly->numOfVars; ++i)
		{
			maple << RationalNTL(num[i], den[i]);
			if ( i != poly->numOfVars -1)
				maple << ", ";
		}
		maple << "], ";
		//done adding the vertex.

		//check the polytope is simple.
		assert(lengthListVector(cone->rays) == poly->numOfVars);


		//start of ray list
		maple << "[";
		for(listVector * ray = cone->rays; ray; ray = ray->rest)
		{
			maple << "\n\t[";
			for(int i = 0; i < poly->numOfVars; ++i)
			{
				maple << (ray->first)[i];
				if (i < poly->numOfVars -1)
					maple << ", ";
			}
			maple << "]";
			if ( ray->rest)
				maple << ",";
		}
		maple << "]";
		//end of ray list


		maple << "]";
		if (cone->rest)
			maple << ",";
	}//for each vertex.
	maple << "]: #end of the vertex-cones" << endl;

	//Ok, now print the linear form list.
	//I inlined the printLinForms(linForm) function here because
	// I don't want to print the coefficients with an extra (degree)! term
	BTrieIterator<RationalNTL, ZZ>* it = new BTrieIterator<RationalNTL, ZZ> ();
	term<RationalNTL, ZZ>* temp;
	it->setTrie(linForm.myForms, linForm.varCount);
	it->begin();

	maple << "linearForms:= [";
	for(temp = it->nextTerm(); temp; )
	{
		ZZ degFactorial;
		degFactorial = 1;
		for(int j = temp->degree; j > 1; --j)
			degFactorial *= j;

		maple << "[" << temp->coef / degFactorial << ", [" << temp->degree << ", [";
		for (int j = 0; j < temp->length; j++)
		{
			maple << temp->exps[j];
			if (j + 1 < temp->length)
			{
				maple << ", ";
			}
		}
		maple << "]]]";
		temp = it->nextTerm();
		if (temp)
			maple << ", \n\t";
	}
	maple << "]: #end of the linear forms. \n\n";
	delete it;

	//end of printing the linear form list.

	//now print the function we are going to call.
	if(realDilations)
	{
		if (numTopCoefficients >= 0)
			maple << "epoly:=printIncrementalEhrhartPolynomial(n,N,simpleCones,linearForms," << poly->numOfVars << ", true, "<< numTopCoefficients << ", \"" << saveTopEhrhartPolynomial.c_str() << "\"):" << endl;
		else
			maple << "epoly:=printIncrementalEhrhartPolynomial(n,N,simpleCones,linearForms," << poly->numOfVars << ", true, -1, \"" << saveTopEhrhartPolynomial.c_str() << "\"):" << endl;
	}
	else
	{
		if (numTopCoefficients >= 0)
			maple << "epoly:=printIncrementalEhrhartPolynomial(n,N,simpleCones,linearForms," << poly->numOfVars << ", false, "<< numTopCoefficients << ", \"" << saveTopEhrhartPolynomial.c_str() << "\"):" << endl;
		else
			maple << "epoly:=printIncrementalEhrhartPolynomial(n,N,simpleCones,linearForms," << poly->numOfVars << ", false, -1, \"" << saveTopEhrhartPolynomial.c_str() <<"\"):" << endl;
	}//integer dilations.

	maple << "\n\n";
	maple << "printf(\"\\n\\n\\nThe Ehrhart polynomial=%a\\n\", epoly);\n" << endl;
	maple << "printf(\"Evaluation at 1 is %a\\n\", eval(subs({N=1,n=1, MOD=latteMod},epoly)));\n" << endl;
	maple.close();

	system_with_error_check(MAPLE_PATH + string(" -q compute-top-ehrhart.mpl"));
}

/**
 * Calls a maple script to compute the unweighted Ehrhart polynomial.
 */
void TopEhrhart::computeTopEhrhartPolynomial()
{
	//insert the linear form (0,0,0,0...)^0 = 1
	linFormSum forms;
	FormLoadConsumer<RationalNTL> *consumer =
			new FormLoadConsumer<RationalNTL> ();
	consumer->setFormSum(forms);

	consumer->setDimension(poly->numOfVars);
	vec_ZZ coefs;
	coefs.SetLength(poly->numOfVars);
	for (int i = 0; i < poly->numOfVars; ++i)
		coefs[i]=0;

	RationalNTL coefficient;
	coefficient = 1;

	consumer->ConsumeLinForm(coefficient, 0, coefs);
	delete consumer;

	computeTopEhrhartPolynomial(forms);
}

