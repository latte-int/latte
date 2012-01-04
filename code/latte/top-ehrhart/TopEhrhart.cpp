/*
 * TopEhrhart.cpp
 *
 *  Created on: Jan 2, 2012
 *      Author: Koeppe, Brandon
 */

#include "TopEhrhart.h"

TopEhrhart::TopEhrhart(Polyhedron * polyhedron, BarvinokParameters & para, int numTopCoeff, bool real) :
	poly(polyhedron), parameters(para), numTopCoefficients(numTopCoeff), realDilations(real)
{
	assert(poly != NULL);
	assert(poly->cones != NULL);
	assert(poly->cones->rays != NULL);
	assert(poly->homogenized == false);
	assert(poly->dualized == false);
}

TopEhrhart::~TopEhrhart()
{
	// TODO Auto-generated destructor stub
}

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


void TopEhrhart::computeTopEhrhartPolynomial(const linFormSum & linForm)
{
	//check the polytope is simple.

	ofstream maple("compute-top-ehrhart.mpl");

	maple << "read(\"" << relocated_pathname(MAPLE_SCRIPT_DIR) << "/"
			<< "Conebyconeapproximations_08_11_2010.mpl" << "\"):\n\n";

	maple << "\n seed:=randomize();" << endl;

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

		//start of ray list
		maple << "[";
		assert(lengthListVector(cone->rays) == poly->numOfVars);
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
	maple << "]; #end of the vertex-cones" << endl;

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
	maple << "];\n\n";
	delete it;

	//end of printing the linear form list.

	//now print the function we are going to call.
	if(realDilations)
	{
		if (numTopCoefficients >= 0)
			maple << "epoly:=printIncrementalEhrhartPolynomial(n,N,simpleCones,linearForms," << poly->numOfVars << ", true, "<< numTopCoefficients << "):" << endl;
		else
			maple << "epoly:=printIncrementalEhrhartPolynomial(n,N,simpleCones,linearForms," << poly->numOfVars << ", true, -1):" << endl;
	}
	else
	{
		if (numTopCoefficients >= 0)
			maple << "epoly:=printIncrementalEhrhartPolynomial(n,N,simpleCones,linearForms," << poly->numOfVars << ", false, "<< numTopCoefficients << "):" << endl;
		else
			maple << "epoly:=printIncrementalEhrhartPolynomial(n,N,simpleCones,linearForms," << poly->numOfVars << ", false, -1):" << endl;
	}//integer dilations.

	maple << "printf(\"The Ehrhart polynomial=%a\", epoly);\n" << endl;
	maple.close();

	system_with_error_check(MAPLE_PATH + string(" compute-top-ehrhart.mpl"));
}
