/*
 * Perturbation.cpp
 *
 *  Created on: Nov 4
 *      Author: Brandon Dutra and Gregory Pinto
 */

#include <cassert>
#include "Perturbation.h"
#include "print.h"

using namespace std;

void LinearPerturbationContainer::setListCones(int dim,
		listCone * simpleConeList)
{

	coneTerms.resize(lengthListCone(simpleConeList));
	currentPerturbation.SetLength(dim);
	dimention = dim;

	//build the vector of cones.
	listCone* ptr = simpleConeList;
	for (unsigned int i = 0; i < coneTerms.size(); ++i)
	{
		coneTerms[i].setSimplicialCone(ptr, dim);
		ptr = ptr->rest;
	}//for every cone.
}

bool LinearPerturbationContainer::tryCurrentPerturbation(const vec_ZZ &l)
{
	divideByZero = false;
	for (unsigned int i = 0; i < coneTerms.size(); ++i)
		if ((divideByZero = coneTerms[i].computeDotProducts(
				currentPerturbation, l)))
			break;
	return divideByZero;
}//tryCurrentPerturbation

bool LinearPerturbationContainer::tryNoPerturbation(const vec_ZZ &l)
{
	divideByZero = false;
	for (unsigned int i = 0; i < coneTerms.size(); ++i)
		if (coneTerms[i].computeDotProducts(l))
			divideByZero = true; //don't break b/c want to find all inner products w/o pertrubation.
	return divideByZero;
}

void LinearPerturbationContainer::findPerturbation(const vec_ZZ &l)
{

	ZZ gcdValue;
	divideByZero = false;
	if (tryNoPerturbation(l) == false)
	{
		//cout << "findPerturbation::" << currentPerturbation << endl;
		//for(long i = 0; i < coneTerms.size(); ++i)
		//	coneTerms[i].printTerm();
		//cout <<" end printing perturbed system" << endl;
		return; //no errors. doing nothing worked!
	}

	currentPerturbation = l;
	currentPerturbation[rand() % currentPerturbation.length()] += 1+ rand() % 10;

	//project p onto l: <p,l>/||l|| * l but we do not care about the scale, so don't divide.
	currentPerturbation = currentPerturbation - (currentPerturbation * l) * l;
	//currentPerturbation is now normal to l.

	//see if we can scale the currentPerturbation down
	gcdValue = currentPerturbation[0];
	for (long i = 1; i < currentPerturbation.length(); ++i)
	{
		gcdValue = GCD(gcdValue, currentPerturbation[i]);
	}
	if (gcdValue > 1)
		for (long i = 1; i < currentPerturbation.length(); ++i)
			currentPerturbation[i] = currentPerturbation[i] / gcdValue;

	while (tryCurrentPerturbation(l) == true)
	{
		//we divided by zero again. that is <l+e, r>=0 for some ray and for for some cone while <l,v>!=0. :(

		//try a new perturbation.
		currentPerturbation[rand() % currentPerturbation.length()] += 1+rand()
				% 10;

		//see if we can scale things down.
		gcdValue = currentPerturbation[0];
		for (long i = 1; i < currentPerturbation.length(); ++i)
		{
			gcdValue = GCD(gcdValue, currentPerturbation[i]);
		}
		if (gcdValue > 1)
			for (long i = 1; i < currentPerturbation.length(); ++i)
				currentPerturbation[i] = currentPerturbation[i] / gcdValue;
	}//make a new perturbation and try again.

	//cout << "findPerturbation::" << currentPerturbation << endl;
	//for(long i = 0; i < coneTerms.size(); ++i)
	//	coneTerms[i].printTerm();
	//cout <<" end printing perturbed system" << endl;


}//findPerturbation

/**integrates the polytope over l.
 *
 */
RationalNTL LinearPerturbationContainer::integratePolytope(int m)
{
	RationalNTL totalSum; //contains the sum. init. to zero.

	for (unsigned int i = 0; i < coneTerms.size(); ++i)
		coneTerms[i].integrateTerm(totalSum, m, dimention);

	if (dimention %2)
		totalSum.changeSign();
	return totalSum; //integral of l over the polytope!
}

//--------------------------------------
//--------------------------------------
//START OF LinearLawrenceIntegration
//--------------------------------------
//--------------------------------------

LinearLawrenceIntegration::LinearLawrenceIntegration() :
	divideByZero(0), simplicialCone(NULL)
{

}

LinearLawrenceIntegration::LinearLawrenceIntegration(listCone * cone) :
	divideByZero(0), simplicialCone(cone)
{

}

void LinearLawrenceIntegration::setSimplicialCone(listCone *cone, int dim)
{
	simplicialCone = cone;
	rayDotProducts.resize(dim);


}

/** Finds the innerproducts of the vertex and rays for each term.
 *  If a residue calculation is needed, returns true, otherwise
 *  we did not divide by zero and false is returned.
 */
bool LinearLawrenceIntegration::computeDotProducts(const vec_ZZ &l)
{
	//polytope should be dilated...should be integer!
	const vec_ZZ &vertex = simplicialCone->vertex->vertex->numerators();

	//update the vertex dot products. if <v,l>=0, we are done processing this cone
	//and so report that we did not divide by zero.
	numeratorDotProduct.constant = vertex * l;
	if (numeratorDotProduct.constant == 0)
	{
		divideByZero = false;
		numeratorDotProduct.epsilon = 0;
		numeratorDotProduct.power = 0;
		return false;
	}

	//for every ray, find <l, r>. report true if we divided by zero.
	unsigned int i = 0;
	divideByZero = false;
	//cout << "printing cone" << endl;
	//printCone(simplicialCone, l.length());


	//now find the abs. value of the det.
	mat_ZZ rayMatrix;
	rayMatrix.SetDims(rayDotProducts.size(), rayDotProducts.size());//dimention.
	for (listVector * ray = simplicialCone->rays; ray; ray = ray->rest, ++i)
	{
		//cout << "i=" << i << endl;

		rayDotProducts[i].constant = l * (ray->first);
		rayDotProducts[i].epsilon = 0;
		rayDotProducts[i].power = 0;

		if (rayDotProducts[i].constant == 0)
			divideByZero = true;

		for(unsigned int j = 0; j < rayDotProducts.size(); ++j)
			rayMatrix[i][j]  = ray->first[j];
	}//for each ray.
	determinant = abs(NTL::determinant(rayMatrix));
	//cout << "going to return " << divideByZero << endl;
	return divideByZero;
}

/** Finds the innerproducts of the vertex and rays for each term with the perturbation term ONLY.
 *  If a residue calculation is needed, returns true, otherwise
 *  we did not divide by zero and false is returned.
 */
bool LinearLawrenceIntegration::computeDotProducts(const vec_ZZ &e,
		const vec_ZZ &l)
{
	if (divideByZero == false)
		return false; //we do not need to do further process. plugging in numbers will work.

	//find that rays that dotted to zero and see if this new perturbation will work.
	unsigned int i = 0;
	for (listVector * ray = simplicialCone->rays; ray; ray = ray->rest, ++i)
	{
		rayDotProducts[i].epsilon = e * ray->first;
		if (rayDotProducts[i].constant == 0 && rayDotProducts[i].epsilon == 0)
			return true; //darn, this perturbation did not fix our problem. stop processing.
	}
	//if we got here then:
	//	1) <l,r>=0 for some rays r. and for those rays, <l+e,r>!=0.
	//  2) <l,vertex>!=0.

	//now update the numerator term.
	const vec_ZZ & vertex = simplicialCone->vertex->vertex->numerators();
	numeratorDotProduct.epsilon = e * vertex;
	//numeratorDotProduct.constant is not zero otherwise divideByZero would have been false.

	//we are done checking this cone's dot products.
	//NOTE::we still did not remove repeated dot products...the powers are all still zero.
	return false;
}

/**
 * Integrates l over this cone. Adds the result to totalSum
 */
void LinearLawrenceIntegration::integrateTerm(RationalNTL &totalSum, int m,
		int dim)
{
	ZZ numerator, denominator;

	//cout << endl;
	denominator = 1;
	if (numeratorDotProduct.constant == 0)
		return; //integral over cone is zero!!!

	if (divideByZero == false)
	{
		numerator = power(numeratorDotProduct.constant, m + dim);
		for (unsigned int i = 0; i < dim; ++i)
			denominator *= rayDotProducts[i].constant;

		totalSum.add(numerator * determinant, denominator);
		//cout << "--compute redidue not needed."
		//	 << "det * term = " << determinant << "* " << RationalNTL(numerator, denominator);
		return;
	}//just plug in numbers.


	//need to find residue after increasing the power of e...or the constant term in the series expansion
	//first, sort/merge repeated terms. so (c1x + c0)*(c1x + c0)=(c1x + c0)^2
	//(actually, we do (c1x + c0)*(c1x + c0)=(c1x + c0)^2(c1x + c0)^-1 and we just ignore negative powers.
	//  so we treat the power as a flag to denote now the term has been processed. 0 = not processed, n = regular power, -1=repeated term.

	//cout << "before update power: ";
	//printTerm(true);
	updatePowers(); //also the location of the (0+e)^{m_0} term is in array index 0.
	//cout << "after update power: ";
	//printTerm();

	//cout << "--going to call computeResidueLawrence(" << dim << ", " << m <<" this, " << numerator << ", " <<  denominator << endl;
	computeResidueLawrence(dim, m, *this, numerator, denominator);
	//cout << "--returned from computeResidueLawrence ok" << endl;
	//cout << "det * term = " << determinant << "* " << RationalNTL(numerator, denominator);
	totalSum.add(numerator * determinant, denominator);

	//todo: write a friend function in residue.cpp (not in /valuation/) that will compute the residue using
	//our data structure.

}//integrateTerm


//prints out the current term.
//Afer updatePowers(), the term of  the zero constants has a special form.
//Intead of (0 + ce)^m, the zero constants should be thought as 0 + c (e^m).
void LinearLawrenceIntegration::printTerm(bool printVertexRays) const
{
	cout << "(" << numeratorDotProduct.constant << "+ " << numeratorDotProduct.epsilon << "e)^" << numeratorDotProduct.power
	     << "/ ";
	for(unsigned int i = 0; i < rayDotProducts.size(); ++i)
		cout << "(" << rayDotProducts[i].constant << " + " << rayDotProducts[i].epsilon << "e)^" << rayDotProducts[i].power << " ";
	if ( true == printVertexRays)
	{
		cout << endl;
		cout << "  vertex: " << simplicialCone->vertex->vertex->numerators() << endl;
		for(listVector *v = simplicialCone->rays; v; v = v->rest)
			cout << "  ray: " << v->first << endl;
	}
	cout << endl;
}//printTerm()


/** Merges powers and returns location of the (0+e) term.
 *  Assumes a residue really need to be called. That is, assumes a (0+e) term exists and will return it's location...or else will return -1.
 *  After this, if power = 0: there should be no 0 powers!
 *  					   k: then the power of this term is k.
 *  					  -1: There were repeated terms and this term has been merged with another term.
 *  example: (c1x + c0)*(b1x + b0)*(c1x + c0)=(c1x + c0)^2*(b1x + b0)*(c1x + c0)^-1
 */
void LinearLawrenceIntegration::updatePowers()
{
	int locationOfe = -1;

	//cout << "rayDotProd.size()" << rayDotProducts.size() << endl;
	//printTerm();

	unsigned int i, j;


	for (i = 0; i < rayDotProducts.size(); ++i)
	{
		if (rayDotProducts[i].power < 0)
			continue; //already processed this term.

		//if found a (0+c*e) term.
		if (rayDotProducts[i].constant == 0)
		{
			if (locationOfe == -1)
			{
				locationOfe = i;
				rayDotProducts[locationOfe].power = 1;
			} //if first time we found a (0+c*e) term.
			else
			{

				rayDotProducts[locationOfe].epsilon	*= rayDotProducts[i].epsilon;
				rayDotProducts[locationOfe].power++;
				rayDotProducts[i].power = -1;
			}//add to the (0+c*e) term.
			continue;
		}//handle (0+c*e) terms differently.


		int numberRepeatedTerms = 1;//count ourself, and start j at i+1
		for (j = i + 1; j < rayDotProducts.size(); ++j)
		{
			if (rayDotProducts[j].constant == rayDotProducts[i].constant
					&& rayDotProducts[j].epsilon == rayDotProducts[i].epsilon
					&& rayDotProducts[j].power == 0)
			{
				++numberRepeatedTerms;
				rayDotProducts[j].power = -1;//don't consider this term again.
			}//if (c0 + c1e) == (b0 + b1e) and (b0 + b1e) has not already been processed.
		}//for j.
		rayDotProducts[i].power = numberRepeatedTerms;
	}//for i.

	//printTerm();
	assert(locationOfe >= 0); //otherwise we would not have called this function!

	if (locationOfe > 0)
	{
		linearPerturbation temp = rayDotProducts[0];
		rayDotProducts[0] = rayDotProducts[locationOfe];
		rayDotProducts[locationOfe] = temp;
	}//swap locationOfE with index 0

}//updatePowers


