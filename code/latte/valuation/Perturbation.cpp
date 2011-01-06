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


/**
 * Allocates space for each term in the integration formula.
 * @dim: dimention of the polytope (number of slots in each vertex).
 * @simpleConeList: linked list of every simple cone.
 */
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

/**
 * Computes the dot products of <e, l> where e is the current pertrubation.
 * If we still divide by zero, we break and return true;
 * @l: the same linear form  tryNoPerturbation() was called with.
 */

bool LinearPerturbationContainer::tryCurrentPerturbation(const vec_ZZ &l)
{
	divideByZero = false;
	for (unsigned int i = 0; i < coneTerms.size(); ++i)
		if ((divideByZero = coneTerms[i].computeDotProducts(
				currentPerturbation, l)))
			break;
	return divideByZero;
}//tryCurrentPerturbation


/**
 * For every vertex and ray, compute <v,l> and <r,l> and save the results.
 * If we divide by zero, return true. However, do not break early.
 * That is, if we divide by zero, we keep computing the rest of the dot products because we need
 * this information in any case.
 *
 * @l: linear form.
 *
 */
bool LinearPerturbationContainer::tryNoPerturbation(const vec_ZZ &l)
{
	divideByZero = false;
	for (unsigned int i = 0; i < coneTerms.size(); ++i)
		if (coneTerms[i].computeDotProducts(l))
			divideByZero = true; //don't break b/c want to find all inner products w/o pertrubation.
	return divideByZero;
}

/**
 * Finds a perturbation s.t we do not divide by zero.
 * First, we try no perturbation and hope for the best.
 * If we divide by zero, we set e to l + (0,0...,0,randon number,0,...,0), and take the orthogonal projection on l.
 * If this still does not work, we just take e to be e + (0,0...,0,randon number,0,...,0) and keep trying.
 *
 */
void LinearPerturbationContainer::findPerturbation(const vec_ZZ &l)
{

	ZZ gcdValue;
	divideByZero = false;
	int numberTimesDiviedByZero = 1;
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
		cout << "findPerturbation(): we divided by zero, trying new perturbation for the " << ++numberTimesDiviedByZero << "th time." << endl;
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

/**
 * Integrates the polytope over l.
 * Assumes we already computed the dot products with l
 * and we are no longer dividing by zero.
 *
 * Instead of saving the terms <r, -l-e>, we save <r,l+e>
 * It is here that we multiply the (-1)^d back in.
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

/** Finds the dot products of the vertex and rays for each term.
 *  If a residue calculation is needed, returns true, otherwise
 *  we did not divide by zero and false is returned.
 *
 *  We also find the abs. value of the determinant for this cone's rays.
 *
 *  For the <v,l>^M+d part, we only save <v,l> because we do not know M at this point.
 *  Instead of saving <r, l>, we factor the (-1)^d out and compute <r, l>.
 *
 *  We set the coeff. of the epsilon term and power term to zero to denote that the powers have not yet been processed.
 *
 *  @l: linear form.
 *  @return: true=we divided by zero wen <v,l>!= 0 (and so we need to find a perturbation for this term).
 *           false=we did not divide by zero or <v,l> = 0.
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

/** Finds the dot products of the vertex and rays for each term with the perturbation term ONLY.
*
* @e: perturbation vector: (a1, a1, ..., an)*epsilon
* @l: linear form. Oh, it is not currently used.
* @return: true:with this perturbation, we are still dividing by zero.
*          false: this perturbation worked or we do not need a perturbation.
*
* We already found <v, l> and <r, l>, if the zero perturbation does not make us divide by zero, we simply return false
* Otherwise, we find <v, e> and <r,e>
*
*
* We still ignore the powers of each term. The powers are all still zero.
*
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
 *
 * @totalSum: output parameter
 * @m: the power of the linear form
 * @dim: dimension of the polytope.
 *
 * It is here that we merge/update the powers in the denominator.
 *
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

}//integrateTerm


//prints out the current term.
//After updatePowers(), the term of the zero constants has a special form.
//Instead of (0 + ce)^m, the zero constants should be thought as 0 + c (e^m).
//This function is really just for debugging.
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
 *  Assumes a residue really need to be called. That is, assumes a (0+e) term exists.
 *  After this, if power = 0: there should be no 0 powers!
 *  					   k: then the power of this term is k.
 *  					  -1: There were repeated terms and this term has been merged with another term.
 *                            Thus, -1 is a lazy deletion.
 *  We swap the location of the (0+ce) term to be at location [0].
 *  example: (c1x + c0)*  (b1x + b0)*  (c1x + c0)    (0+ 2e) *(0+3e)
 *          =(c1x + c0)^2*(b1x + b0)^1*(c1x + c0)^-1*(0+ 6e)^2 *(0+3e)^-1 (before the swap)
 *          =(0 + 6e)^2*  (b1x + b0)^1*(c1x + c0)^-1*(c1x + c0)^2 *(0+3e)^-1 (after the swap)
 *              \_> note that the semantics of our data structure is (constant + epsilon*e)^power.
 *                  But for the zero constant terms, the semantic should be 0 + epsilon* (e^power).
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


