/*
 * Perturbation.cpp
 *
 *  Created on: Nov 4
 *      Author: Brandon Dutra and Gregory Pinto
 */

#include "Perturbation.h"
#include "print.h"

using namespace std;

void LinearPerturbationContainer::setListCones(int dim,
		listCone * simpleConeList)
{

	coneTerms.resize(lengthListCone(simpleConeList));
	currentPerturbation.SetLength(dim);

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
		if ((divideByZero = coneTerms[i].computeDotProducts(currentPerturbation, l)))
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
	if (tryNoPerturbation(l) == false)
		return; //no errors. doing nothing worked!

	currentPerturbation = l;
	currentPerturbation[rand() % currentPerturbation.length()] += rand() % 10;

	//project p onto l: <p,l>/||l|| * l but we do not care about the scale, so don't divide.
	currentPerturbation = currentPerturbation - (currentPerturbation * l) * l;

	//see if we can scale the currentPerturbation down
	gcdValue = currentPerturbation[0];
	for (long i = 1; i < currentPerturbation.length(); ++i)
	{
		gcdValue = GCD(gcdValue, currentPerturbation[i]);
	}
	if (gcdValue > 1)
		for (long i = 1; i < currentPerturbation.length(); ++i)
			currentPerturbation[i] = currentPerturbation[i]/gcdValue;

	while (tryCurrentPerturbation(l) == true)
	{
		//we divided by zero again. that is <l+e, r>=0 for some ray and for for some cone while <l,v>!=0. :(

		//try a new perturbation.
		currentPerturbation[rand() % currentPerturbation.length()] += rand() % 10;

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

	//currentPerturbation is now normal to l.


}//findPerturbation

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
	for(listVector * ray = simplicialCone->rays; ray; ray = ray->rest, ++i)
	{
		//cout << "i=" << i << endl;

		rayDotProducts[i].constant = l * (ray->first);
		rayDotProducts[i].epsilon = 0;
		rayDotProducts[i].power = 0;


		if ( rayDotProducts[i].constant == 0)
			divideByZero = true;

	}//for each ray.
	//cout << "going to return " << divideByZero << endl;
	return divideByZero;
}


/** Finds the innerproducts of the vertex and rays for each term with the perturbation term ONLY.
 *  If a residue calculation is needed, returns true, otherwise
 *  we did not divide by zero and false is returned.
 */
bool LinearLawrenceIntegration::computeDotProducts(const vec_ZZ &e, const vec_ZZ &l)
{
	if (divideByZero == false)
		return false; //we do not need to do further process. plugging in numbers will work.

	//find that rays that dotted to zero and see if this new perturbation will work.
	unsigned int i = 0;
	for(listVector * ray = simplicialCone->rays; ray; ray = ray->rest, ++i)
	{
		rayDotProducts[i].epsilon = e * ray->first;
		if ( rayDotProducts[i].constant == 0 && rayDotProducts[i].epsilon == 0)
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

