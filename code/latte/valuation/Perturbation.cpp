/*
 * Perturbation.cpp
 *
 *  Created on: Nov 4
 *      Author: Brandon Dutra and Gregory Pinto
 */

#include "Perturbation.h"

using namespace std;



void LinearPerturbationContainer::setListCones(int dim, listCone * simpleConeList)
{

	coneTerms.resize(lengthListCone(simpleConeList));
	currentPerturbation.setLength(dim);

	//build the vector of cones.
	listCone* ptr = simpleConeList;
	for (i = 0; i < coneTerms.size(); ++i)
	{
		coneTerms[i].setSimplicialCone(ptr);
		ptr = ptr->rest;
	}//for every cone.
}


bool LinearPerturbationContainer::tryNoPerturbation()
{
	divideByZero = false;
	for (unsigned int i = 0; i < coneTerms.size(); ++i)
		if ( (divideByZero = coneTerms[i].computeDotProducts()))
				break;
	return divideByZero;
}


void LinearPerturbationContainer::findPerturbation(const vec_ZZ &l)
{

	if (tryNoPerturbation() == false )
		return; //no errors. doing nothing worked!

	currentPerturbation = l;
	currentPerturbation[rand()%currentPerturbation.size()] += rand()%10;

	currentPerturbation = currentPerturbation - (currentPerturbation*l)/(l*l)) * l;
	//currentPerturbation is now normal to l.



}//findPerturbation

///--------------------------------------

LinearLawrenceIntegration::LinearLawrenceIntegration():divideByZero(0), simplicialCone(NULL)
{

}

LinearLawrenceIntegration::LinearLawrenceIntegration(listCone * cone):divideByZero(0), simplicialCone(cone)
{

}


void LinearLawrenceIntegration::setSimplicialCone(listCone *cone)
{
	simplicialCone = cone;
}

bool LinearLawrenceIntegration::computeDotProducts(vec_ZZ e)
{

	return false;
}







