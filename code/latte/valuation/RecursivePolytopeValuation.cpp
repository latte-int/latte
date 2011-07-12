/*
 * RecursivePolytopeValuation.cpp
 *
 *  Created on: May 22, 2011
 *      Author: bedutra
 */


#include "RecursivePolytopeValuation.h"
#include "sstream"

RecursivePolytopeValuation::RecursivePolytopeValuation():maxRecursiveLevel(3), recursiveLevel(0), minAffineDimension(0)
{
	// TODO Auto-generated constructor stub

}

RecursivePolytopeValuation::~RecursivePolytopeValuation() {
	// TODO Auto-generated destructor stub
}

RationalNTL RecursivePolytopeValuation::findVolume(ReadPolyhedronDataRecursive & readPolyhedron, BarvinokParameters * parm)
{

	RationalNTL ans;

	//set up the linear forms.
	int degree = 0;
	vec_ZZ exp;
	exp.SetLength(parm->Number_of_Variables);
	for(int i = 0; i < parm->Number_of_Variables; ++i)
		exp[i] = i+1;

#if 0
	exp[0]=1;
	//exp[1]=0;
	//exp[2]=0;

	insertLinForm(RationalNTL(1,1), degree, exp, linform);

	readPolyhedron.latticeInverse();

	PolytopeValuation pv(readPolyhedron.getPolyhedron(), *parm);
	pv.setLatticeInverse(readPolyhedron.getLatticeInverse(), readPolyhedron.getLatticeInverseDilation());
	pv.setFullDimension(readPolyhedron.getFullDimensionCount());
	cout << " getFull dim = " << readPolyhedron.getFullDimensionCount() << endl;
	ans = pv.findIntegral(linform);


	destroyLinForms(linform);
#endif


	ans = findVolume_recursive(readPolyhedron, parm, 0, exp);
	cout << "volume = " << ans << endl;
	cout << "volume = " << readPolyhedron.volumeCorrection(ans) << endl;

	exit(1);
	//+++++++++++++++++++++++++++++

}//findVolume


RationalNTL RecursivePolytopeValuation::findVolume_recursive(ReadPolyhedronDataRecursive & readPolyhedron, BarvinokParameters * parm, int power, vec_ZZ & l)
{
	++recursiveLevel;
	cout << "******starting find vol rec: " << recursiveLevel << endl;
	RationalNTL ans;
	stringstream out;

	if(recursiveLevel >= maxRecursiveLevel || readPolyhedron.getFullDimensionCount() <= minAffineDimension )
	{
		//set up the linear forms.
		linFormSum linform;
		linform.termCount = 0;
		linform.varCount = parm->Number_of_Variables;

		cout << "Integrating form " << l << " to the power " << power << endl;
		ZZ powerFactorial;
		powerFactorial = PolytopeValuation::factorial(power);
		insertLinForm(RationalNTL(powerFactorial,1), power, l, linform);
		Polyhedron *Poly = readPolyhedron.findTangentCones();
		readPolyhedron.latticeInverse();

		PolytopeValuation pv(Poly, *parm);
		pv.setLatticeInverse(readPolyhedron.getLatticeInverse(), readPolyhedron.getLatticeInverseDilation());
		pv.setFullDimension(readPolyhedron.getFullDimensionCount());
		ans = pv.findIntegral(linform);


		destroyLinForms(linform);
		delete Poly;

		--recursiveLevel;

		//cout << "base case done, ans=" << ans << " num of var=" << parm->Number_of_Variables << endl;
		//exit(1);
		return ans;
	}//base case. don't recurse. Find the volume of this polytope.

	//now we need to use stokes.
	++power;

	int numRows = readPolyhedron.getNumberRows();
	out << "(";
	for(int i = 1; i <= numRows; ++i)
	{
		RationalNTL newSum;
		RationalNTL lDotNormal;

		ReadPolyhedronDataRecursive newMatrix(readPolyhedron);

		readPolyhedron.getFacetPolytope(i, newMatrix, l, lDotNormal);

		//readPolyhedron.findTangentCones();
		//readPolyhedron.latticeInverse();

		if ( lDotNormal.getNumerator() != 0)
		{
cout << "start here: val/recPolyVal.cpp:120 " << endl;
exit(1); //error ntl to zz error.
			RationalNTL normalCorrection = newMatrix.getNormalFactor();
			cout << "go to find vol recur " << i << "/" << numRows << endl;
			newSum = findVolume_recursive(newMatrix, parm, power, l);
			cout << "***Level recursiveLevel " << recursiveLevel << ", row " << i << "/" << numRows << ", int" << newSum << " ." << endl;
			ans.add(lDotNormal * newSum * normalCorrection);
			out << " + " << lDotNormal << "*" << newSum << "*" << normalCorrection;
		}
		else
			cout << "***Level recursiveLevel " << recursiveLevel << ", row " << i << "/" << numRows << ", dot is zero" << " ." << endl;


	}//for every facet
	out << ") / (" << power << "* " << l*l << ") ";

	ans.div(to_ZZ(power));
	ans.div(to_ZZ(l*l)); // <y,l> , dot product.
	cout << "******leaving find vol rec: " << recursiveLevel << " ans=" << ans << endl;
	cout << "out = " << out.str().c_str();
	--recursiveLevel;
	return ans;
}


void RecursivePolytopeValuation::setMaxRecursiveLevel(int i) { maxRecursiveLevel = i;}
void RecursivePolytopeValuation::setMinDimension(int i) {minAffineDimension = i;}
