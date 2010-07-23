#define COEFF_MAX 10000

#include "PolyRep.h"
#include "newIntegration.h"
#include <iostream>
#include <fstream>
#include <NTL/ZZ.h>
#include <cstring>

using namespace std;

int main(int argc, char *argv[])
{
	if (argc < 4) {cout<<"Usage: "<< argv[0] <<" fileIn fileOut IntegrateMode (polynomial/linear)" <<endl; return 1;}; //usage of the function
	_linFormSum forms;
	_monomialSum myPoly;
	simplexZZ mySimplex;
	ifstream myStream (argv[1]);
	ofstream hisStream (argv[2]);
	if (!myStream.is_open()) {cout<<"error opening file "<< argv[1] <<",please make sure it is spelled correctly!"<<endl; return 1; };
	ZZ a,b;
	string line,line2;
	int polyCount=0;
	while (!myStream.eof())
	{
		getline(myStream,line);delSpace(line);if (line.length()==0) break;
		getline(myStream,line2);delSpace(line2);
		polyCount++;
		if (!strcmp(argv[3],"polynomial")) 
		{
			a=0;b=0;
			_loadMonomials(myPoly, line);
			convertToSimplex(mySimplex, line2);
			_integrateMonomialSum(a, b, myPoly, mySimplex);
			hisStream<<a<<endl<<b<<endl;
		}
		else if (!strcmp(argv[3],"linear"))
		{
			a=0;b=0;
			_loadLinForms(forms, line);
			convertToSimplex(mySimplex, line2);
			LBlockIterator<RationalNTL>* it_ = new LBlockIterator<RationalNTL>();
			it_->setLists(forms.lHead, forms.cHead, forms.varCount, forms.termCount);
			integrateLinFormSum(a, b, it_, mySimplex);
			hisStream<<a<<endl<<b<<endl;
		};
	};
	myStream.close();
	hisStream.close();
	return 0;
}
