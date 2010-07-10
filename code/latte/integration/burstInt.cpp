#define COEFF_MAX 10000

#include "newIntegration.h"
#include "PolyTrie.h"
#include <iostream>
#include <fstream>
#include <NTL/ZZ.h>
#include <cstring>

using namespace std;

//This main function takes in a number of polynomial-simplex or liear-form-simplex pairs and integrate them.

int main(int argc, char *argv[])
{
	if (argc < 4) {cout<<"Usage: "<< argv[0] <<" fileIn fileOut IntegrateMode (polynomial/linear)" <<endl; return 1;}; //usage of the function
	linFormSum forms;
	monomialSum myPoly;
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
			loadMonomials(myPoly, line);
			convertToSimplex(mySimplex, line2);
			integrateMonomialSum(a, b, myPoly, mySimplex);
			hisStream<<a<<endl<<b<<endl;
		}
		else if (!strcmp(argv[3],"linear"))
		{
			a=0;b=0;
			loadLinForms(forms, line);
			BTrieIterator<ZZ, ZZ>* it = new BTrieIterator<ZZ, ZZ>();			
			it->setTrie(forms.myForms, forms.varCount);
			convertToSimplex(mySimplex, line2);
			integrateLinFormSum(a, b, it, mySimplex);
			hisStream<<a<<endl<<b<<endl;
		};
	};
	myStream.close();
	hisStream.close();
	return 0;
}
