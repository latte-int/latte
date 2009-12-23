#define COEFF_MAX 10000

#include "PolyRep.h"

#include <iostream>
#include <fstream>
#include <NTL/ZZ.h>

using namespace std;

int main(int argc, char *argv[])
{
	if (argc < 2) { cout << "Usage: ./integrate filename" << endl; return 1; }
	string line;
	polynomial myPoly;
	ifstream myStream (argv[1]);
	if (!myStream.is_open()) { cout << "Error opening file " << argv[1] << ", please make sure it is spelled correctly." << endl; return 1; }
	while (!myStream.eof())
	{
		myStream >> line;
		loadPolynomial(myPoly, line);
		cout << "The following polynomial has " << myPoly.termCount << " terms containing " << myPoly.varCount << " variables: " << endl;
		cout << printPolynomial(myPoly) << endl;
		destroyPolynomial(myPoly);
	}
	myStream.close();
	/*
	//now, some testing
	srand ( time(NULL) );
	SetSeed( to_ZZ(rand()) );
	ZZ degree = RandomBnd(300);
	int numVars = rand() % 50;
	int termCount = rand() % 1000;
	
	cout << "Creating random polynomial of degree " << degree << " with " << numVars << " variables." << endl;
	stringstream output (stringstream::in | stringstream::out);
	output << "[";
	do
	{
		for (int i = 0; i < BLOCK_SIZE && termCount > -1; i++)
		{
			output << "[" << RandomBnd(COEFF_MAX) << ",[";
			for (int j = 0; j < myPoly.varCount; j++)
			{
				output << expTmp->data[i][j];
				if (j + 1 < myPoly.varCount)
				{ output << ","; }
			}
			output << "]]";
			if (termCount != 0)
			{ output << ","; }
			termCount--;
		}
		coeffTmp = coeffTmp->next; expTmp = expTmp->next;
	}
	while (coeffTmp != NULL);
	output << "]";
	loadPolynomial(myPoly, output.str());
	*/
	return 0; 
}
