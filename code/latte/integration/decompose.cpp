#include "PolyRep.h"

#include <iostream>
#include <fstream>

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
	return 0; 
}
