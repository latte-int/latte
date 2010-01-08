#define COEFF_MAX 10000

#include "PolyRep.h"
#include "newIntegration.h"

#include <iostream>
#include <fstream>
#include <NTL/ZZ.h>

using namespace std;

int main(int argc, char *argv[])
{
	if (argc < 2) { cout << "Usage: ./integrate fileIn" << endl; return 1; }
	string line,line2;
	polynomial myPoly;
	linearPoly lForm;
	ifstream myStream (argv[1]);
//	ofstream outStream(argv[2]);
	if (!myStream.is_open()) { cout << "Error opening file " << argv[1] << ", please make sure it is spelled correctly." << endl; return 1; };
		myStream >> line;
		loadPolynomial(myPoly, line);
		lForm.termCount = 0;
		lForm.varCount = myPoly.varCount;
		
		cout << "Decomposing";
		for (int i = 0; i < myPoly.termCount; i++)
		{
			cout << ".";
			decompose(myPoly, lForm, i);
		}
		cout << endl;
 		cout << "The polynomial is:" <<line<<endl;
		cout << "Integrating by decomposition" << endl;
		myStream >> line2;	//line2 is the simplex for example, [[0,0],[1,2],[2,3]]
		print_Integrate(myPoly.varCount,lForm,line2);
		destroyForm(lForm);
		destroyPolynomial(myPoly);
	myStream.close();
//	outStream.close();
	return 0; 
}
