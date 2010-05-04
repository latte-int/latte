#define COEFF_MAX 10000

#include "PolyRep.h"

#include <iostream>
#include <fstream>
#include <NTL/ZZ.h>

using namespace std;

int main(int argc, char *argv[])
{
	if (argc < 3) { cout << "Usage: ./integrate fileIn fileOut" << endl; return 1; }
	string line;
	monomialSum myPoly;
	linFormSum lForm;
	ifstream myStream (argv[1]);
	ofstream outStream(argv[2]);
	float startTime, loadTime, decomposeTime;
	loadTime = decomposeTime = 0.0f;
	int count = 0;
	if (!myStream.is_open()) { cout << "Error opening file " << argv[1] << ", please make sure it is spelled correctly." << endl; return 1; }
	while (!myStream.eof())
	{
		getline(myStream, line, '\n');
		if (!line.empty())
		{
		startTime = time(NULL);
		loadMonomials(myPoly, line);
		loadTime += (time(NULL) - startTime);
		
		cout << "The following polynomial has " << myPoly.termCount << " terms containing " << myPoly.varCount << " variables: " << endl;
		cout << printPolynomial(myPoly) << endl;
		
		lForm.termCount = 0;
		lForm.varCount = myPoly.varCount;
		
		cout << "Decomposing";
		for (int i = 0; i < myPoly.termCount; i++)
		{
			cout << ".";
			startTime = time(NULL);
			decompose(myPoly, lForm, i);
			decomposeTime += (time(NULL) - startTime);
		}
		cout << endl;
		//cout << "About to print linear form to file" << endl;
		outStream << printLinForms(lForm) << endl; //print to output file
		
		//cout << "Maple expression is: " << endl;
		//outStream << printMapleForm(lForm) << endl; //let's print this instead
		destroyLinForms(lForm);
		destroyMonomials(myPoly);
		count++;
		}
	}
	cout << "Total time to load " << count << " polynomials: " << loadTime << ", avg is " << loadTime / count << endl;
	cout << "Total time to decompose " << count << " polynomials: " << decomposeTime << ", avg is " << decomposeTime / count << endl;

	myStream.close();
	outStream.close();
	return 0; 
}
