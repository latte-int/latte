/*
Provides interactive console tool for decomposition and/or integration over arbitrary simplices
*/

#include "PolyTrie.h"
#include "newIntegration.h"
#include "../timing.h"
#include <iostream>
#include <fstream>
#include <NTL/ZZ.h>
#include <string>

using namespace std;

int main(int argc, char *argv[])
{
	bool options = true;
	bool decomposing = false;
	ifstream inFile; istream inStream (std::cin.rdbuf());
	ofstream outFile; ostream outStream (std::cout.rdbuf());
	for (int i = 1; i < argc; i++)
	{
		if (options)
		{
			if (strncmp(argv[i], "-", 1) != 0)
			{
				options = false;
				inFile.open(argv[i], ios::in);
				if (!inFile.is_open()) { cout << "Error: cannot open " << argv[i] << ", please try again." << endl; return 1; }
				inStream.rdbuf(inFile.rdbuf());
			}
			else
			{
				if (strcmp(argv[i], "-m") == 0 || strcmp(argv[i], "-monomial"))
				{
					decomposing = true;
				}
			}
		}
		else
		{
			outFile.open(argv[i], ios::out);
			if (!outFile.is_open()) { cout << "Error: cannot open " << argv[i] << ", please try again." << endl; return 1; }
			outStream.rdbuf(outFile.rdbuf());
		}
	}
	
	bool polynomial = true; //file is assumed to alternate between monomial / linform sum and the associated simplex
	float lastTime;
	float sampleTime;
	string line;
	monomialSum monomials;
	linFormSum forms;
	int polyCount = 0;
	int irregularForms = 0;
	float loadTime, decomposeTime, integrateTime;
	loadTime = decomposeTime = integrateTime = 0.0f;
	Timer myTimer("Integration timer");
	FormIntegrateConsumer<ZZ> *integrator;
	if (!decomposing) { integrator = new FormIntegrateConsumer<ZZ>(); }

	BTrieIterator<ZZ, int>* it = new BTrieIterator<ZZ, int>();
	while (!inStream.eof())
	{
		getline(inStream, line, '\n');
		if (!line.empty())
		{
			if (polynomial) //reading polynomial
			{
				sampleTime = myTimer.get_seconds();
				if (decomposing) //input is sum of monomials that we decompose into sum of linear forms
				{
					cout << "Loading monomials..." << endl;
					if (options) { outStream << line << endl; }
					lastTime = myTimer.get_seconds();
					myTimer.start();
					loadMonomials(monomials, line);
					myTimer.stop();
					loadTime += (myTimer.get_seconds() - lastTime);
					outStream << "Loaded monomials in " << myTimer.get_seconds() - lastTime << "s." << endl;

					cout << printMonomials(monomials) << endl;
	
					if (monomials.termCount == 0 || monomials.varCount == 0)
					{
						cout << "Error: loaded invalid monomial sum." << endl;
						return 1;
					}

					forms.termCount = 0;
					forms.varCount = monomials.varCount;
		
					cout << "Decomposing into sum of linear forms..." << endl;
					
					it->setTrie(monomials.myMonomials, monomials.varCount);
					
					lastTime = myTimer.get_seconds();
					myTimer.start();
					decompose(it, forms);
					myTimer.stop();
					decomposeTime += (myTimer.get_seconds() - lastTime);
					outStream << "Decomposition finished in " << myTimer.get_seconds() - lastTime << "s." << endl;
					
					if (forms.termCount == 0 || forms.varCount == 0)
					{
						cout << "Error: no terms in decomposition to sum of linear forms.";
						return 1;	
					}
					destroyMonomials(monomials);
				}
				else //input is just linear forms
				{
					cout << "Reading sum of powers of linear forms..." << endl;
					if (options) { outStream << line << endl; }
					
					/*lastTime = myTimer.get_seconds();
					myTimer.start();
					loadLinForms(forms, line);
					myTimer.stop();
					loadTime += (myTimer.get_seconds() - lastTime);
					outStream << "Loaded forms in " << myTimer.get_seconds() - lastTime << "s." << endl;

					cout << printLinForms(forms) << endl;
	
					if (forms.varCount == 0)
					{
						cout << "Error: loaded invalid form sum." << endl;
						return 1;
					}
					destroyLinForms(forms);*/
					integrator->setFormSum(line);
				}
				polynomial = false;
			}
			else //reading simplex
			{
				simplexZZ mySimplex;
				cout << "Reading simplex..." << endl;
				convertToSimplex(mySimplex, line);
				//integrate here
				ZZ numerator, denominator;
				cout << "Integrating..." << endl;
				if (decomposing)
				{
					lastTime = myTimer.get_seconds();
					myTimer.start();
					integrateLinFormSum(numerator, denominator, forms, mySimplex);
					myTimer.stop();
					integrateTime += (myTimer.get_seconds() - lastTime);
					destroyLinForms(forms);
				}
				else
				{
					integrator->setSimplex(mySimplex);
					lastTime = myTimer.get_seconds();
					myTimer.start();
					parseLinForms(integrator, integrator->getFormSum());
					myTimer.stop();
					integrator->getResults(numerator, denominator);
					integrateTime += (myTimer.get_seconds() - lastTime);
				}
				outStream << "Integrated in " << myTimer.get_seconds() - lastTime << "s." << endl;
				outStream << "Calculated integral is \\frac{" << numerator << "}{" << denominator << "}" << endl;
				polyCount++;
				polynomial = true;
			}
		}
	}
	delete it;
	
	outStream << "Processed " << polyCount;
	if (decomposing)
	{
		outStream << " monomial sums." << endl;
		outStream << "Average load time: " << (loadTime)/(polyCount + 0.0) << endl;
		outStream << "Average decomposition time:" << (decomposeTime)/(polyCount + 0.0) << endl;
	}
	else
	{
		outStream << " sums of powers of linear forms." << endl;
	}
	outStream << "Average integration time:" << (integrateTime)/(polyCount + 0.0) << endl;
	
	if (!decomposing) { delete integrator; }
	if (inFile.is_open()) { inFile.close(); }
	if (outFile.is_open()) { outFile.close(); }

	return 0; 
}
