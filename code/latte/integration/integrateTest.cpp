#include "PolyRep.h"
#include "newIntegration.h"
#include "../timing.h"
#include <iostream>
#include <fstream>
#include <NTL/ZZ.h>

using namespace std;

NTL_CLIENT

int main(int argc, char *argv[])
{
	if (argc < 3) { cout << "Usage: " << argv[0] << " fileIn fileOut [decompose]" << endl; return 1; }
	bool decomposing = true; //decomposing by default
	bool polynomial = true; //file is assumed to alternate between polynomial and the simplex
	if (argc == 4) { decomposing = (strcmp(argv[3], "1") == 0); };
	string line;
	monomialSum monomials;
	linFormSum forms;
	ifstream myStream (argv[1]);
	ofstream outStream(argv[2]);
	if (!myStream.is_open()) { cout << "Error opening file " << argv[1] << ", please make sure it is spelled correctly." << endl; return 1; }
	int polyCount = 0;
	int dimension;
	ZZ degree = to_ZZ(-1);
	int irregularForms = 0;
	float loadTime, decomposeTime, integrateTime = 0.0f;
	Timer myTimer("Simplex Integration");
	while (!myStream.eof())
	{
		getline(myStream, line, '\n');
		if (!line.empty())
		{
			if (polynomial)
			{
				if (decomposing) //input is sum of monomials that we decompose into sum of linear forms
				{
					myTimer.start();
					loadMonomials(monomials, line);
					myTimer.stop();
					loadTime += myTimer.get_seconds();

					forms.termCount = 0;
					dimension = forms.varCount = monomials.varCount;
		
					cout << "Decomposing " << printMonomials(monomials);
					for (int i = 0; i < monomials.termCount; i++)
					{
						cout << ".";
						myTimer.start();
						decompose(monomials, forms, i);
						myTimer.stop();
						decomposeTime += myTimer.get_seconds();
					}
					cout << endl;
					outStream << printLinForms(forms) << endl;
					if (degree == to_ZZ(-1))
					{
						degree = 0;
						for (int i = 0; i < monomials.varCount; i++)
						{
							degree += monomials.eHead->data[0][i];
						}
					}
					destroyMonomials(monomials);
				}
				else //input is just linear forms
				{
					myTimer.start();
					loadLinForms(forms, line);
					myTimer.stop();
					loadTime += myTimer.get_seconds();
				}
				polynomial = false;
				//cout << "Loaded into " << forms.termCount << " linear forms" << endl;
			}
			else
			{
				//integrate here
				ZZ numerator, denominator;
				myTimer.start();
				integrateFlatVector(numerator, denominator, forms, line);
				myTimer.stop();
				if (IsZero(denominator)) //irregular
				{	
					irregularForms++;
				}
				integrateTime += myTimer.get_seconds();
				outStream << "[" << numerator << "," << denominator << "]" << endl;
				//cout << "Integral over " << line << " is " << numerator << "/" << denominator << endl;
				destroyLinForms(forms);
				polyCount++;
				polynomial = true;
			}
		}
	}
	if (decomposing) { cout << "Dimension " << dimension << ", degree " << degree << ". " << irregularForms << " forms were irregular." << endl; }
	cout << "Total time to load " << polyCount << " polynomials: " << loadTime << ", avg. is " << loadTime / polyCount << endl;
	if (decomposing) { cout << "Total time to decompose " << polyCount << " polynomials: " << decomposeTime << ", avg. is " << decomposeTime / polyCount << endl; }
	cout << "Total time to integrate " << polyCount << " polynomials: " << integrateTime << ", avg. is " << integrateTime / polyCount << endl;
	cout << "Total time is " << (decomposing ? integrateTime + decomposeTime : integrateTime) << ", avg. is " << (decomposing ? integrateTime + decomposeTime : integrateTime) / polyCount << endl;

	myStream.close();
	outStream.close();
	return 0; 
}
