#include "PolyTrie.h"
#include "newIntegration.h"
#include "../timing.h"
#include <iostream>
#include <fstream>
#include <NTL/ZZ.h>

using namespace std;

void timedOut()
{

}

int main(int argc, char *argv[])
{
	if (argc < 2)
	{
		cout << "Usage: " << argv[0] << " [options] inFile outFile" << endl; return 1;
	}
	bool options = true;
	bool decomposing = false;
	bool benchmarking = false;
	char benchFile[128];
	float myTimeout = 0.0f;
	ifstream inFile;
	ofstream outFile;
	
	for (int i = 1; i < argc; i++)
	{
		if (options)
		{
			if (strncmp(argv[i], "-", 1) != 0)
			{
				options = false;
				if (argc < (i + 2)) //has to be one more file after this
				{
					cout << "Usage: " << argv[0] << " [options] inFile outFile" << endl; return 1;
				}
				inFile.open(argv[i], ios::in);
				if (!inFile.is_open()) { cout << "Error: cannot open " << argv[i] << ", please try again." << endl; return 1; }
			}
			else
			{
				if (strcmp(argv[i], "-m") == 0) //monomial
				{
					decomposing = true;
				}
				else if (strcmp(argv[i], "-b") == 0) //benchmark
				{
					i++;
					strncpy (benchFile, argv[i], 128);
					benchmarking = true;
				}
				else if (strcmp(argv[i], "-t") == 0) //timeout
				{
					i++;
					myTimeout = atof(argv[i]);
				}
			}
		}
		else
		{
			outFile.open(argv[i], ios::out);
			if (!outFile.is_open()) { cout << "Error: cannot open " << argv[i] << ", please try again." << endl; return 1; }
		}
	}
	if (options) { cout << "Usage: " << argv[0] << " [options] inFile outFile" << endl; return 1; }

	bool polynomial = true; //file is assumed to alternate between polynomial and its simplex

	float tempTime, iterationTime;
	string line;
	monomialSum monomials;
	linFormSum forms;

	int polyCount = 0;
	float loadTime, decomposeTime, integrateTime;
	loadTime = decomposeTime = integrateTime = 0.0f;
	Timer myTimer("Integration timer");
	FormIntegrateConsumer<ZZ> *integrator;
	if (!decomposing) { integrator = new FormIntegrateConsumer<ZZ>(); }

	while (!inFile.eof())
	{
		getline(inFile, line, '\n');
		if (!line.empty())
		{
			
			if (polynomial) //reading form
			{
				iterationTime = myTimer.get_seconds();
				if (decomposing) //input is sum of monomials that we decompose into sum of linear forms
				{
					tempTime = myTimer.get_seconds();
					myTimer.start();
					loadMonomials(monomials, line);
					myTimer.stop();
					loadTime += (myTimer.get_seconds() - tempTime);

					if (monomials.termCount == 0 || monomials.varCount == 0)
					{
						cout << "Error: loaded invalid monomial sum." << endl;
						return 1;
					}

					forms.termCount = 0;
					forms.varCount = monomials.varCount;
		
					tempTime = myTimer.get_seconds();
					myTimer.start();
					//for (int i = 0; i < monomials.termCount; i++)
					//{
						//cout << ".";
					//	decompose(monomials, forms, i);
					//}
					decompose(monomials, forms);
					myTimer.stop();
					decomposeTime += (myTimer.get_seconds() - tempTime);
					//cout << endl;
					
					if (forms.termCount == 0)
					{
						cout << "Error: no terms in decomposition to sum of linear forms.";
						return 1;	
					}
					
					if (!benchmarking)
					{ outFile << printLinForms(forms) << endl; }

					destroyMonomials(monomials);
				}
				else //input is just linear forms
				{
					integrator->setFormSum(line);
				}
				polynomial = false;
			}
			else //reading simplex
			{
				simplexZZ mySimplex;
				//cout << "reading simplex" << endl;
				//if (decomposing) { destroyLinForms(forms); }
				convertToSimplex(mySimplex, line);
				ZZ numerator, denominator;
				if (decomposing)
				{
					tempTime = myTimer.get_seconds();
					myTimer.start();
					integrateLinFormSum(numerator, denominator, forms, mySimplex);
					myTimer.stop();
					integrateTime += (myTimer.get_seconds() - tempTime);
					destroyLinForms(forms);
				}
				else
				{
					integrator->setSimplex(mySimplex);
					tempTime = myTimer.get_seconds();
					myTimer.start();
					parseLinForms(integrator, integrator->getFormSum());
					myTimer.stop();
					integrator->getResults(numerator, denominator);
					integrateTime += (myTimer.get_seconds() - tempTime);
				}
				outFile << "[" << numerator << "," << denominator << "]" << endl;
				polyCount++;
				polynomial = true;
				if (benchmarking) { cout << "Sample took " << myTimer.get_seconds() - iterationTime << "s." << endl; }
				if (myTimeout > 0.001 && (myTimer.get_seconds() - iterationTime) > myTimeout) //we timed out
				{
					delete integrator;
					inFile.close();
					
					outFile.seekp(0, ios_base::beg);
					outFile.write("Error\n", 6);
					outFile.close();
					
					//FILE* myFile = fopen(,"w"); //overwriting results file
					//fprintf(myFile, "Error");
					//fclose(myFile);
											
					cout << "Integration timed out." << endl;
					if (benchmarking)
					{
						FILE* myFile = fopen(strtok(benchFile, " "),"a");
						fprintf(myFile, "%10s", "--");
						fclose(myFile);
					}
					return 1;
				}
			}
		}
	}

	if (benchmarking)
	{
		cout << "Total time " << (decomposing ? loadTime + integrateTime + decomposeTime : integrateTime) << endl;
		cout << "       avg " << (decomposing ? loadTime + integrateTime + decomposeTime : integrateTime) / (polyCount + 0.0) << endl;
		FILE* benchmarks = fopen(strtok(benchFile, " "),"a");
		fprintf(benchmarks, "%10.2f", (decomposing ? loadTime + integrateTime + decomposeTime : integrateTime) / polyCount);
		fclose(benchmarks);
	}
	
	if (!decomposing) { delete integrator; }
	inFile.close();
	outFile.close();
	return 0; 
}
