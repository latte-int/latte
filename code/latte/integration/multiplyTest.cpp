#include "multiply.h"
#include <iostream>
#include <fstream>
#include <NTL/ZZ.h>

using namespace std;

int main(int argc, char *argv[])
{
	if (argc < 3) { cout << "Usage: " << argv[0] << " fileIn fileOut" << endl; return 1; }
	bool first = true;
	string line;
	monomialSum firstPoly, secondPoly, product;
	ifstream myStream (argv[1]);
	ofstream outStream(argv[2]);
	if (!myStream.is_open()) { cout << "Error opening file " << argv[1] << ", please make sure it is spelled correctly." << endl; return 1; }
	int polyCount = 0;
	string testForms;
	int *low, *high;
	
	MonomialLoadConsumer<ZZ>* myLoader = new MonomialLoadConsumer<ZZ>();	
	
	while (!myStream.eof())
	{
		getline(myStream, line, '\n');
		if (!line.empty())
		{
			if (first) //first polynomial
			{
				firstPoly.termCount = 0;
				myLoader->setMonomialSum(firstPoly);
				parseMonomials(myLoader, line);

				if (firstPoly.termCount == 0 || firstPoly.varCount == 0)
				{
					cout << "Error: loaded invalid monomial sum." << endl;
					return 1;
				}
				first = false;
			}
			else //second polynomial
			{
				secondPoly.termCount = 0;
				myLoader->setMonomialSum(secondPoly);
				parseMonomials(myLoader, line);

				if (secondPoly.termCount == 0 || secondPoly.varCount == 0)
				{
					cout << "Error: loaded invalid monomial sum." << endl;
					return 1;
				}
				first = true;
				
				low = new int[secondPoly.varCount];
				high = new int[secondPoly.varCount];
				for (int i = 0; i < secondPoly.varCount; i++)
				{
					low[i] = 0;
					high[i] = 100;
				}
				product.varCount = secondPoly.varCount;
				multiply<ZZ>(firstPoly, secondPoly, product, low, high);
				delete[] low;
				delete[] high;

				outStream << printMonomials(product) << endl;

				destroyMonomials(firstPoly);
				destroyMonomials(secondPoly);
				destroyMonomials(product);
			}
		}
	}

	delete myLoader;
	myStream.close();
	outStream.close();
	return 0; 
}
