#include "multiply.h"
#include "PolyTrie.h"
#include "PolyRep.h"
#include "../timing.h"
#include <iostream>
#include <fstream>
#include <NTL/ZZ.h>

using namespace std;

int main(int argc, char *argv[])
{
	string line;
        int option = 0;
        if (strcmp(argv[1], "-m") == 0) { option = 1; } //multiply
        else if (strcmp(argv[1], "-d") == 0) { option = 2; } //decompose
        int dimension = -1;
	monomialSum myPoly, myProduct;
	_monomialSum _myPoly, _myProduct;
        linFormSum myForms;
        _linFormSum _myForms;
	ifstream myStream (argv[2]);
	if (!myStream.is_open()) { cout << "Error opening file " << argv[2] << ", please make sure it is spelled correctly." << endl; return 1; }
	int polyCount = 0;
	string testForms;
	int *low, *high;
	
	float myTime;
	Timer myTimer("");
        
        float oldMult, oldDecomp, newMult, newDecomp;
        oldMult = oldDecomp = newMult = newDecomp = 0.0f;
	
	MonomialLoadConsumer<ZZ>* myLoader = new MonomialLoadConsumer<ZZ>();
	_MonomialLoadConsumer<ZZ>* _myLoader = new _MonomialLoadConsumer<ZZ>();
	
	while (!myStream.eof())
	{
		getline(myStream, line, '\n');
		if (!line.empty())
		{
		    myPoly.termCount = 0;
		    myLoader->setMonomialSum(myPoly);
		    parseMonomials(myLoader, line);

		    if (myPoly.termCount == 0 || myPoly.varCount == 0)
		    {
			cout << "Error: loaded invalid monomial sum." << endl;
		    	return 1;
		    }
                    else { myForms.varCount = myProduct.varCount = myPoly.varCount; if (dimension == -1) { dimension = myPoly.varCount; } }
				
		    _myPoly.termCount = 0;
		    _myLoader->setMonomialSum(_myPoly);
		    _parseMonomials(_myLoader, line);

		    if (_myPoly.termCount == 0 || _myPoly.varCount == 0)
		    {
			cout << "Error: loaded invalid monomial sum." << endl;
			return 1;
		    }
                    else { _myForms.varCount = _myProduct.varCount = _myPoly.varCount; }
                    
                    //square the polynomials
                    if (option == 1)
                    {
                        low = new int[myProduct.varCount];
                        high = new int[myProduct.varCount];
                        for (int i = 0; i < myProduct.varCount; i++)
                        {
                            low[i] = INT_MIN;
                            high[i] = INT_MAX;
                        }
                        
                        myTime = myTimer.get_seconds();
                        myTimer.start();
                        _multiply<ZZ>(_myPoly, _myPoly, _myProduct, low, high);
                        myTimer.stop();
                        myTime = myTimer.get_seconds() - myTime;
                        oldMult += myTime;
                        //cout << "Old Algorithm @ dimension " << _myProduct.varCount << ":" <<  myTime << "s. " << endl;
                        _destroyMonomials(_myProduct);
                                    
                        myTime = myTimer.get_seconds();
                        myTimer.start();
                        multiply<ZZ>(myPoly, myPoly, myProduct, low, high);
                        myTimer.stop();
                        myTime = myTimer.get_seconds() - myTime;
                        newMult += myTime;
                        //cout << "New Algorithm @ dimension " << myProduct.varCount << ":" <<  myTime << "s. " << endl;
                        destroyMonomials(myProduct);

			   delete [] low;
			   delete [] high;
                    }
                    else if (option == 2)
                    {
                        _myForms.termCount = 0;
                        myTime = myTimer.get_seconds();
                        myTimer.start();
cout << "Decomposing" << endl;
                        for (int i = 0; i < _myPoly.termCount; i++)
                        {
cout << ". ";
                            _decompose(_myPoly, _myForms, i);
                        }
cout << endl;
                        myTimer.stop();
                        myTime = myTimer.get_seconds() - myTime;
                        oldDecomp += myTime;
                        _destroyLinForms(_myForms);                        
                        
                        myForms.termCount = 0;
                        myTime = myTimer.get_seconds();
                        myTimer.start();
                        decompose(myPoly, myForms);
                        myTimer.stop();
                        myTime = myTimer.get_seconds() - myTime;
                        newDecomp += myTime;
                        destroyLinForms(myForms);
                    }
		    destroyMonomials(myPoly);
                    _destroyMonomials(_myPoly);
		}
	}

	delete myLoader;
	delete _myLoader;
	myStream.close();
       cout << "Comparing @ dimension " << dimension << " (old vs. new) " << endl;
       if (option == 1) { cout << "Squaring polynomials: " << oldMult << " vs. " << newMult << endl; }
       else if (option == 2) { cout << "Decomposing: " << oldDecomp << " vs. " << newDecomp << endl; }
	return 0; 
}
