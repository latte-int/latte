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
	
	BTrieIterator<ZZ, int>* it = new BTrieIterator<ZZ, int>();
	BTrieIterator<ZZ, int>* it2 = new BTrieIterator<ZZ, int>();
	
	BlockIterator<ZZ, int>* _it = new BlockIterator<ZZ, int>();
	BlockIterator<ZZ, int>* _it2 = new BlockIterator<ZZ, int>();
	
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
			//cout << myTime << " vs. ";
			
			_it->setLists(_myPoly.eHead, _myPoly.cHead, _myPoly.varCount, _myPoly.termCount);
			_it2->setLists(_myPoly.eHead, _myPoly.cHead, _myPoly.varCount, _myPoly.termCount);
			myTime = myTimer.get_seconds();
                        myTimer.start();
                        multiply<ZZ>(_it, _it2, _myProduct, low, high);
                        myTimer.stop();
                        myTime = myTimer.get_seconds() - myTime;
			//cout << myTime << endl;
                        //oldMult += myTime;
                        //cout << "Old Algorithm @ dimension " << _myProduct.varCount << ":" <<  myTime << "s. " << endl;
                        _destroyMonomials(_myProduct);
                                    
				    
			it->setTrie(myPoly.myMonomials, myPoly.varCount);
			it2->setTrie(myPoly.myMonomials, myPoly.varCount);
                        myTime = myTimer.get_seconds();
                        myTimer.start();
                        multiply<ZZ>(it, it2, myProduct, low, high);
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

                        for (int i = 0; i < _myPoly.termCount; i++)
                        {
                            _decompose(_myPoly, _myForms, i);
                        }
                        myTimer.stop();
                        myTime = myTimer.get_seconds() - myTime;
                        oldDecomp += myTime;
                        _destroyLinForms(_myForms);                        
                        
                        myForms.termCount = 0;
			it->setTrie(myPoly.myMonomials, myPoly.varCount);
                        myTime = myTimer.get_seconds();
                        myTimer.start();
                        decompose(it, myForms);
                        myTimer.stop();
                        myTime = myTimer.get_seconds() - myTime;
                        newDecomp += myTime;
                        destroyLinForms(myForms);
                    }
		    //destroyMonomials(myPoly);
                    //_destroyMonomials(_myPoly);
		}
	}

	delete myLoader;
	delete _myLoader;
	delete it; delete it2;
	delete _it; delete _it2;
	myStream.close();
       cout << "Comparing @ dimension " << dimension << " (old vs. new) " << endl;
       if (option == 1) { cout << "Squaring polynomials: " << oldMult << " vs. " << newMult << endl; }
       else if (option == 2) { cout << "Decomposing: " << oldDecomp << " vs. " << newDecomp << endl; }
	return 0; 
}
