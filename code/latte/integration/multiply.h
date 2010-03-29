#ifndef MULTIPLY_H
#define MULTIPLY_H
#include "PolyTrie.h"


class TermMultiplier: public TrieIterator<ZZ, int> {
    public:
    burstTrie<ZZ, int>* myParent;
    term<ZZ, int>* myTerm;
    int* min;
    int* max;
        
    void init(burstTrie<ZZ, int>* parent, term<ZZ, int>* term, int* myMin, int* myMax)
    {
	cout << "INIT: Term Multiplier" << endl;
	if (!parent) { cout << "NULL PARENT" << endl; }
	myParent = parent;
	myTerm = term;
	cout << "Init with term: [" << *myTerm->coef << ", [";
	for (int i = 0; i < myTerm->expCount; i++)
	{ cout << myTerm->exponents[i] << ", "; }
	cout << "]]" << endl;
	min = myMin;
	max = myMax;
    }
    
    int addExp(int a, int b, int i)
    {
	if ((a + b) > max[i]) { return max[i]; }
	if ((a + b) < min[i]) { return min[i]; }
	return a + b;
    }
    
    void consumeTerm(term<ZZ, int>* thisTerm) //insert thisTerm*myTerm into myParent
    {
	cout << "Multiplying terms..." << endl;
	cout << "Enumerated term: [" << *thisTerm->coef << ", [";
	for (int i = 0; i < thisTerm->expCount; i++)
	{ cout << thisTerm->exponents[i] << ", "; }
	cout << "]]" << endl;
	term<ZZ, int>* newTerm = (term<ZZ, int>*) malloc(sizeof(term<ZZ, int>));
	newTerm->coef = new ZZ((*thisTerm->coef) * (*myTerm->coef));
	cout << "New coef is " << *newTerm->coef << endl;
	newTerm->expCount = dimension;
	cout  << "exp count is " << newTerm->expCount << endl;
	newTerm->exponents = new int[dimension];
	cout << "Depth is " << curDepth << endl;
	int i = 0;
	for (; i < curDepth; i++) //sorted values
	{
		cout << "1. Adding " << curPrefix[i] << endl;
		newTerm->exponents[i] = addExp(myTerm->exponents[i], curPrefix[i], i);
	}
	for (; i < thisTerm->expCount + curDepth; i++) //stored exps
	{
		cout << "2. Adding " << thisTerm->exponents[i - curDepth] << endl;
		newTerm->exponents[i] = addExp(myTerm->exponents[i], thisTerm->exponents[i - curDepth], i);
	}
	for ( ; i < dimension; i++) //trailing zeroes
	{
		cout << "3. Adding 0" << endl;
		newTerm->exponents[i] = addExp(myTerm->exponents[i], 0, i);
	}
        newTerm->degree = -1;
	cout << "Finished, inserting.." << endl;
	
	trieInsert(myParent, newTerm);
	cout << "Finished" << endl;
    }
};

class TrieMultiplier: public TrieIterator<ZZ, int> {
    public:
    burstTrie<ZZ, int>* myTrie;
    burstTrie<ZZ, int>* result;
    TermMultiplier* myMultiplier;
    int* min;
    int* max;
     
    burstTrie<ZZ, int>* getProduct()
    {
	return result;
    } 
     
    void init(burstTrie<ZZ, int>* first, int* myMin, int* myMax)
    {
	cout << "INIT: Trie Multiplier" << endl;
        myTrie = first;
	
	cout << "Creating product trie" << endl;
	result = (burstTrie<ZZ, int>*) malloc (sizeof(burstTrie<ZZ, int>));
	cout << "Creating range" << endl;
	result->range = new int[2];
	cout << "Setting upper and lower bound to 0" << endl;
	result->range[0] = 0;
	result->range[1] = 0;
	cout << "Creating first container" << endl;
	result->firstCont = createContainer<ZZ, int>();
	
	myMultiplier = new TermMultiplier();
	myMultiplier->curPrefix = NULL;
	myMultiplier->setDimension(dimension);
	min = myMin;
	max = myMax;
    }
    
    void consumeTerm(term<ZZ, int>* myTerm) //for each term in the first monomial sum, multiply it by everything in the second
    {
	cout << "Multiplying trie by term .. " << endl;
	term<ZZ, int>* newTerm = (term<ZZ, int>*) malloc(sizeof(term<ZZ, int>));
	newTerm->coef = new ZZ(*myTerm->coef);
	cout << "New coef is " << *newTerm->coef << endl;
	newTerm->expCount = dimension;
	cout  << "exp count is " << newTerm->expCount << endl;
	newTerm->exponents = new int[dimension];
	int i = 0;
	for (; i < curDepth; i++) //sorted values
	{
		cout << "1. Setting to " << curPrefix[i] << endl;
		newTerm->exponents[i] = curPrefix[i];
	}
	for (; i < myTerm->expCount + curDepth; i++) //stored exps
	{
		cout << "2. Setting to " << myTerm->exponents[i - curDepth] << endl;
		newTerm->exponents[i] = myTerm->exponents[i - curDepth];
	}
	for (; i < dimension; i++) //trailing zeroes
	{
		newTerm->exponents[i] = 0;
	}
        newTerm->degree = -1;
	myMultiplier->init(result, newTerm, min, max);
	cout << "Created term .." << endl;
	myMultiplier->enumTrie(myTrie);
	cout << "Destroying.." << endl;
	destroyTerm(newTerm);
    }
    
    void finish()
    {
	delete myMultiplier;
	myMultiplier = NULL;
    }
};


// Multipies two monomial sums, storing the result in the third one
// Any values stored in result will be overwritten
// result is every term in the product of two monomial sums whose exponents are greater than min and lower than max.
// min, max point to int arrays of length result.varCount
template <class T>
void multiply(monomialSum& first, monomialSum& second, monomialSum& result, int* min, int* max)
{
	if (first.termCount == 0 || second.termCount == 0) { cout << "Only one monomial sum given, aborting."; return; }
	cout << "Dimensions are " << first.varCount << " x " << second.varCount << " = " << result.varCount << endl;
	cout << "Printing tree one" << endl;
	TriePrinter<T, int>* myPrinter = new TriePrinter<T, int>();
	myPrinter->setDimension(first.varCount);
	cout << myPrinter->printTrie(first.myMonomials) << endl;
	cout << "Printing two" << endl;
	cout << myPrinter->printTrie(second.myMonomials) << endl;
	delete myPrinter;
	cout << "Multiplying ..." << endl;
	TrieMultiplier* genie = new TrieMultiplier();
	genie->setDimension(result.varCount);
	genie->init(first.myMonomials, min, max);
	genie->enumTrie(second.myMonomials);
	cout << "Finished" << endl;
	result.myMonomials = genie->getProduct();
	
	genie->finish();
	delete genie;
	//cout << "Yo" << endl;
	/*eBlock* firstExp = first.eHead; cBlock<T>* firstCoef = first.cHead;
	eBlock* secondExp = second.eHead; cBlock<T>* secondCoef = second.cHead;
	
	int* exponents = new int[result.varCount];
	bool valid;
	result.termCount = 0;
	for (int i = 0; i < first.termCount; i++)
	{
		if (i > 0 && i % BLOCK_SIZE == 0) //this block is done, get next one
		{
			firstExp = firstExp->next; firstCoef = firstCoef->next;
		}
		for (int j = 0; j < second.termCount; j++)
		{
			if (j > 0 && j % BLOCK_SIZE == 0) //this block is done, get next one
			{
				secondExp = secondExp->next; secondCoef = secondCoef->next;
			}
			valid = true;
			for (int k = 0; k < result.varCount; k++)
			{
				exponents[k] = firstExp->data[(i % BLOCK_SIZE)*first.varCount + k] + secondExp->data[(j % BLOCK_SIZE)*second.varCount + k];
				if (exponents[k] < min[k] || exponents[k] > max[k]) {valid = false; break; }
			}
			if (valid) //all exponents are within range
			{
				insertMonomial<T>(firstCoef->data[i % BLOCK_SIZE] * secondCoef->data[j % BLOCK_SIZE], exponents, result);
			}
		}
		
	}
	delete [] exponents;*/
}

#endif
