#include "PolyTrie.h"
#include <stdio.h>
#include <sstream>

//Loads a string by parsing it as a sum of monomials
//monomial sum: c_{1}*(x_{1}^e_{1}...x_{varCount}^e_{varCount}) + ...
//nested lists: [[c_{1}, [e_{1}, e_{2}, ..., e_{varCount}]], .. ]
void loadMonomials(monomialSum &monomials, const string &line)
{
	monomials.termCount = 0;
	MonomialLoadConsumer<ZZ>* myLoader = new MonomialLoadConsumer<ZZ>();	
	myLoader->setMonomialSum(monomials);
	parseMonomials(myLoader, line);
}

void parseMonomials(MonomialConsumer<ZZ>* consumer, const string &line)
{
	int varCount = 0;
	for (int i = 0; line[i] != ']'; i++)
	{ varCount += (line[i] == ','); }
	if (varCount < 1)
	{ cout << "There are " << varCount << " variables, bailing." << endl; return; }
	consumer->setDimension(varCount);
	
	int termIndex, lastPos, expIndex, flag;
	termIndex = lastPos = flag = 0; //0 means we expect coefficient, 1 means we expect exponent vector
	
	int *exponents = new int[varCount];
	ZZ coefficient;

	for (int i = 1; i < line.length() - 1; i++) //ignore outermost square brackets
	{
		if (line [i] == '[')
		{
		switch (flag)
		{
			case 0: //coefficient
				lastPos = i + 1; 
				for (; line[i] != ','; i++);
				coefficient = to_ZZ(line.substr(lastPos, i - lastPos).c_str());
				flag = 1;
				break;
			case 1: //exponent vector
				expIndex = 0;
				for (i++; line[i] != ']'; i++)
				{
					if (line[i] != ' ')
					{
						lastPos = i;
						for (; line[i] != ',' && line[i] != ']'; i++);
						exponents[expIndex++] = atoi(line.substr(lastPos, i - lastPos).c_str());
					}
				}
				consumer->ConsumeMonomial(coefficient, exponents);
				flag = 0;
				break;
			default: //error
				cout << "Flag is " << flag << ", bailing." << endl;
				return;
		}
		}
	}

	delete [] exponents;
}

//watch the magic happen
template <class T>
void insertMonomial(const T& coefficient, int* exponents, monomialSum& monomials)
{
	cout << "Inserting " << coefficient << ", [";
	for (int i = 0; i < monomials.varCount; i++)
	{
		cout << exponents[i] << ", ";
	}
	cout << "]" << endl;
	
	burstTrie<T, int>* curTrie;
	
	if (monomials.termCount == 0) //need to construct the first burst trie (sorted on the first variable) and first container 
	{
		cout << "Creating trie" << endl;
		monomials.myMonomials = (burstTrie<T, int>*) malloc (sizeof(burstTrie<T, int>));
		curTrie = monomials.myMonomials;
		curTrie->range = new int[2];
		cout << "Setting upper and lower bound to " << exponents[0] << endl;
		curTrie->range[1] = curTrie->range[0] = exponents[0];
		cout << "Creating first container" << endl;
		curTrie->firstCont = createContainer<T, int>();
	}
	else
	{
		curTrie = monomials.myMonomials;
	}
	
	int termLength = monomials.varCount - 1;
	while (termLength >= 1 && exponents[termLength] == 0) { termLength--; }
	cout << "Exponent vector length is " << termLength << endl;
	cout << "Sorting on: " << exponents[0] << endl;
	cout << "Bounds are [" << curTrie->range[0] << ", " << curTrie->range[1] << "]" << endl;
	
	term<T, int>* curTerm = (term<T, int>*) malloc (sizeof(term<T, int>));
	curTerm->coef = new T(coefficient);
	curTerm->expCount = termLength + 1;
	curTerm->exponents = new int[curTerm->expCount];
	for (int i = 0; i < curTerm->expCount; i++)
	{ curTerm->exponents[i] = exponents[i]; } //copy exponent values
	curTerm->degree = -1;
		
	trieInsert(curTrie, curTerm);
	destroyTerm(curTerm);
		
	monomials.termCount++;
}

//Prints a nested list representation of our sum of monomials
//monomial sum: c_{1}*(x_{1}^e_{1}...x_{varCount}^e_{varCount}) + ...
//nested lists: [[c_{1}, [e_{1}, e_{2}, ..., e_{varCount}]], .. ]
string printMonomials(const monomialSum &myPoly)
{
	TriePrinter<ZZ, int>* myPrinter = new TriePrinter<ZZ, int>();
	myPrinter->setDimension(myPoly.varCount);
	string output = myPrinter->printTrie(myPoly.myMonomials);
	delete myPrinter;
	return "[" + output + "]";
}

//Deallocates space and nullifies internal pointers and counters
void destroyMonomials(monomialSum &myPoly)
{
	destroyTrie(myPoly.myMonomials);
	myPoly.myMonomials = NULL;
	myPoly.termCount = myPoly.varCount = 0;
}

void loadLinForms(linFormSum &forms, const string line)
{
	forms.termCount = 0;
	FormLoadConsumer<ZZ>* myLoader = new FormLoadConsumer<ZZ>();	
	myLoader->setFormSum(forms);
	parseLinForms(myLoader, line);
}
//Loads a string by parsing it as a sum of linear forms
//linear form: (c_{1} / d_{1}!)[(p_{1}*x_{1} + ... p_{varCount}*x_{varCount})^d_{1}] + ...
//nested list: [[c_{1}, [d_{1}, [p_{1}, p_{2}, ..., p_{varCount}]], .. ]
void parseLinForms(FormSumConsumer<ZZ>* consumer, const string& line)
{
	int termIndex = 0;
	int lastPos = 0;
	int varCount = 0;
	int k;
	int flag = 0; //0 means we expect coefficient, 1 means we expect degree, 2 means we expect coefficient vector
	
	for (int i = 0; line[i] != ']'; i++)
	{ varCount += (line[i] == ','); } 
	//varCount is now the number of commas in a linear form - there is 1 less variable;
	varCount--;
	if (varCount < 1)
	{ cout << "There are " << varCount << " variables, bailing." << endl; return; }
	consumer->setDimension(varCount);
	
	vec_ZZ coefs;
	coefs.SetLength(varCount);
	int degree;
	ZZ coefficient;

	for (int i = 1; i < line.length() - 1; i++) //ignore outermost square brackets
	{
		if (line [i] == '[')
		{
		switch (flag)
		{
			case 0: //coefficient
				lastPos = i + 1; 
				for (; line[i] != ','; i++);
				coefficient = to_ZZ(line.substr(lastPos, i - lastPos).c_str());
				flag = 1;
				break;
			case 1: //degree
				lastPos = i + 1;
				for (; line[i] != ','; i++);
				degree = atoi(line.substr(lastPos, i - lastPos).c_str());
				flag = 2;
				break;
			case 2: //coefficient vector
				k = 0;
				for (i++; line[i] != ']'; i++)
				{
					if (line[i] != ' ')
					{
						lastPos = i;
						for (; line[i] != ',' && line[i] != ']'; i++);
						coefs[k++] = to_ZZ(line.substr(lastPos, i - lastPos).c_str());
					}
				}
				for (int j = 1; j <= degree; j++) { coefficient *= j; } //in linFormSum, coefficient is assumed to be divided by the factorial of the form degree
				consumer->ConsumeLinForm(coefficient, degree, coefs);
				flag = 0;
				break;
			default: //error
				cout << "Flag is " << flag << ", bailing." << endl;
				return;
		}
		}
	}
}

// Attempts to find a linear form in formSum with same degree and coefficients as those passed in
// if found, the linear form coefficient in formSum is incremented by coef
// if not found, a new linear form term is added and formSum's termCount is incremented
// if formSum is empty, this sets formSum's lHead and cHead variables
template <class T>
void insertLinForm(const T& coef, int degree, const vec_ZZ& coeffs, linFormSum& formSum) //sort on degree first or last?
{
	cout << "Inserting " << coef << ", [";
	for (int i = 0; i < formSum.varCount; i++)
	{
		cout << coeffs[i] << ", ";
	}
	cout << "]" << endl;
	
	burstTrie<T, ZZ> *curTrie;
	
	if (formSum.termCount == 0) //need to construct the first burst trie (sorted on the first variable) and first container 
	{
		cout << "Creating trie" << endl;
		formSum.myForms = (burstTrie<T, ZZ>*) malloc (sizeof(burstTrie<T, ZZ>));
		curTrie = formSum.myForms;
		cout << "Creating range" << endl;
		curTrie->range = new ZZ[2];
		cout << "Setting upper and lower bound to " << coeffs[0] << endl;
		curTrie->range[0] = coeffs[0];
		curTrie->range[1] = coeffs[0];
		cout << "Creating first container" << endl;
		curTrie->firstCont = createContainer<T, ZZ>();
	}
	else
	{
		curTrie = formSum.myForms;
	}
	
	int termLength = formSum.varCount - 1;
	while (termLength >= 1 && IsZero(coeffs[termLength])) { termLength--; }
	cout << "Term of length " << termLength << endl;
	cout << "Sorting on: " << coeffs[0] << endl;
	cout << "Bounds are [" << curTrie->range[0] << ", " << curTrie->range[1] << "]" << endl;
	
	term<T, ZZ>* curTerm = (term<T, ZZ>*) malloc (sizeof(term<T, ZZ>));
	curTerm->coef = new T(coef);
	curTerm->expCount = termLength + 1;
	curTerm->exponents = new ZZ[curTerm->expCount];
	for (int i = 0; i < curTerm->expCount; i++)
	{ curTerm->exponents[i] = coeffs[i]; } //copy exponent values
	curTerm->degree = degree;
		
	trieInsert(curTrie, curTerm);
	destroyTerm(curTerm);
	formSum.termCount++;
}

//Prints a nested list representation of our sum of linear forms
//linear form: (c_{1} / d_{1}!)[(p_{1}*x_{1} + ... p_{varCount}*x_{varCount})^d_{1}] + ...
//nested list: [[c_{1}, [d_{1}, [p_{1}, p_{2}, ..., p_{varCount}]], .. ]
string printLinForms(const linFormSum &myForm)
{
	TriePrinter<ZZ, ZZ>* myPrinter = new TriePrinter<ZZ, ZZ>();
	myPrinter->setDimension(myForm.varCount);
	string output = myPrinter->printTrie(myForm.myForms);
	delete myPrinter;
	return "[" + output + "]";
}

//Deallocates space and nullifies internal pointers and counters
void destroyLinForms(linFormSum &myPoly)
{
	destroyTrie(myPoly.myForms);
	myPoly.myForms = NULL;
	myPoly.termCount = myPoly.varCount = 0;
}

//INPUT: monomial specified by myPoly.coefficientBlocks[mIndex / BLOCK_SIZE].data[mIndex % BLOCK_SIZE]
//	and myPoly.exponentBlocks[mIndex / BLOCK_SIZE].data[mIndex % BLOCK_SIZE]
//OUTPUT: lForm now also contains the linear decomposition of this monomial 
//	note: all linear form coefficients assumed to be divided by their respective |M|!, and the form is assumed to be of power M
void decompose(monomialSum &myPoly, linFormSum &lForm, int mIndex)
{
	Decomposer* myDecomposer = new Decomposer();
	lForm.myForms = (burstTrie<ZZ, ZZ>*) malloc(sizeof(burstTrie<ZZ, ZZ>));
	myDecomposer->init(&myPoly, &lForm);
	myDecomposer->setDimension(myPoly.varCount);
	myDecomposer->enumTrie(myPoly.myMonomials);
	delete myDecomposer;
}
