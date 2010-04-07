#include "PolyTrie.h"
#include <stdio.h>
#include <sstream>

#define PT_DEBUG 0

//Loads a string by parsing it as a sum of monomials
//monomial sum: c_{1}*(x_{1}^e_{1}...x_{varCount}^e_{varCount}) + ...
//nested lists: [[c_{1}, [e_{1}, e_{2}, ..., e_{varCount}]], .. ]
void loadMonomials(monomialSum &monomials, const string &line)
{
	monomials.termCount = 0;
	MonomialLoadConsumer<ZZ>* myLoader = new MonomialLoadConsumer<ZZ>();	
	myLoader->setMonomialSum(monomials);
	parseMonomials(myLoader, line);
	delete myLoader;
}

void parseMonomials(MonomialConsumer<ZZ>* consumer, const string &line)
{
	int varCount = 0;
	for (int i = 0; line[i] != ']'; i++)
	{ varCount += (line[i] == ','); }
	if (varCount < 1)
	{ cout << "line: `" << line << "'" << endl; cout << "There are " << varCount << " variables, bailing." << endl; return; }
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
void insertMonomial(const ZZ& coefficient, int* exponents, monomialSum& monomials)
{	
	BurstTrie<ZZ, int>* curTrie;
	
	if (monomials.termCount == 0) //need to construct the first burst trie (sorted on the first variable) and first container 
	{
		if (PT_DEBUG) { cout << "Creating trie" << endl; }
		monomials.myMonomials = new BurstTrie<ZZ, int>();
		curTrie = monomials.myMonomials;
	}
	else
	{
		curTrie = monomials.myMonomials;
	}
	
	//int termLength = monomials.varCount - 1;
	//while (termLength >= 1 && exponents[termLength] == 0) { termLength--; }
	//if (PT_DEBUG) { cout << "Exponent vector length is " << termLength << endl;
	//cout << "Sorting on: " << exponents[0] << endl;
	//cout << "Bounds are [" << curTrie->range[0] << ", " << curTrie->range[1] << "]" << endl; }
	
	/*term<T, int>* curTerm = (term<T, int>*) malloc (sizeof(term<T, int>));
	curTerm->coef = new T(coefficient);
	curTerm->expCount = termLength + 1;
	curTerm->exponents = new int[curTerm->expCount];
	for (int i = 0; i < curTerm->expCount; i++)
	{ curTerm->exponents[i] = exponents[i]; } //copy exponent values
	curTerm->degree = -1;*/
		
	//trieInsert(curTrie, curTerm);
	//destroyTerm(curTerm);
	curTrie->insertTerm(coefficient, exponents, 0, monomials.varCount, -1);
		
	monomials.termCount++;
}

//Prints a nested list representation of our sum of monomials
//monomial sum: c_{1}*(x_{1}^e_{1}...x_{varCount}^e_{varCount}) + ...
//nested lists: [[c_{1}, [e_{1}, e_{2}, ..., e_{varCount}]], .. ]
string printMonomials(const monomialSum &myPoly)
{
	BurstTerm<ZZ, int>* temp = new BurstTerm<ZZ, int>(myPoly.varCount);
	BurstTrie<ZZ, int>* myTrie = myPoly.myMonomials;
	myTrie->begin();
	int i = 0;
	stringstream output (stringstream::in | stringstream::out);
	while (myTrie->nextTerm(temp))
	{
		if (output.str() != "")
		{ output << ", "; }
		output << "[" << temp->coef << ", [";
		for (int j = 0; j < temp->length; j++)
		{
			output << temp->exps[j];
			if (j + 1 < temp->length)
			{ output << ", "; }
		}
		output << "]]";
	}
	delete temp;
	return "[" + output.str() + "]";
}

//Deallocates space and nullifies internal pointers and counters
void destroyMonomials(monomialSum &myPoly)
{
	delete myPoly.myMonomials;
	myPoly.myMonomials = NULL;
	myPoly.termCount = myPoly.varCount = 0;
}

void loadLinForms(linFormSum &forms, const string line)
{
	forms.termCount = 0;
	FormLoadConsumer<ZZ>* myLoader = new FormLoadConsumer<ZZ>();	
	myLoader->setFormSum(forms);
	parseLinForms(myLoader, line);
	delete myLoader;
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
	{ cout << "line: `" << line << "'" << endl; cout << "There are " << varCount << " variables, bailing." << endl; return; }
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
void insertLinForm(const ZZ& coef, int degree, const vec_ZZ& coeffs, linFormSum& formSum) //sort on degree first or last?
{
	BurstTrie<ZZ, ZZ> *curTrie;
	
	if (formSum.termCount == 0) //need to construct the first burst trie (sorted on the first variable) and first container 
	{
		formSum.myForms = new BurstTrie<ZZ, ZZ>();
		curTrie = formSum.myForms;
	}
	else
	{
		curTrie = formSum.myForms;
	}
	
	/*int termLength = formSum.varCount - 1;
	while (termLength >= 1 && IsZero(coeffs[termLength])) { termLength--; }
	if (PT_DEBUG) { cout << "Term of length " << termLength << endl;
	cout << "Sorting on: " << coeffs[0] << endl;
	cout << "Bounds are [" << curTrie->range[0] << ", " << curTrie->range[1] << "]" << endl; }
	
	term<T, ZZ>* curTerm = (term<T, ZZ>*) malloc (sizeof(term<T, ZZ>));
	curTerm->coef = new T(coef);
	curTerm->expCount = termLength + 1;
	curTerm->exponents = new ZZ[curTerm->expCount];
	for (int i = 0; i < curTerm->expCount; i++)
	{ curTerm->exponents[i] = coeffs[i]; } //copy exponent values
	curTerm->degree = degree;
		
	trieInsert(curTrie, curTerm);
	destroyTerm(curTerm);*/
	ZZ* exps = new ZZ[formSum.varCount];
	for (int i = 0; i < formSum.varCount; i++) { exps[i] = coeffs[i]; }
	curTrie->insertTerm(coef, exps, 0, formSum.varCount, degree);
	delete [] exps;
	formSum.termCount++;
}

//Prints a nested list representation of our sum of linear forms
//linear form: (c_{1} / d_{1}!)[(p_{1}*x_{1} + ... p_{varCount}*x_{varCount})^d_{1}] + ...
//nested list: [[c_{1}, [d_{1}, [p_{1}, p_{2}, ..., p_{varCount}]], .. ]
string printLinForms(const linFormSum &myForm)
{
	BurstTerm<ZZ, ZZ>* temp = new BurstTerm<ZZ, ZZ>(myForm.varCount);
	BurstTrie<ZZ, ZZ>* myTrie = myForm.myForms;
	myTrie->begin();
	int i = 0;
	stringstream output (stringstream::in | stringstream::out);
	while (myTrie->nextTerm(temp))
	{
		if (output.str() != "")
		{ output << ", "; }
		output << "[" << temp->coef << ", [" << temp->degree << ", [";
		for (int j = 0; j < temp->length; j++)
		{
			output << temp->exps[j];
			if (j + 1 < temp->length)
			{ output << ", "; }
		}
		output << "]]]";
	}
	delete temp;
	return "[" + output.str() + "]";
}

//Deallocates space and nullifies internal pointers and counters
void destroyLinForms(linFormSum &myPoly)
{
	delete myPoly.myForms;
	myPoly.myForms = NULL;
	myPoly.termCount = myPoly.varCount = 0;
}

//INPUT: monomial specified by myPoly.coefficientBlocks[mIndex / BLOCK_SIZE].data[mIndex % BLOCK_SIZE]
//	and myPoly.exponentBlocks[mIndex / BLOCK_SIZE].data[mIndex % BLOCK_SIZE]
//OUTPUT: lForm now also contains the linear decomposition of this monomial 
//	note: all linear form coefficients assumed to be divided by their respective |M|!, and the form is assumed to be of power M
void decompose(monomialSum &myPoly, linFormSum &lForm)
{
	BurstTerm<ZZ, int>* temp = new BurstTerm<ZZ, int>(myPoly.varCount);
	lForm.myForms = new BurstTrie<ZZ, ZZ>();
	BurstTrie<ZZ, int>* myTrie = myPoly.myMonomials;
	myTrie->begin();
	
	while (myTrie->nextTerm(temp))
	{
		decompose(myPoly, lForm, temp);
	}
	
	delete temp;
	//delete myPoly.myMonomials;
}

void decompose(monomialSum &myPoly, linFormSum &lForm, BurstTerm<ZZ, int>* myTerm)
{
	vec_ZZ myExps; myExps.SetLength(lForm.varCount);
	
	if (myTerm->length == 0) //constant
	{
		for (int j = 0; j < lForm.varCount; j++) { myExps[j] = 0; }
		insertLinForm(myTerm->coef, 0, myExps, lForm);
		return;
	}
	
	ZZ formsCount = to_ZZ(myTerm->exps[0] + 1); //first exponent
	int totalDegree = myTerm->exps[0];
	for (int i = 1; i < myTerm->length; i++)
	{
		formsCount *=  myTerm->exps[i] + 1;
		totalDegree +=  myTerm->exps[i];
	}
	formsCount--;
	cout << "At most " << formsCount << " linear forms will be required for this decomposition." << endl;
	cout << "Total degree is " << totalDegree << endl;
	
	int* p = new int[myPoly.varCount];
	int* counter = new int[myPoly.varCount];
	ZZ* binomCoeffs = new ZZ[myPoly.varCount]; //for calculating the product of binomial coefficients for each linear form

	ZZ temp;
	int g;
	int myIndex;
	for (int i = 0; i < myPoly.varCount; i++) { counter[i] = 0; binomCoeffs[i] = to_ZZ(1); }
	for (ZZ i = to_ZZ(1); i <= formsCount; i++)
	{
		//cout << "i is " << i << endl;
		counter[0] += 1;
		for (myIndex = 0; counter[myIndex] > myTerm->exps[myIndex]; myIndex++)
		{
			counter[myIndex] = 0;
			binomCoeffs[myIndex] = to_ZZ(1);
			counter[myIndex+1] += 1;
		}
		binomCoeffs[myIndex] *= myTerm->exps[myIndex] - counter[myIndex] + 1; // [n choose k] = [n choose (k - 1) ] * (n - k + 1)/k
		binomCoeffs[myIndex] /= counter[myIndex];
		//cout << "counter is: " << counter << endl;
		
		//find gcd of all elements and calculate the product of binomial coefficients for the linear form
		g = counter[0];
		int parity = totalDegree - counter[0];
		p[0] = counter[0];
		temp = binomCoeffs[0];
		if (myPoly.varCount > 1)
		{
			for (int k = 1; k < myPoly.varCount; k++)
			{
				p[k] = counter[k];
				g = GCD (g, p[k]);
				parity -= p[k];
				temp *= binomCoeffs[k];
			}
		}
		
		//calculate coefficient
		temp *= myTerm->coef;
		if ((parity % 2) == 1) { temp *= -1; } // -1 ^ [|M| - (p[0] + p[1] + ... p[n])], checks for odd parity using modulo 2
		
		if (g != 1)
		{
			for (int k = 0; k < myPoly.varCount; k++)
			{
				p[k] /= g;
			}
			temp *= power_ZZ(g, totalDegree);
		}
		//cout << "coefficient is " << temp << endl;
		for (int i = 0; i < myPoly.varCount; i++) { myExps[i] = p[i]; }
		insertLinForm(temp, totalDegree, myExps, lForm);
	}
	delete [] p;
	delete [] counter;
	delete [] binomCoeffs;
}
