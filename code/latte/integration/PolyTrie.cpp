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
	MonomialLoadConsumer<RationalNTL>* myLoader = new MonomialLoadConsumer<
			RationalNTL> ();
	myLoader->setMonomialSum(monomials);
	parseMonomials(myLoader, line);
	delete myLoader;
}

void parseMonomials(MonomialConsumer<RationalNTL>* consumer, const string &line)
{
	int varCount = 0;
	for (int i = 0; line[i] != ']'; i++)
	{
		varCount += (line[i] == ',');
	}
	if (varCount < 1)
	{
		cout << "line: `" << line << "'" << endl;
		cout << "There are " << varCount << " variables, bailing." << endl;
		return;
	}
	consumer->setDimension(varCount);

	int termIndex, lastPos, expIndex, flag;
	termIndex = lastPos = flag = 0; //0 means we expect coefficient, 1 means we expect exponent vector

	int *exponents = new int[varCount];
	RationalNTL coefficient;

	for (int i = 1; i < line.length() - 1; i++) //ignore outermost square brackets
	{
		if (line[i] == '[')
		{
			switch (flag)
			{
				case 0: //coefficient
					lastPos = i + 1;
					for (; line[i] != ','; i++)
						;
					coefficient = RationalNTL(
							line.substr(lastPos, i - lastPos).c_str());
					flag = 1;
					break;
				case 1: //exponent vector
					expIndex = 0;
					for (i++; line[i] != ']'; i++)
					{
						if (line[i] != ' ')
						{
							lastPos = i;
							for (; line[i] != ',' && line[i] != ']'; i++)
								;
							exponents[expIndex++] = atoi(line.substr(lastPos, i
									- lastPos).c_str());
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

	delete[] exponents;
}

//watch the magic happen
void insertMonomial(const RationalNTL& coefficient, int* exponents,
		monomialSum& monomials)
{
	BurstTrie<RationalNTL, int>* curTrie;

	if (monomials.termCount == 0) //need to construct the first burst trie (sorted on the first variable) and first container 
	{
		if (PT_DEBUG)
		{
			cout << "Creating trie" << endl;
		}
		monomials.myMonomials = new BurstTrie<RationalNTL, int> ();
		curTrie = monomials.myMonomials;
	} else
	{
		curTrie = monomials.myMonomials;
	}

	curTrie->insertTerm(coefficient, exponents, 0, monomials.varCount, -1);

	monomials.termCount++;
}

//Prints a nested list representation of our sum of monomials
//monomial sum: c_{1}*(x_{1}^e_{1}...x_{varCount}^e_{varCount}) + ...
//nested lists: [[c_{1}, [e_{1}, e_{2}, ..., e_{varCount}]], .. ]
string printMonomials(const monomialSum &myPoly)
{
	BTrieIterator<RationalNTL, int>* it =
			new BTrieIterator<RationalNTL, int> ();
	term<RationalNTL, int>* temp;
	it->setTrie(myPoly.myMonomials, myPoly.varCount);
	it->begin();
	int i = 0;
	stringstream output(stringstream::in | stringstream::out);
	temp = it->nextTerm();
	do
	{
		if (output.str() != "")
		{
			output << ", ";
		}
		output << "[" << temp->coef << ", [";
		for (int j = 0; j < temp->length; j++)
		{
			output << temp->exps[j];
			if (j + 1 < temp->length)
			{
				output << ", ";
			}
		}
		output << "]]";
		temp = it->nextTerm();
	} while (temp);
	delete it;
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
	FormLoadConsumer<RationalNTL>* myLoader =
			new FormLoadConsumer<RationalNTL> ();
	myLoader->setFormSum(forms);
	parseLinForms(myLoader, line);
	delete myLoader;
}
//Loads a string by parsing it as a sum of linear forms
//linear form: (c_{1} / d_{1}!)[(p_{1}*x_{1} + ... p_{varCount}*x_{varCount})^d_{1}] + ...
//nested list: [[c_{1}, [d_{1}, [p_{1}, p_{2}, ..., p_{varCount}]], .. ]
void parseLinForms(FormSumConsumer<RationalNTL>* consumer, const string& line)
{
	int termIndex = 0;
	int lastPos = 0;
	int varCount = 0;
	int k;
	int flag = 0; //0 means we expect coefficient, 1 means we expect degree, 2 means we expect coefficient vector

	//cout << "parseLinForms: line = " << line.c_str() << endl;
	for (int i = 0; line[i] != ']'; i++)
	{
		varCount += (line[i] == ',');
	}
	//varCount is now the number of commas in a linear form - there is 1 less variable;
	varCount--;
	if (varCount < 1)
	{
		cout << "line: `" << line << "'" << endl;
		cout << "There are " << varCount << " variables, bailing." << endl;
		return;
	}
	consumer->setDimension(varCount);

	vec_ZZ coefs;
	coefs.SetLength(varCount);
	int degree;
	RationalNTL coefficient;

	for (int i = 1; i < line.length() - 1; i++) //ignore outermost square brackets
	{
		if (line[i] == '[')
		{
			ZZ degreeFactorial;
			switch (flag)
			{
				case 0: //coefficient
					lastPos = i + 1;
					for (; line[i] != ','; i++)
						;
					coefficient = RationalNTL(
							(line.substr(lastPos, i - lastPos).c_str()));
					flag = 1;
					break;
				case 1: //degree
					lastPos = i + 1;
					for (; line[i] != ','; i++)
						;
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
							for (; line[i] != ',' && line[i] != ']'; i++)
								;

							coefs[k++] = to_ZZ(
									line.substr(lastPos, i - lastPos).c_str());

						}
					}
					degreeFactorial = 1;
					for (int j = 1; j <= degree; j++)
					{
						degreeFactorial *= j;
						//coefficient *= j;
					} //in linFormSum, coefficient is assumed to be divided by the factorial of the form degree
					coefficient *= degreeFactorial;
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
void insertLinForm(const RationalNTL& coef, int degree, const vec_ZZ& coeffs,
		linFormSum& formSum) //sort on degree first or last?
{
	BurstTrie<RationalNTL, ZZ> *curTrie;
	//cout << "inserting into linear form with " << formSum.varCount << " variables" << endl;
	if (formSum.termCount == 0) //need to construct the first burst trie (sorted on the first variable) and first container 
	{
		formSum.myForms = new BurstTrie<RationalNTL, ZZ> ();
		curTrie = formSum.myForms;
	} else
	{
		curTrie = formSum.myForms;
	}

	ZZ* exps = new ZZ[formSum.varCount];
	for (int i = 0; i < formSum.varCount; i++)
	{
		exps[i] = coeffs[i];
	}
	curTrie->insertTerm(coef, exps, 0, formSum.varCount, degree);

	delete[] exps;
	formSum.termCount++;
}

//Prints a nested list representation of our sum of linear forms
//linear form: (c_{1} / d_{1}!)[(p_{1}*x_{1} + ... p_{varCount}*x_{varCount})^d_{1}] + ...
//nested list: [[c_{1}, [d_{1}, [p_{1}, p_{2}, ..., p_{varCount}]], .. ]
string printLinForms(const linFormSum &myForm)
{
	BTrieIterator<RationalNTL, ZZ>* it = new BTrieIterator<RationalNTL, ZZ> ();
	term<RationalNTL, ZZ>* temp;
	it->setTrie(myForm.myForms, myForm.varCount);
	it->begin();

	stringstream output(stringstream::in | stringstream::out);
	temp = it->nextTerm();
	do
	{
		if (output.str() != "")
		{
			output << ", ";
		}
		output << "[" << temp->coef << ", [" << temp->degree << ", [";
		for (int j = 0; j < temp->length; j++)
		{
			output << temp->exps[j];
			if (j + 1 < temp->length)
			{
				output << ", ";
			}
		}
		output << "]]]";
		temp = it->nextTerm();
	} while (temp);
	delete it;
	return "[" + output.str() + "]";
}

//Deallocates space and nullifies internal pointers and counters
void destroyLinForms(linFormSum &myPoly)
{
	if (myPoly.myForms)
		delete myPoly.myForms;
	myPoly.myForms = NULL;
	myPoly.termCount = myPoly.varCount = 0;
}

//INPUT: monomial specified by myPoly.coefficientBlocks[mIndex / BLOCK_SIZE].data[mIndex % BLOCK_SIZE]
//	and myPoly.exponentBlocks[mIndex / BLOCK_SIZE].data[mIndex % BLOCK_SIZE]
//OUTPUT: lForm now also contains the linear decomposition of this monomial 
//	note: all linear form coefficients assumed to be divided by their respective |M|!, and the form is assumed to be of power M
void decompose(BTrieIterator<RationalNTL, int>* it, linFormSum &lForm)
{
	//cout << "decomposing " << lForm.varCount << " variables" << endl;
	term<RationalNTL, int>* temp;
	//BTrieIterator<ZZ, int>* it = new BTrieIterator<ZZ, int>();


	//it->setTrie(myPoly.myMonomials, myPoly.varCount);
	it->begin();

	//cout << "in decompost()::\n";
	temp = it->nextTerm();

	do
	{
		//cout << "monomial " << temp->coef;
		//for(int i = 0; i < temp->length; ++i)
		//	cout << temp->exps[i];
		//cout << "\n";
		decompose(temp, lForm);
		temp = it->nextTerm();
	} while (temp);

	//delete myPoly.myMonomials;
}

void decompose(term<RationalNTL, int>* myTerm, linFormSum &lForm)
{
	vec_ZZ myExps;
	myExps.SetLength(lForm.varCount);

	//April 27. 2011 Brandon: I don't think this if is ever true. If I insert [[[1,[0,0]]], I get a monomial of length 2 still.
	//This or the "string to polynomial" function is not  correct. I think we can delete this if statement.
	if (myTerm->length == 0) //constant
	{
		for (int j = 0; j < lForm.varCount; j++)
		{
			myExps[j] = 0;
		}
		insertLinForm(myTerm->coef, 0, myExps, lForm);
		return;
	}


	ZZ formsCount = to_ZZ(myTerm->exps[0] + 1); //first exponent
	int totalDegree = myTerm->exps[0];
	for (int i = 1; i < myTerm->length; i++)
	{
		formsCount *= myTerm->exps[i] + 1;
		totalDegree += myTerm->exps[i];
	}
	formsCount--;

	//If this is zero, then the term is a constant [c,[0,0,0,0...0]]
	if ( formsCount == 0)
	{
		for (int j = 0; j < lForm.varCount; j++)
		{
			myExps[j] = 0;
		}
		insertLinForm(myTerm->coef, 0, myExps, lForm);
		return;
	}//if formsCount

	//cout << "At most " << formsCount << " linear forms will be required for this decomposition." << endl;
	//cout << "Total degree is " << totalDegree << endl;

	int* p = new int[lForm.varCount];
	int* counter = new int[lForm.varCount];
	ZZ* binomCoeffs = new ZZ[lForm.varCount]; //for calculating the product of binomial coefficients for each linear form

	RationalNTL temp;
	int g;
	int myIndex;
	for (int i = 0; i < lForm.varCount; i++)
	{
		counter[i] = 0;
		binomCoeffs[i] = to_ZZ(1);
	}
	for (ZZ i = to_ZZ(1); i <= formsCount; i++)
	{
		//cout << "i is " << i << endl;
		counter[0] += 1;
		for (myIndex = 0; counter[myIndex] > myTerm->exps[myIndex]; myIndex++)
		{
			counter[myIndex] = 0;
			binomCoeffs[myIndex] = to_ZZ(1);
			counter[myIndex + 1] += 1;
		}
		binomCoeffs[myIndex] *= myTerm->exps[myIndex] - counter[myIndex] + 1; // [n choose k] = [n choose (k - 1) ] * (n - k + 1)/k
		binomCoeffs[myIndex] /= counter[myIndex];
		//cout << "counter is: " << counter << endl;

		//find gcd of all elements and calculate the product of binomial coefficients for the linear form
		g = counter[0];
		int parity = totalDegree - counter[0];
		p[0] = counter[0];
		temp = binomCoeffs[0];
		if (lForm.varCount > 1)
		{
			for (int k = 1; k < lForm.varCount; k++)
			{
				p[k] = counter[k];
				g = GCD(g, p[k]);
				parity -= p[k];
				temp *= binomCoeffs[k];
			}
		}

		//calculate coefficient
		temp *= myTerm->coef;
		if ((parity % 2) == 1)
		{
			temp *= to_ZZ(-1);
		} // -1 ^ [|M| - (p[0] + p[1] + ... p[n])], checks for odd parity using modulo 2

		if (g != 1)
		{
			for (int k = 0; k < lForm.varCount; k++)
			{
				p[k] /= g;
			}
			temp *= power_ZZ(g, totalDegree);
		}
		//cout << "coefficient is " << temp << endl;
		for (int i = 0; i < lForm.varCount; i++)
		{
			myExps[i] = p[i];
		}
		insertLinForm(temp, totalDegree, myExps, lForm);
	}
	delete[] p;
	delete[] counter;
	delete[] binomCoeffs;
}
