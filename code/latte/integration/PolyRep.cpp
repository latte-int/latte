#include "PolyRep.h"
#include <stdio.h>
#include <sstream>

//Loads a string by parsing it as a sum of monomials
//monomial sum: c_{1}*(x_{1}^e_{1}...x_{varCount}^e_{varCount}) + ...
//nested lists: [[c_{1}, [e_{1}, e_{2}, ..., e_{varCount}]], .. ]

void _loadMonomials(_monomialSum &monomials, const string &line)
{
	monomials.termCount = 0;
	_MonomialLoadConsumer<RationalNTL>* myLoader = new _MonomialLoadConsumer<
			RationalNTL> ();
	myLoader->setMonomialSum(monomials);
	_parseMonomials(myLoader, line);
}

void _parseMonomials(_MonomialConsumer<RationalNTL>* consumer,
		const string &line)
{
	int varCount = 0;
	for (int i = 0; line[i] != ']'; i++)
	{
		varCount += (line[i] == ',');
	}
	if (varCount < 1)
	{
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
					coefficient = RationalNTL(line.substr(lastPos, i - lastPos).c_str());
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

//Prints a nested list representation of our sum of monomials
//monomial sum: c_{1}*(x_{1}^e_{1}...x_{varCount}^e_{varCount}) + ...
//nested lists: [[c_{1}, [e_{1}, e_{2}, ..., e_{varCount}]], .. ]
string _printMonomials(const _monomialSum &myPoly)
{
	stringstream output(stringstream::in | stringstream::out);
	output << "[";
	eBlock* expTmp = myPoly.eHead;
	cBlock<RationalNTL>* coeffTmp = myPoly.cHead;
	int termCount = 0;
	do
	{
		for (int i = 0; i < BLOCK_SIZE && termCount < myPoly.termCount; i++)
		{
			output << "[" << coeffTmp->data[i] << ",[";
			for (int j = (i * myPoly.varCount); j < ((i + 1) * myPoly.varCount); j++)
			{
				output << expTmp->data[j];
				if (j + 1 < ((i + 1) * myPoly.varCount))
				{
					output << ",";
				}
			}
			output << "]]";
			if (termCount + 1 < myPoly.termCount)
			{
				output << ",";
			}
			termCount++;
		}
		coeffTmp = coeffTmp->next;
		expTmp = expTmp->next;
	} while (coeffTmp != NULL);
	output << "]";
	return output.str();
}

//Deallocates space and nullifies internal pointers and counters
void _destroyMonomials(_monomialSum &myPoly)
{
	eBlock* expTmp = myPoly.eHead;
	cBlock<RationalNTL>* coeffTmp = myPoly.cHead;
	eBlock* oldExp = NULL;
	cBlock<RationalNTL>* oldCoeff = NULL;
	do
	{
		oldExp = expTmp;
		oldCoeff = coeffTmp;
		expTmp = expTmp->next;
		coeffTmp = coeffTmp->next;
		free(oldExp);
		free(oldCoeff);
	} while (coeffTmp != NULL);
	myPoly.eHead = NULL;
	myPoly.cHead = NULL;
	myPoly.termCount = myPoly.varCount = 0;
}

void _loadLinForms(_linFormSum &forms, const string line)
{
	forms.termCount = 0;
	_FormLoadConsumer<RationalNTL>* myLoader = new _FormLoadConsumer<
			RationalNTL> ();
	myLoader->setFormSum(forms);
	_parseLinForms(myLoader, line);
}
//Loads a string by parsing it as a sum of linear forms
//linear form: (c_{1} / d_{1}!)[(p_{1}*x_{1} + ... p_{varCount}*x_{varCount})^d_{1}] + ...
//nested list: [[c_{1}, [d_{1}, [p_{1}, p_{2}, ..., p_{varCount}]], .. ]
void _parseLinForms(_FormSumConsumer<RationalNTL>* consumer, const string& line)
{
	int termIndex = 0;
	int lastPos = 0;
	int varCount = 0;
	int k;
	int flag = 0; //0 means we expect coefficient, 1 means we expect degree, 2 means we expect coefficient vector

	for (int i = 0; line[i] != ']'; i++)
	{
		varCount += (line[i] == ',');
	}
	//varCount is now the number of commas in a linear form - there is 1 less variable;
	varCount--;
	if (varCount < 1)
	{
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
			int degreeFactorial;
			switch (flag)
			{
				case 0: //coefficient
					lastPos = i + 1;
					for (; line[i] != ','; i++)
						;
					coefficient = RationalNTL(line.substr(lastPos, i - lastPos).c_str());
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
					} //in _linFormSum, coefficient is assumed to be divided by the factorial of the form degree
					coefficient *= to_ZZ(degreeFactorial);
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

//Prints a nested list representation of our sum of linear forms
//linear form: (c_{1} / d_{1}!)[(p_{1}*x_{1} + ... p_{varCount}*x_{varCount})^d_{1}] + ...
//nested list: [[c_{1}, [d_{1}, [p_{1}, p_{2}, ..., p_{varCount}]], .. ]
string _printLinForms(const _linFormSum &myForm)
{
	stringstream output(stringstream::in | stringstream::out);
	output << "[";
	lBlock* formTmp = myForm.lHead;
	cBlock<RationalNTL>* coeffTmp = myForm.cHead;
	for (int i = 0; i < myForm.termCount; i++)
	{
		if (i > 0 && i % BLOCK_SIZE == 0)
		{
			formTmp = formTmp->next;
			coeffTmp = coeffTmp->next;
		}

		output << "[" << coeffTmp->data[i % BLOCK_SIZE] << ", ["
				<< formTmp->degree[i % BLOCK_SIZE] << ", [";
		for (int j = 0; j < myForm.varCount; j++)
		{
			output << formTmp->data[i % BLOCK_SIZE][j];
			if (j + 1 < myForm.varCount)
			{
				output << ", ";
			}
		}
		output << "]]]";
		if (i + 1 < myForm.termCount)
		{
			output << ", ";
		}
	}
	output << "]";
	return output.str();
}

//Deallocates space and nullifies internal pointers and counters
void _destroyLinForms(_linFormSum &myPoly)
{
	lBlock* expTmp = myPoly.lHead;
	cBlock<RationalNTL>* coeffTmp = myPoly.cHead;
	lBlock* oldExp = NULL;
	cBlock<RationalNTL>* oldCoeff = NULL;
	int termCount = 0;
	do
	{
		oldExp = expTmp;
		oldCoeff = coeffTmp;
		expTmp = expTmp->next;
		coeffTmp = coeffTmp->next;
		free(oldExp);
		free(oldCoeff);
	} while (coeffTmp != NULL);
	myPoly.lHead = NULL;
	myPoly.cHead = NULL;
	myPoly.termCount = myPoly.varCount = 0;
}

//INPUT: monomial specified by myPoly.coefficientBlocks[mIndex / BLOCK_SIZE].data[mIndex % BLOCK_SIZE]
//	and myPoly.exponentBlocks[mIndex / BLOCK_SIZE].data[mIndex % BLOCK_SIZE]
//OUTPUT: lForm now also contains the linear decomposition of this monomial 
//	note: all linear form coefficients assumed to be divided by their respective |M|!, and the form is assumed to be of power M
void _decompose(_monomialSum &myPoly, _linFormSum &lForm, int mIndex)
{
	eBlock* expTmp = myPoly.eHead;
	cBlock<RationalNTL>* coeffTmp = myPoly.cHead;
	for (int i = 0; i < (mIndex / BLOCK_SIZE); i++)
	{
		expTmp = expTmp->next;
		coeffTmp = coeffTmp->next;
	}

	bool constantTerm = true;
	for (int i = (mIndex * myPoly.varCount); i < ((mIndex + 1)
			* myPoly.varCount); i++)
	{
		if (expTmp->data[i] != 0)
		{
			constantTerm = false;
			break;
		}
	}
	vec_ZZ myExps;
	myExps.SetLength(lForm.varCount);
	if (constantTerm) //exponents are all 0, this is a constant term - linear form is already known
	{
		for (int j = 0; j < lForm.varCount; j++)
		{
			myExps[j] = 0;
		}
		_insertLinForm<RationalNTL> (coeffTmp->data[mIndex % BLOCK_SIZE], 0, myExps,
				lForm);
		return;
	}

	ZZ formsCount = to_ZZ(expTmp->data[(mIndex % BLOCK_SIZE) * myPoly.varCount]
			+ 1);
	int totalDegree = expTmp->data[(mIndex % BLOCK_SIZE) * myPoly.varCount];
	for (int i = 1; i < myPoly.varCount; i++)
	{
		formsCount *= expTmp->data[(mIndex % BLOCK_SIZE) * myPoly.varCount + i]
				+ 1;
		totalDegree
				+= expTmp->data[(mIndex % BLOCK_SIZE) * myPoly.varCount + i];
	}
	formsCount--;
	//cout << "At most " << formsCount << " linear forms will be required for this decomposition." << endl;
	//cout << "Total degree is " << totalDegree << endl;

	int* p = new int[myPoly.varCount];
	int* counter = new int[myPoly.varCount];
	ZZ* binomCoeffs = new ZZ[myPoly.varCount]; //for calculating the product of binomial coefficients for each linear form

	RationalNTL temp;
	int g;
	bool found;
	int myIndex;
	for (int i = 0; i < myPoly.varCount; i++)
	{
		counter[i] = 0;
		binomCoeffs[i] = to_ZZ(1);
	}
	for (ZZ i = to_ZZ(1); i <= formsCount; i++)
	{
		//cout << "i is " << i << endl;
		counter[0] += 1;
		for (myIndex = 0; counter[myIndex] > expTmp->data[(mIndex % BLOCK_SIZE)
				* myPoly.varCount + myIndex]; myIndex++)
		{
			counter[myIndex] = 0;
			binomCoeffs[myIndex] = to_ZZ(1);
			counter[myIndex + 1] += 1;
		}
		binomCoeffs[myIndex] *= expTmp->data[(mIndex % BLOCK_SIZE)
				* myPoly.varCount + myIndex] - counter[myIndex] + 1; // [n choose k] = [n choose (k - 1) ] * (n - k + 1)/k
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
				g = GCD(g, p[k]);
				parity -= p[k];
				temp *= binomCoeffs[k];
			}
		}

		//calculate coefficient
		temp *= coeffTmp->data[mIndex % BLOCK_SIZE];
		if ((parity % 2) == 1)
		{
			temp *= to_ZZ(-1);
		} // -1 ^ [|M| - (p[0] + p[1] + ... p[n])], checks for odd parity using modulo 2

		if (g != 1)
		{
			for (int k = 0; k < myPoly.varCount; k++)
			{
				p[k] /= g;
			}
			temp *= power_ZZ(g, totalDegree);
		}
		//cout << "coefficient is " << temp << endl;
		for (int i = 0; i < myPoly.varCount; i++)
		{
			myExps[i] = p[i];
		}
		_insertLinForm<RationalNTL> (temp, totalDegree, myExps, lForm);
	}
	delete[] p;
	delete[] counter;
	delete[] binomCoeffs;
}
