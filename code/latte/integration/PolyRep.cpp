#include "PolyRep.h"
#include <stdio.h>
#include <sstream>

//Loads a string by parsing it as a sum of monomials
//monomial sum: c_{1}*(x_{1}^e_{1}...x_{varCount}^e_{varCount}) + ...
//nested lists: [[c_{1}, [e_{1}, e_{2}, ..., e_{varCount}]], .. ]
void loadMonomials(monomialSum &monomials, const string line)
{
	int termIndex = 0;
	int lastPos = 0;
	int varCount = 0;
	int k;
	int flag = 0; //0 means we expect coefficient, 1 means we expect exponent vector
	
	for (int i = 0; line[i] != ']'; i++)
	{ varCount += (line[i] == ','); } 
	//varCount is now the number of commas in the monomial, same as number of variables
	if (varCount < 1)
	{ cout << "There are " << varCount << " variables, bailing." << endl; return; }

	//at least one record requires at least one block
	monomials.eHead = (eBlock*) malloc (sizeof(eBlock));
	monomials.eHead->next = NULL;
	monomials.cHead = (cBlock*) malloc (sizeof(cBlock));
	monomials.cHead->next = NULL;

	eBlock* expBlock = monomials.eHead;
	cBlock* coefBlock = monomials.cHead;
	expBlock->data = new vec_ZZ[BLOCK_SIZE];
	coefBlock->data = new ZZ[BLOCK_SIZE];

	for (int i = 1; i < line.length() - 1; i++) //ignore outermost square brackets
	{
		if (line [i] == '[')
		{
		switch (flag)
		{
			case 0: //coefficient
				lastPos = i + 1; 
				for (; line[i] != ','; i++);
				coefBlock->data[termIndex % BLOCK_SIZE] = to_ZZ(line.substr(lastPos, i - lastPos).c_str());
				flag = 1;
				break;
			case 1: //exponent vector
				expBlock->data[termIndex % BLOCK_SIZE].SetLength(varCount);
				k = 0;
				for (i++; line[i] != ']'; i++)
				{
					if (line[i] != ' ')
					{
						lastPos = i;
						for (; line[i] != ',' && line[i] != ']'; i++);
						expBlock->data[termIndex % BLOCK_SIZE][k++] = to_ZZ(line.substr(lastPos, i - lastPos).c_str());
					}
				}
				termIndex++;
				flag = 0;
				break;
			default: //error
				cout << "Flag is " << flag << ", bailing." << endl;
				return;
		}
		}
	}

	monomials.termCount = termIndex;
	monomials.varCount = varCount;
	cout << "Loaded " << monomials.termCount << " monomials of dimension " << monomials.varCount << endl;
}

//Prints a nested list representation of our sum of monomials
//monomial sum: c_{1}*(x_{1}^e_{1}...x_{varCount}^e_{varCount}) + ...
//nested lists: [[c_{1}, [e_{1}, e_{2}, ..., e_{varCount}]], .. ]
string printMonomials(const monomialSum &myPoly)
{
	stringstream output (stringstream::in | stringstream::out);
	output << "[";
	eBlock* expTmp = myPoly.eHead; cBlock* coeffTmp = myPoly.cHead;
	int termCount = 0;
	do
	{
		for (int i = 0; i < BLOCK_SIZE && termCount < myPoly.termCount; i++)
		{
			output << "[" << coeffTmp->data[i] << ",[";
			for (int j = 0; j < myPoly.varCount; j++)
			{
				output << expTmp->data[i][j];
				if (j + 1 < myPoly.varCount)
				{ output << ","; }
			}
			output << "]]";
			if (termCount + 1 < myPoly.termCount)
			{ output << ","; }
			termCount++;
		}
		coeffTmp = coeffTmp->next; expTmp = expTmp->next;
	}
	while (coeffTmp != NULL);
	output << "]";
	return output.str();
}

//Deallocates space and nullifies internal pointers and counters
void destroyMonomials(monomialSum &myPoly)
{
	eBlock* expTmp = myPoly.eHead; cBlock* coeffTmp = myPoly.cHead;
	eBlock* oldExp = NULL; cBlock* oldCoeff = NULL;
	do
	{
		oldExp = expTmp; oldCoeff = coeffTmp;
		expTmp = expTmp->next;
		coeffTmp = coeffTmp->next;
		delete [] oldExp;
		delete [] oldCoeff;
	}
	while (coeffTmp != NULL);
	myPoly.eHead = NULL;
	myPoly.cHead = NULL;
	myPoly.termCount = myPoly.varCount = 0;
}

//Loads a string by parsing it as a sum of linear forms
//linear form: (c_{1} / d_{1}!)[(p_{1}*x_{1} + ... p_{varCount}*x_{varCount})^d_{1}] + ...
//nested list: [[c_{1}, [d_{1}, [p_{1}, p_{2}, ..., p_{varCount}]], .. ]
void loadLinForms(linFormSum &forms, const string line)
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

	//at least one record requires at least one block
	forms.lHead = (lBlock*) malloc (sizeof(lBlock));
	forms.lHead->next = NULL;
	forms.cHead = (cBlock*) malloc (sizeof(cBlock));
	forms.cHead->next = NULL;

	lBlock* formBlock = forms.lHead;
	cBlock* coefBlock = forms.cHead;
	formBlock->data = new vec_ZZ[BLOCK_SIZE];
	formBlock->degree = new ZZ[BLOCK_SIZE];
	coefBlock->data = new ZZ[BLOCK_SIZE];

	for (int i = 1; i < line.length() - 1; i++) //ignore outermost square brackets
	{
		if (line [i] == '[')
		{
		switch (flag)
		{
			case 0: //coefficient
				lastPos = i + 1; 
				for (; line[i] != ','; i++);
				coefBlock->data[termIndex % BLOCK_SIZE] = to_ZZ(line.substr(lastPos, i - lastPos).c_str());
				flag = 1;
				break;
			case 1: //degree
				lastPos = i + 1;
				for (; line[i] != ','; i++);
				formBlock->degree[termIndex % BLOCK_SIZE] = to_ZZ(line.substr(lastPos, i - lastPos).c_str());
				flag = 2;
				break;
			case 2: //coefficient vector
				formBlock->data[termIndex % BLOCK_SIZE].SetLength(varCount);
				k = 0;
				for (i++; line[i] != ']'; i++)
				{
					if (line[i] != ' ')
					{
						lastPos = i;
						for (; line[i] != ',' && line[i] != ']'; i++);
						formBlock->data[termIndex % BLOCK_SIZE][k++] = to_ZZ(line.substr(lastPos, i - lastPos).c_str());
					}
				}
				termIndex++;
				flag = 0;
				break;
			default: //error
				cout << "Flag is " << flag << ", bailing." << endl;
				return;
		}
		}
	}

	forms.termCount = termIndex;
	forms.varCount = varCount;
	cout << "Loaded " << forms.termCount << " terms of dimension " << forms.varCount << endl;
}

//Prints a nested list representation of our sum of linear forms
//linear form: (c_{1} / d_{1}!)[(p_{1}*x_{1} + ... p_{varCount}*x_{varCount})^d_{1}] + ...
//nested list: [[c_{1}, [d_{1}, [p_{1}, p_{2}, ..., p_{varCount}]], .. ]
string printLinForms(const linFormSum &myForm)
{
	//cout << "Printing form with " << myForm.termCount << " terms." << endl;
	stringstream output (stringstream::in | stringstream::out);
	output << "[";
	lBlock* formTmp = myForm.lHead; cBlock* coeffTmp = myForm.cHead;
	for (int i = 0; i < myForm.termCount; i++) //loop through each monomial
	{
		if (i > 0 && i % BLOCK_SIZE == 0)
		{
			formTmp = formTmp->next;
			coeffTmp = coeffTmp->next;
		}
		
		output << "[" << coeffTmp->data[i % BLOCK_SIZE] << ", [" << formTmp->degree[i % BLOCK_SIZE] << ", [";
		for (int j = 0; j < myForm.varCount; j++)
		{
			output << formTmp->data[i % BLOCK_SIZE][j];
			if (j + 1 < myForm.varCount)
			{ output << ", "; }
		}
		output << "]]]";
		if (i + 1 < myForm.termCount)
		{ output << ", "; }
	}
	output << "]";
	return output.str();
}

//Deallocates space and nullifies internal pointers and counters
void destroyLinForms(linFormSum &myPoly)
{
	lBlock* expTmp = myPoly.lHead; cBlock* coeffTmp = myPoly.cHead;
	lBlock* oldExp = NULL; cBlock* oldCoeff = NULL;
	int termCount = 0;
	do
	{
		oldExp = expTmp; oldCoeff = coeffTmp;
		expTmp = expTmp->next;
		coeffTmp = coeffTmp->next;
		delete [] oldExp;
		delete [] oldCoeff;
	}
	while (coeffTmp != NULL);
	myPoly.lHead = NULL;
	myPoly.cHead = NULL;
	myPoly.termCount = myPoly.varCount = 0;
}


//INPUT: monomial specified by myPoly.coefficientBlocks[mIndex / BLOCK_SIZE].data[mIndex % BLOCK_SIZE]
//	and myPoly.exponentBlocks[mIndex / BLOCK_SIZE].data[mIndex % BLOCK_SIZE]
//OUTPUT: lForm now also contains the linear decomposition of this monomial 
//	note: all linear form coefficients assumed to be divided by their respective |M|!, and the form is assumed to be of power M
void decompose(monomialSum &myPoly, linFormSum &lForm, int mIndex)
{
	//cout << "Decomposing monomial " << mIndex << " into linear forms: " << printPolynomial(myPoly) << endl;

	eBlock* expTmp = myPoly.eHead; cBlock* coeffTmp = myPoly.cHead;
	for (int i = 0; i < (mIndex / BLOCK_SIZE); i++) { expTmp = expTmp->next; coeffTmp = coeffTmp->next; }
	if (IsZero(expTmp->data[mIndex % BLOCK_SIZE])) //exponents are all 0, this is a constant term - linear form is already known
	{
		//search for a constant term linear form
		lBlock* linForm; cBlock* linCoeff;
		ZZ temp = coeffTmp->data[mIndex % BLOCK_SIZE];
		bool found;
		int myIndex;
		if (lForm.termCount > 0)
		{
			linForm = lForm.lHead; linCoeff = lForm.cHead; found = false;
			
			//check for compatible form here
			for (int i = 0; !found && i < lForm.termCount; i++)
			{
				if (i > 0 && i % BLOCK_SIZE == 0)
				{
					linForm = linForm->next, linCoeff = linCoeff->next;
				}
				
				if (IsZero(linForm->degree[i % BLOCK_SIZE]))
				{
					if (IsZero(linForm->data[i % BLOCK_SIZE]))
					{
						myIndex = i; found = true; break;
					}
				}
				
				if (i == lForm.termCount) { break; }
			}

			if (!found) //nothing found
			{
				//cout << "Didn't find compatible linear form, inserting term " << lForm.termCount + 1 << endl;
				if (lForm.termCount % BLOCK_SIZE == 0) //need to allocate a new block for coeffs and exponents
				{
					//cout << "Allocating new block" << endl;
					linCoeff->next = (cBlock*) malloc (sizeof(cBlock));
					linForm->next = (lBlock*) malloc (sizeof(lBlock));
					//cout << "alloc finished" << endl;
					linForm = linForm->next; linCoeff = linCoeff->next;
					linForm->next = NULL; linCoeff->next = NULL;
					//cout << "set up pointers" << endl;
					linForm->data = new vec_ZZ[BLOCK_SIZE]; linForm->degree = new ZZ[BLOCK_SIZE]; linCoeff->data = new ZZ[BLOCK_SIZE];
				}
				linForm->data[lForm.termCount % BLOCK_SIZE].SetLength(myPoly.varCount);
				//cout << "set length" << endl;
				VectorCopy(linForm->data[lForm.termCount % BLOCK_SIZE], expTmp->data[mIndex % BLOCK_SIZE], myPoly.varCount);
				//cout << "copied vector" << endl;
				linForm->degree[lForm.termCount % BLOCK_SIZE] = to_ZZ(0);
				//cout << "set degree" << endl;
				linCoeff->data[lForm.termCount % BLOCK_SIZE] = temp;
				//cout << "Added following term: [" << temp << ", [" << totalDegree << ", " << linForm->data[lForm.termCount % BLOCK_SIZE] << "]]" << endl;
				lForm.termCount++;
			}
			else //found this linear form
			{
				//cout << "Found compatible linear form at " << myIndex << endl;
				linCoeff->data[myIndex % BLOCK_SIZE] += temp;
				//cout << "Modified term: [" << linCoeff->data[myIndex % BLOCK_SIZE] << ", [" << linForm->degree[myIndex  % BLOCK_SIZE] << ", " << linForm->data[myIndex % BLOCK_SIZE] << "]]" << endl;
			}
		}
		else //no constant linear form, creating one
		{
			//cout << "NULL head" << endl;
			lForm.cHead = (cBlock*) malloc (sizeof(cBlock));
			lForm.lHead = (lBlock*) malloc (sizeof(lBlock));
			linForm = lForm.lHead; linCoeff = lForm.cHead;
			linForm->next = NULL; linCoeff->next = NULL;
			linForm->data = new vec_ZZ[BLOCK_SIZE];  linForm->degree = new ZZ[BLOCK_SIZE]; linCoeff->data = new ZZ[BLOCK_SIZE];
			linForm->data[lForm.termCount % BLOCK_SIZE].SetLength(myPoly.varCount);
			VectorCopy(linForm->data[lForm.termCount % BLOCK_SIZE], expTmp->data[mIndex % BLOCK_SIZE], myPoly.varCount);
			linForm->degree[lForm.termCount % BLOCK_SIZE] = 0;
			linCoeff->data[lForm.termCount % BLOCK_SIZE] = temp;
			//cout << "Added following term: [" << temp << ", [" << totalDegree << ", " << linForm->data[lForm.termCount % BLOCK_SIZE] << "]]" << endl;
			lForm.termCount++;
		}
		return;
	}
	ZZ formsCount = expTmp->data[mIndex % BLOCK_SIZE][0] + 1;
	ZZ totalDegree = expTmp->data[mIndex % BLOCK_SIZE][0];
	for (int i = 1; i < myPoly.varCount; i++)
	{
		formsCount *= expTmp->data[mIndex % BLOCK_SIZE][i] + 1;
		totalDegree += expTmp->data[mIndex % BLOCK_SIZE][i];
	}
	formsCount--;
	//cout << "At most " << formsCount << " linear forms will be required for this decomposition." << endl;
	//cout << "Total degree is " << totalDegree << endl;

	vec_ZZ p, counter;
	ZZ temp;
	ZZ g;
	bool found;
	int myIndex;
	counter.SetLength(myPoly.varCount);
	p.SetLength(myPoly.varCount);
	for (ZZ i = to_ZZ(1); i <= formsCount; i++)
	{
		//cout << "i is " << i << endl;
		counter[0] += to_ZZ(1);
		for (int j = 0; counter[j] > expTmp->data[mIndex % BLOCK_SIZE][j]; j++)
		{
			counter[j] = to_ZZ(0);
			counter[j+1] += to_ZZ(1);
		}
		//cout << "counter is: " << counter << endl;
		
		//find gcd of all elements
		g = counter[0];
		ZZ parity = totalDegree - counter[0];
		p[0] = counter[0];
		if (myPoly.varCount > 1)
		{
			for (int k = 1; k < myPoly.varCount; k++)
			{
				p[k] = counter[k];
				g = GCD (g, p[k]);
				parity -= p[k];
			}
		}
		//cout << "GCD is " << g << endl;
		
		//calculate coefficient
		temp = coeffTmp->data[mIndex % BLOCK_SIZE];
		if (IsOdd(parity)) { temp *= -1; } // -1 ^ [|M| - (p[0] + p[1] + ... p[n])]
		for (int j = 0; j < myPoly.varCount; j++) //calculate binomial coefficients
		{
			for (ZZ k = expTmp->data[mIndex % BLOCK_SIZE][j]; k > p[j]; k--) { temp *= k; } // M[j]! / p[j]!
			for (ZZ k = expTmp->data[mIndex % BLOCK_SIZE][j] - p[j]; k > 1; k--) { temp /= k; } // 1 / (M[j] - p[j])!
		}
		
		if (!IsOne(g))
		{
			for (int k = 0; k < myPoly.varCount; k++)
			{
				p[k] /= g;
			}
			temp *= Power(g, totalDegree);
		}
		//cout << "Normalized p is: " << p << endl;
		//cout << "coefficient is " << temp << endl;
		// is there a vector in lForm's exponent block equal to p?
		// 	yes: add our coefficient to the existing one, we're done
		//	no : allocate space for a new vec_ZZ, insert ours with our coefficient, done
		lBlock* linForm; cBlock* linCoeff;
		if (lForm.termCount > 0)
		{
			//cout << lForm.termCount << " linear forms present" << endl;
			linForm = lForm.lHead; linCoeff = lForm.cHead; found = false;
			
			//check for compatible form here
			for (int i = 0; !found && i < lForm.termCount; i++)
			{
				if (i > 0 && i % BLOCK_SIZE == 0)
				{
					linForm = linForm->next, linCoeff = linCoeff->next;
				}
				
				if (linForm->degree[i % BLOCK_SIZE] == totalDegree)
				{
					if (linForm->data[i % BLOCK_SIZE] == p)
					{
						myIndex = i; found = true; break;
					}
				}
				
				if (i == lForm.termCount) { break; }
			}

			if (!found) //nothing found
			{
				//cout << "Didn't find compatible linear form, inserting term " << lForm.termCount + 1 << endl;
				if (lForm.termCount % BLOCK_SIZE == 0) //need to allocate a new block for coeffs and exponents
				{
					//cout << "Allocating new block" << endl;
					linCoeff->next = (cBlock*) malloc (sizeof(cBlock));
					linForm->next = (lBlock*) malloc (sizeof(lBlock));
					//cout << "alloc finished" << endl;
					linForm = linForm->next; linCoeff = linCoeff->next;
					linForm->next = NULL; linCoeff->next = NULL;
					//cout << "set up pointers" << endl;
					linForm->data = new vec_ZZ[BLOCK_SIZE]; linForm->degree = new ZZ[BLOCK_SIZE]; linCoeff->data = new ZZ[BLOCK_SIZE];
				}
				linForm->data[lForm.termCount % BLOCK_SIZE].SetLength(myPoly.varCount);
				//cout << "set length" << endl;
				VectorCopy(linForm->data[lForm.termCount % BLOCK_SIZE], p, myPoly.varCount);
				//cout << "copied vector" << endl;
				linForm->degree[lForm.termCount % BLOCK_SIZE] = totalDegree;
				//cout << "set degree" << endl;
				linCoeff->data[lForm.termCount % BLOCK_SIZE] = temp;
				//cout << "Added following term: [" << temp << ", [" << totalDegree << ", " << linForm->data[lForm.termCount % BLOCK_SIZE] << "]]" << endl;
				lForm.termCount++;
			}
			else //found this linear form
			{
				//cout << "Found compatible linear form at " << myIndex << endl;
				linCoeff->data[myIndex % BLOCK_SIZE] += temp;
				//cout << "Modified term: [" << linCoeff->data[myIndex % BLOCK_SIZE] << ", [" << linForm->degree[myIndex  % BLOCK_SIZE] << ", " << linForm->data[myIndex % BLOCK_SIZE] << "]]" << endl;
			}
		}
		else
		{
			//cout << "NULL head" << endl;
			lForm.cHead = (cBlock*) malloc (sizeof(cBlock));
			lForm.lHead = (lBlock*) malloc (sizeof(lBlock));
			linForm = lForm.lHead; linCoeff = lForm.cHead;
			linForm->next = NULL; linCoeff->next = NULL;
			linForm->data = new vec_ZZ[BLOCK_SIZE];  linForm->degree = new ZZ[BLOCK_SIZE]; linCoeff->data = new ZZ[BLOCK_SIZE];
			linForm->data[lForm.termCount % BLOCK_SIZE].SetLength(myPoly.varCount);
			VectorCopy(linForm->data[lForm.termCount % BLOCK_SIZE], p, myPoly.varCount);
			linForm->degree[lForm.termCount % BLOCK_SIZE] = totalDegree;
			linCoeff->data[lForm.termCount % BLOCK_SIZE] = temp;
			//cout << "Added following term: [" << temp << ", [" << totalDegree << ", " << linForm->data[lForm.termCount % BLOCK_SIZE] << "]]" << endl;
			lForm.termCount++;
		}
	}
	p.kill();
	counter.kill();
}

ZZ Power(const ZZ&a, const ZZ& e)
{
	if (e == 0) return to_ZZ(1);

	long k = NumBits(e);
     
	ZZ res;
	res = 1;
     
	for (long i = k-1; i >= 0; i--) {
	   res = (res*res);
	   if (bit(e, i) == 1) res = (res*a);
	}
	return res;
}

/*
string printMapleExpr(const linFormSum &myForm)
{
	stringstream output (stringstream::in | stringstream::out);
	lBlock* formTmp = myForm.lHead; cBlock* coeffTmp = myForm.cHead;
	cout << "there are " << myForm.termCount << " terms." << endl;
	for (int i = 0; i < myForm.termCount; i++)
	{
		if (i > 0 && i % BLOCK_SIZE == 0)
		{
			formTmp = formTmp->next;
			coeffTmp = coeffTmp->next;
		}
		
		output << "(" << coeffTmp->data[i % BLOCK_SIZE] << "/factorial(" << formTmp->degree[i % BLOCK_SIZE] << "))*("; // divide coefficient by |M|!
		for (int j = 0; j < myForm.varCount; j++)
		{
			output << formTmp->data[i % BLOCK_SIZE][j] << "*x[" << j + 1 << "]";
			if (j + 1 < myForm.varCount)
			{ output << " + "; }
		}
		output << ")^" << formTmp->degree[i % BLOCK_SIZE];
		if (i + 1 < myForm.termCount)
		{ output << " + "; }
	}
	return output.str();
}
*/
