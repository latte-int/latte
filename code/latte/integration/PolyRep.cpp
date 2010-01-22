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
	cout << varCount << " variables." << endl;

	//at least one record requires at least one block
	monomials.eHead = (eBlock*) malloc (sizeof(eBlock));
	monomials.eHead->next = NULL;
	monomials.cHead = (cBlock*) malloc (sizeof(cBlock));
	monomials.cHead->next = NULL;

	eBlock* expBlock = monomials.eHead;
	cBlock* coefBlock = monomials.cHead;
	expBlock->data = new int[varCount * BLOCK_SIZE];
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
				k = varCount * (termIndex % BLOCK_SIZE); //0 for the 1st term, varCount for the 2nd, etc.
				for (i++; line[i] != ']'; i++)
				{
					if (line[i] != ' ')
					{
						lastPos = i;
						for (; line[i] != ',' && line[i] != ']'; i++);
						expBlock->data[k++] = atoi(line.substr(lastPos, i - lastPos).c_str());
						
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
	//cout << "Loaded " << monomials.termCount << " monomials of dimension " << monomials.varCount << endl;
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
			for (int j = (i * myPoly.varCount); j < ((i + 1) * myPoly.varCount); j++)
			{
				output << expTmp->data[j];
				if (j + 1 < ((i + 1) * myPoly.varCount))
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
				formBlock->degree[termIndex % BLOCK_SIZE] = atoi(line.substr(lastPos, i - lastPos).c_str());
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
}

//Prints a nested list representation of our sum of linear forms
//linear form: (c_{1} / d_{1}!)[(p_{1}*x_{1} + ... p_{varCount}*x_{varCount})^d_{1}] + ...
//nested list: [[c_{1}, [d_{1}, [p_{1}, p_{2}, ..., p_{varCount}]], .. ]
string printLinForms(const linFormSum &myForm)
{
	stringstream output (stringstream::in | stringstream::out);
	output << "[";
	lBlock* formTmp = myForm.lHead; cBlock* coeffTmp = myForm.cHead;
	for (int i = 0; i < myForm.termCount; i++)
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
	eBlock* expTmp = myPoly.eHead; cBlock* coeffTmp = myPoly.cHead;
	for (int i = 0; i < (mIndex / BLOCK_SIZE); i++) { expTmp = expTmp->next; coeffTmp = coeffTmp->next; }
	bool constantTerm = true;
	for (int i = (mIndex * myPoly.varCount); i < ((mIndex + 1) * myPoly.varCount); i++) { if (expTmp->data[i] != 0) { constantTerm = false; break; }}
	if (constantTerm) //exponents are all 0, this is a constant term - linear form is already known
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
				
				if (linForm->degree[i % BLOCK_SIZE] == 0)
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
				if (lForm.termCount % BLOCK_SIZE == 0) //need to allocate a new block for coeffs and exponents
				{
					linCoeff->next = (cBlock*) malloc (sizeof(cBlock));
					linForm->next = (lBlock*) malloc (sizeof(lBlock));
					linForm = linForm->next; linCoeff = linCoeff->next;
					linForm->next = NULL; linCoeff->next = NULL;
					linForm->data = new vec_ZZ[BLOCK_SIZE];
					linCoeff->data = new ZZ[BLOCK_SIZE];
				}
				linForm->data[lForm.termCount % BLOCK_SIZE].SetLength(myPoly.varCount); //should be initialized to default ZZ, which is 0
				linForm->degree[lForm.termCount % BLOCK_SIZE] = 0;
				linCoeff->data[lForm.termCount % BLOCK_SIZE] = temp;
				lForm.termCount++;
			}
			else //found this linear form
			{
				linCoeff->data[myIndex % BLOCK_SIZE] += temp;
			}
		}
		else //no constant linear form, creating one
		{
			lForm.cHead = (cBlock*) malloc (sizeof(cBlock));
			lForm.lHead = (lBlock*) malloc (sizeof(lBlock));
			linForm = lForm.lHead; linCoeff = lForm.cHead;
			linForm->next = NULL; linCoeff->next = NULL;
			linForm->data[lForm.termCount % BLOCK_SIZE].SetLength(myPoly.varCount); //should be initialized to default ZZ, which is 0
			linCoeff->data = new ZZ[BLOCK_SIZE];
			linForm->degree[lForm.termCount % BLOCK_SIZE] = 0;
			linCoeff->data[lForm.termCount % BLOCK_SIZE] = temp;
			lForm.termCount++;
		}
		return;
	}
	ZZ formsCount = to_ZZ(expTmp->data[(mIndex % BLOCK_SIZE) * myPoly.varCount] + 1);
	int totalDegree = expTmp->data[(mIndex % BLOCK_SIZE) * myPoly.varCount];
	for (int i = 1; i < myPoly.varCount; i++)
	{
		formsCount *= expTmp->data[(mIndex % BLOCK_SIZE) * myPoly.varCount + i] + 1;
		totalDegree += expTmp->data[(mIndex % BLOCK_SIZE) * myPoly.varCount + i];
	}
	formsCount--;
	//cout << "At most " << formsCount << " linear forms will be required for this decomposition." << endl;
	//cout << "Total degree is " << totalDegree << endl;

	int* p = new int[myPoly.varCount];
	int* counter = new int[myPoly.varCount];
	ZZ* binomCoeffs = new ZZ[myPoly.varCount]; //for calculating the product of binomial coefficients for each linear form

	ZZ temp;
	int g;
	bool found;
	int myIndex;
	for (int i = 0; i < myPoly.varCount; i++) { counter[i] = 0; binomCoeffs[i] = to_ZZ(1); }
	for (ZZ i = to_ZZ(1); i <= formsCount; i++)
	{
		//cout << "i is " << i << endl;
		counter[0] += 1;
		for (myIndex = 0; counter[myIndex] > expTmp->data[(mIndex % BLOCK_SIZE) * myPoly.varCount + myIndex]; myIndex++)
		{
			counter[myIndex] = 0;
			binomCoeffs[myIndex] = to_ZZ(1);
			counter[myIndex+1] += 1;
		}
		binomCoeffs[myIndex] *= expTmp->data[(mIndex % BLOCK_SIZE) * myPoly.varCount + myIndex] - counter[myIndex] + 1; // [n choose k] = [n choose (k - 1) ] * (n - k + 1)/k
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
		temp *= coeffTmp->data[mIndex % BLOCK_SIZE];
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
					found = true;
					for (int j = 0; j < myPoly.varCount; j++)
					{
						if (linForm->data[i % BLOCK_SIZE][j] != p[j])
						{ found = false; break; }
					}
					if (found)
					{
						myIndex = i; break;
					}
				}
				
				if (i == lForm.termCount) { break; }
			}

			if (!found) //nothing found
			{
				if (lForm.termCount % BLOCK_SIZE == 0) //need to allocate a new block for coeffs and exponents
				{
					linCoeff->next = (cBlock*) malloc (sizeof(cBlock));
					linForm->next = (lBlock*) malloc (sizeof(lBlock));
					linForm = linForm->next; linCoeff = linCoeff->next;
					linForm->next = NULL; linCoeff->next = NULL;
					linForm->data = new vec_ZZ[BLOCK_SIZE];
					linCoeff->data = new ZZ[BLOCK_SIZE];
				}
				linForm->data[lForm.termCount % BLOCK_SIZE].SetLength(myPoly.varCount);
				for (int j = 0; j < myPoly.varCount; j++)
				{
					linForm->data[lForm.termCount % BLOCK_SIZE][j] = p[j];
				}
				linForm->degree[lForm.termCount % BLOCK_SIZE] = totalDegree;
				linCoeff->data[lForm.termCount % BLOCK_SIZE] = temp;
				lForm.termCount++;
			}
			else //found this linear form
			{
				linCoeff->data[myIndex % BLOCK_SIZE] += temp;
			}
		}
		else
		{
			lForm.cHead = (cBlock*) malloc (sizeof(cBlock));
			lForm.lHead = (lBlock*) malloc (sizeof(lBlock));
			linForm = lForm.lHead; linCoeff = lForm.cHead;
			linForm->next = NULL; linCoeff->next = NULL;
			linForm->data = new vec_ZZ[BLOCK_SIZE];
			linCoeff->data = new ZZ[BLOCK_SIZE];
			linForm->data[lForm.termCount % BLOCK_SIZE].SetLength(myPoly.varCount);
			for (int j = 0; j < myPoly.varCount; j++)
			{
				linForm->data[lForm.termCount % BLOCK_SIZE][j] = p[j];
			}
			linForm->degree[lForm.termCount % BLOCK_SIZE] = totalDegree;
			linCoeff->data[lForm.termCount % BLOCK_SIZE] = temp;
			lForm.termCount++;
		}
	}
	delete [] p;
	delete [] counter;
	delete [] binomCoeffs
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
