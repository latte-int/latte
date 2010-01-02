#include "PolyRep.h"
#include <stdio.h>
#include <sstream>

void loadPolynomial(polynomial &myPoly, const string line)
{
	//cout << "Loading polynomial `" << line << "'" << endl;
	int varCount = 1;
	string temp = line.substr(line.find(",[") + 2, line.find("]],") - (line.find(",[") + 2));
	for (int i = 0; i < temp.length(); i++)
	{ varCount += (temp.at(i) == ','); }
	int termCount = 1; //this is 1 + occurences of "]," in line
	for (int i = 2; i < line.length() - 1; i++) //first 2 characters are "[["
	{
		if (line.at(i) == ']' && line.at(i+1) == ',') { termCount++; }
	}
	myPoly.eHead = (eBlock*) malloc (sizeof(eBlock));
	myPoly.eHead->next = NULL;
	myPoly.cHead = (cBlock*) malloc (sizeof(cBlock));
	myPoly.cHead->next = NULL;
	//still need to ensure that terms are sorted by exponent?

	int index = 1; char *exps, *myTok;
	bool exponents = false;
	int termIndex, varIndex;
	eBlock* expTmp = myPoly.eHead; cBlock* coeffTmp = myPoly.cHead;
	expTmp->data = new vec_ZZ[BLOCK_SIZE]; coeffTmp->data = new ZZ[BLOCK_SIZE];
	termIndex = 0;
	while (index < line.length() - 1)
	{
		if (exponents)
		{
			temp = line.substr(index + 1, line.find("]]", index) - index - 1);
			index += temp.length() + 4;
			exps = (char*) malloc (sizeof(char) * temp.length());
			strcpy(exps, temp.c_str());
			myTok = strtok(exps, ",");
			expTmp->data[termIndex % BLOCK_SIZE].SetLength(varCount);
			varIndex = 0;
			while (myTok != NULL)
			{
				expTmp->data[termIndex % BLOCK_SIZE][varIndex++] = to_ZZ(myTok);
				myTok = strtok (NULL, ",");
			}
			exponents = false;
			termIndex++;
			free (exps);
		}
		else
		{
			if (termIndex > 0 && termIndex % BLOCK_SIZE == 0)
			{
				expTmp->next = (eBlock*) malloc (sizeof(eBlock));
				coeffTmp->next = (cBlock*) malloc (sizeof(cBlock));
				expTmp = expTmp->next; coeffTmp = coeffTmp->next;
				expTmp->next = NULL; coeffTmp->next = NULL;
				expTmp->data = new vec_ZZ[BLOCK_SIZE];
				coeffTmp->data = new ZZ[BLOCK_SIZE];
			}
			temp = line.substr(index + 1, line.find(",", index) - index - 1);
			index += temp.length() + 2;
			coeffTmp->data[termIndex % BLOCK_SIZE] = to_ZZ(temp.c_str());
			exponents = true;
		}
	}
	myPoly.termCount = termCount;
	myPoly.varCount = varCount;
}

void decompose(polynomial &myPoly, linearPoly &lForm, int mIndex)
{
	//INPUT: monomial specified by myPoly.coefficientBlocks[mIndex / BLOCK_SIZE].data[mIndex % BLOCK_SIZE]
	//	 			and myPoly.exponentBlocks[mIndex / BLOCK_SIZE].data[mIndex % BLOCK_SIZE]
	//OUTPUT: lForm now also contains the linear decomposition of this monomial (all linear form coefficients assumed to be divided by their respective |M|!, and the form is assumed to be of power M)

	//cout << "Decomposing monomial " << mIndex << " into linear forms: " << printPolynomial(myPoly) << endl;

	eBlock* expTmp = myPoly.eHead; cBlock* coeffTmp = myPoly.cHead;
	for (int i = 0; i < (mIndex / BLOCK_SIZE); i++) { expTmp = expTmp->next; coeffTmp = coeffTmp->next; }
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

string printPolynomial(const polynomial &myPoly)
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

string printForm(const linearPoly &myForm)
{
	cout << "Printing form with " << myForm.termCount << " terms." << endl;
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
		
		for (int k = 0; k < BLOCK_SIZE && i < myForm.termCount; k++, i++)
		{
			output << "[" << coeffTmp->data[i] << ",[" << formTmp->degree[i] << ",[";
			for (int j = 0; j < myForm.varCount; j++)
			{
				output << formTmp->data[k][j];
				if (j + 1 < myForm.varCount)
				{ output << ","; }
			}
			output << "]]]";
			if (i	 + 1 < myForm.termCount)
			{ output << ","; }
		}
	}
	/*
	do
	{
		for (int i = 0; i < BLOCK_SIZE && termCount < myForm.termCount; i++)
		{
			output << "[" << coeffTmp->data[i] << ",[" << formTmp->degree[i] << ",[";
			for (int j = 0; j < myForm.varCount; j++)
			{
				output << formTmp->data[i][j];
				if (j + 1 < myForm.varCount)
				{ output << ","; }
			}
			output << "]]]";
			if (termCount + 1 < myForm.termCount)
			{ output << ","; }
			termCount++;
		}
		coeffTmp = coeffTmp->next; formTmp = formTmp->next;
	}
	while (coeffTmp != NULL);*/
	output << "]";
	return output.str();
}

string printMapleForm(const linearPoly &myForm)
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

void destroyForm(linearPoly &myPoly)
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

void destroyPolynomial(polynomial &myPoly)
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
