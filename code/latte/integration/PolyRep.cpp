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
			}
			temp = line.substr(index + 1, line.find(",", index) - index - 1);
			index += temp.length() + 2;
			coeffTmp->data[termIndex % BLOCK_SIZE] = to_ZZ(temp.c_str());
			exponents = true;
		}
	}
	myPoly.termCount = termCount;
	myPoly.varCount = varCount;
	myPoly.linearForm = false;

	polynomial lForm;
	lForm.termCount = 0;
	lForm.varCount = myPoly.varCount;
	//need to allocate space for exponent blocks / coeff blocks?
	for (int i = 0; i < myPoly.termCount; i++)
	{
		decompose(myPoly, lForm, i);
	}
	lForm.linearForm = true;
}

void decompose(polynomial &myPoly, polynomial &lForm, int mIndex)
{
	//INPUT: monomial specified by myPoly.coefficientBlocks[mIndex / BLOCK_SIZE].data[mIndex % BLOCK_SIZE]
	//	 			and myPoly.exponentBlocks[mIndex / BLOCK_SIZE].data[mIndex % BLOCK_SIZE]
	//OUTPUT: lForm now also contains the linear decomposition of this monomial

	//how to deal with ( 1 / |M|! ) ? maybe store inverse of coefficient?
	cout << "Decomposing monomial " << mIndex << " into linear forms: " << printPolynomial(myPoly) << endl;

	eBlock* expTmp = myPoly.eHead; cBlock* coeffTmp = myPoly.cHead;
	for (int i = 0; i < (mIndex / BLOCK_SIZE); i++) { expTmp = expTmp->next; coeffTmp = coeffTmp->next; }
	ZZ formsCount = expTmp->data[mIndex % BLOCK_SIZE][0] + 1;
	for (int i = 1; i < myPoly.varCount; i++)
	{
		formsCount *= expTmp->data[mIndex % BLOCK_SIZE][i] + 1;
	}
	formsCount--;
	cout << "At most " << formsCount << " linear forms will be required for this decomposition." << endl;

	vec_ZZ p;
	ZZ temp;
	int myIndex;
	bool found;
	p.SetLength(myPoly.varCount);
	for (ZZ i = to_ZZ(1); i <= formsCount; i++)
	{
		p[myPoly.varCount - 1] = i % (expTmp->data[mIndex % BLOCK_SIZE][myPoly.varCount- 1] + 1);
		temp = expTmp->data[mIndex % BLOCK_SIZE][myPoly.varCount - 1] + 1;
		for (int j = myPoly.varCount - 2; j >= 0; j--)
		{
			p[j] = i / temp;
			temp *= (expTmp->data[mIndex % BLOCK_SIZE][j] + 1);
		}
		cout << "p is: " << p << endl;
		//find gcd of all elements?
		temp = p[0];
		if (myPoly.varCount > 1)
		{
			for (int k = 1; k < myPoly.varCount; k++)
			{
				temp = GCD (temp, p[k]);
			}
		}
		cout << "GCD is " << temp << endl;
		if (!IsOne(temp))
		{
			for (int k = 0; k < myPoly.varCount; k++)
			{
				p[k] /= temp;
			}
		}
		cout << "Normalized p is : " << p << endl;
		// if only one p is non-zero, should simplify it by bringing p^|M| outside the form
		// now, calculate and store the coefficient (multiplied by the GCD stored in temp) in front of the simplified linear form (assign it to temp)
		// is there a vector in lForm's exponent block equal to p?
		// 	yes: add our coefficient to the existing one, we're done
		//	no : allocate space for a new vec_ZZ, insert ours with our coefficient, done
		eBlock* linExp; cBlock* linCoeff;
		if (lForm.cHead)
		{
			cout << "Head not null" << endl;
			linExp = lForm.eHead;
			linCoeff = lForm.cHead;

			found = false;
			do
			{
				for (myIndex = 0; myIndex < BLOCK_SIZE && linExp->data[myIndex].length() > 0 && linExp->data[myIndex] != p; myIndex++);
				if (myIndex < BLOCK_SIZE && linExp->data[myIndex].length() > 0)
				{
					found = true;
				}
				else
				{
					if (!linCoeff->next) { break; }
					linCoeff = linCoeff->next; linExp = linExp->next;
				}
			}
			while (!found);
			if (!found) //nothing found
			{
				cout << "Didn't find compatible linear form, inserting term " << lForm.termCount << endl;
				if (lForm.termCount % BLOCK_SIZE == 0) //need to allocate a new block for coeffs and exponents
				{
					cout << "Allocating new block" << endl;
					linCoeff->next = (cBlock*) malloc (sizeof(cBlock));
					linExp->next = (eBlock*) malloc (sizeof(eBlock));
					linExp = linExp->next; linCoeff = linCoeff->next;
					linExp->next = NULL; linCoeff->next = NULL;
				}
				linExp->data[lForm.termCount % BLOCK_SIZE].SetLength(myPoly.termCount);
				cout << "Set Length" << endl;
				VectorCopy(linExp->data[lForm.termCount % BLOCK_SIZE], p, myPoly.termCount);
				cout << "Copied vector" << endl;
				linCoeff->data[lForm.termCount % BLOCK_SIZE] = temp;
				lForm.termCount++;
				cout << "Added term" << endl;
			}
			else //found this linear form
			{
				cout << "Found compatible linear form" << endl;
				linCoeff->data[myIndex % BLOCK_SIZE] += temp;
			}
		}
		else
		{
			cout << "NULL head" << endl;
			lForm.cHead = (cBlock*) malloc (sizeof(cBlock));
			lForm.eHead = (eBlock*) malloc (sizeof(eBlock));
			linExp = lForm.eHead; linCoeff = lForm.cHead;
			linExp->next = NULL; linCoeff->next = NULL;
			VectorCopy(linExp->data[lForm.termCount % BLOCK_SIZE], p, myPoly.termCount);
			linCoeff->data[lForm.termCount % BLOCK_SIZE] = temp;
			lForm.termCount++;
if (!linExp) { cout << "exps null" << endl; }
		}
	}
}

string printPolynomial(const polynomial &myPoly)
{
	stringstream output (stringstream::in | stringstream::out);
	output << "[";
	if (myPoly.linearForm)
	{

	}
	else
	{
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
	}
	output << "]";
	return output.str();
}

void destroyPolynomial(polynomial &myPoly)
{
	if (myPoly.linearForm)
	{
		
	}
	else
	{
		eBlock* expTmp = myPoly.eHead; cBlock* coeffTmp = myPoly.cHead;
		eBlock* oldExp = NULL; cBlock* oldCoeff = NULL;
		int termCount = 0;
		do
		{
			for (int i = 0; i < BLOCK_SIZE && termCount < myPoly.termCount; i++)
			{
				expTmp->data[i].kill();
				coeffTmp->data[i].kill();
				termCount++;
			}
			oldExp = expTmp; oldCoeff = coeffTmp;
			expTmp = expTmp->next;
			coeffTmp = coeffTmp->next;
			free (oldExp);
			free (oldCoeff);
		}
		while (coeffTmp != NULL);
	}
	myPoly.eHead = NULL;
	myPoly.cHead = NULL;
	myPoly.termCount = myPoly.varCount = 0;
}
