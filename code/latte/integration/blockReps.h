#ifndef REPRESENTATIONS_H
#define REPRESENTATIONS_H

#include <NTL/vec_ZZ.h>
#include <NTL/ZZ.h>

NTL_CLIENT

#define BLOCK_SIZE 64

//block-based

//this defines the data structures used
struct eBlock //used to store exponent vector
{
	eBlock* next;
	int* data; //represents 2D array of ints, storing an exponent vector for each monomial
	//in reality, it's 1D array of length dimension*BLOCK_SIZE
};

struct lBlock //used to store linear forms (coefficient vector and degree)
{
	lBlock* next;
	vec_ZZ* data; //represents 2D array of ints, storing a coefficient vector for each linear form
	//in reality, it's 1D array of length dimension*BLOCK_SIZE
	int degree[BLOCK_SIZE]; //total degree of each linear form
};

template <class T> struct cBlock
{
	cBlock* next;
	T* data; //this stores the monomial coefficient for each monomial
};

struct _monomialSum
{
	int termCount, varCount;
	eBlock* eHead; //variable exponents
	cBlock<ZZ>* cHead; //monomial coefficients over ZZ
};

struct _linFormSum
{
	int termCount, varCount;
	lBlock* lHead; //linear forms
	cBlock<ZZ>* cHead; //linear form coefficients over ZZ
};

// Attempts to find a monomial in _monomialSum with same exponents as those passed in
// if found, the monomial coefficient in _monomialSum is multiplied by the one passed in
// if not found, a new monomial term is added and _monomialSum's termCount is incremented
// if _monomialSum is empty, this sets formSum's eHead and cHead variables
template <class T>
void _insertMonomial(const T& coefficient, int* exponents, _monomialSum& monomials)
{
	bool found;
	int myIndex;
	eBlock* myExps; cBlock<T>* myCoeffs;
	
	if (monomials.termCount > 0)
	{
		myExps = monomials.eHead; myCoeffs = monomials.cHead; found = false;
		
		//check for compatible monomial here 
		for (int i = 0; !found && i < monomials.termCount; i++)
		{
			if (i > 0 && i % BLOCK_SIZE == 0)
			{
				myExps = myExps->next, myCoeffs = myCoeffs->next;
			}
			
			found = true;
			for (int j = 0; j < monomials.varCount; j++)
			{
				if (myExps->data[(i % BLOCK_SIZE)*monomials.varCount + j] != exponents[j])
				{ found = false; break; }
			}
			if (found)
			{
				myIndex = i; break;
			}
			
			if (i == monomials.termCount) { break; }
		}
		
		if (!found) //nothing found
		{
			if (monomials.termCount % BLOCK_SIZE == 0) //need to allocate a new block for coeffs and exponents
			{
				myCoeffs->next = (cBlock<T>*) malloc (sizeof(cBlock<T>));
				myExps->next = (eBlock*) malloc (sizeof(eBlock));
				myExps = myExps->next; myCoeffs = myCoeffs->next;
				myExps->next = NULL; myCoeffs->next = NULL;
				myExps->data = new int[monomials.varCount * BLOCK_SIZE];
				myCoeffs->data = new T[BLOCK_SIZE];
			}
			for (int j = 0; j < monomials.varCount; j++)
			{
				myExps->data[(monomials.termCount % BLOCK_SIZE)*monomials.varCount + j] = exponents[j];
			}
			myCoeffs->data[monomials.termCount % BLOCK_SIZE] = coefficient;
			monomials.termCount++;
		}
		else //found this linear form
		{
			myCoeffs->data[myIndex % BLOCK_SIZE] += coefficient;
		}
	}
	else
	{
		monomials.cHead = (cBlock<T>*) malloc (sizeof(cBlock<T>));
		monomials.eHead = (eBlock*) malloc (sizeof(eBlock));
		myExps = monomials.eHead; myCoeffs = monomials.cHead;
		myExps->next = NULL; myCoeffs->next = NULL;
		myExps->data = new int[monomials.varCount * BLOCK_SIZE];
		myCoeffs->data = new T[BLOCK_SIZE];
		for (int j = 0; j < monomials.varCount; j++)
		{
			myExps->data[(monomials.termCount % BLOCK_SIZE)*monomials.varCount + j] = exponents[j];
		}
		myCoeffs->data[monomials.termCount % BLOCK_SIZE] = coefficient;
		monomials.termCount++;
	}
}

// Attempts to find a linear form in formSum with same degree and coefficients as those passed in
// if found, the linear form coefficient in formSum is incremented by coef
// if not found, a new linear form term is added and formSum's termCount is incremented
// if formSum is empty, this sets formSum's lHead and cHead variables
template <class T>
void _insertLinForm(const T& coef, int degree, const vec_ZZ& coeffs, _linFormSum& formSum)
{
	bool found;
	int myIndex;
	lBlock* linForm; cBlock<T>* linCoeff;
	
	if (formSum.termCount > 0)
	{
		//cout << lForm.termCount << " linear forms present" << endl;
		linForm = formSum.lHead; linCoeff = formSum.cHead; found = false;
		
		//check for compatible form here 
		for (int i = 0; !found && i < formSum.termCount; i++)
		{
			if (i > 0 && i % BLOCK_SIZE == 0)
			{
				linForm = linForm->next, linCoeff = linCoeff->next;
			}
			
			if (linForm->degree[i % BLOCK_SIZE] == degree)
			{
				//found = true;
				if (linForm->data[i % BLOCK_SIZE] == coeffs)
				/*for (int j = 0; j < formSum.varCount; j++)
				{
					if (linForm->data[i % BLOCK_SIZE][j] != coeffs[j])
					{ found = false; break; }
				}
				if (found)*/
				{
					found = true; myIndex = i; break;
				}
			}
			
			if (i == formSum.termCount) { break; }
		}
		
		if (!found) //nothing found
		{
			if (formSum.termCount % BLOCK_SIZE == 0) //need to allocate a new block for coeffs and exponents
			{
				linCoeff->next = (cBlock<T>*) malloc (sizeof(cBlock<T>));
				linForm->next = (lBlock*) malloc (sizeof(lBlock));
				linForm = linForm->next; linCoeff = linCoeff->next;
				linForm->next = NULL; linCoeff->next = NULL;
				linForm->data = new vec_ZZ[BLOCK_SIZE];
				linCoeff->data = new T[BLOCK_SIZE];
			}
			linForm->data[formSum.termCount % BLOCK_SIZE].SetLength(formSum.varCount);
			VectorCopy(linForm->data[formSum.termCount % BLOCK_SIZE], coeffs, formSum.varCount);
			linForm->degree[formSum.termCount % BLOCK_SIZE] = degree;
			linCoeff->data[formSum.termCount % BLOCK_SIZE] = coef;
			formSum.termCount++;
		}
		else //found this linear form
		{
			linCoeff->data[myIndex % BLOCK_SIZE] += coef;
		}
	}
	else
	{
		formSum.cHead = (cBlock<T>*) malloc (sizeof(cBlock<T>));
		formSum.lHead = (lBlock*) malloc (sizeof(lBlock));
		linForm = formSum.lHead; linCoeff = formSum.cHead;
		linForm->next = NULL; linCoeff->next = NULL;
		linForm->data = new vec_ZZ[BLOCK_SIZE];
		linCoeff->data = new T[BLOCK_SIZE];
		linForm->data[formSum.termCount % BLOCK_SIZE].SetLength(formSum.varCount);
		VectorCopy(linForm->data[formSum.termCount % BLOCK_SIZE], coeffs, formSum.varCount);
		linForm->degree[formSum.termCount % BLOCK_SIZE] = degree;
		linCoeff->data[formSum.termCount % BLOCK_SIZE] = coef;
		formSum.termCount++;
	}
}

#endif
