#ifndef CONSUMERS_H
#define CONSUMERS_H

#include "blockReps.h"
#include "burstTrie.h"

//trie-based

//define base consumers to enable flexible string parsing
template <class T> class MonomialConsumer {
public:
  // Take monomial and consume it.
  virtual void ConsumeMonomial(const T&, int*) = 0;
  virtual void setDimension(int dimension) = 0;
  virtual int getDimension() = 0;
};

template <class T> class FormSumConsumer {
public:
  // Take linear form and consume it.
  //virtual int ConsumeLinForm(const T&, int, int*) = 0;
  virtual void ConsumeLinForm(const T&, int, const vec_ZZ&) = 0;
  virtual void setDimension(int dimension) = 0;
  virtual int getDimension() = 0;
};

//consumers for loading data structures
template <class T> class MonomialLoadConsumer : public MonomialConsumer<T> {
public:
  MonomialLoadConsumer() {}
  // Take monomial and consume it.
  void ConsumeMonomial(const T& coef, int* exps) { insertMonomial(coef, exps, *monomials); }
  void setMonomialSum(monomialSum& mySum) { monomials = &mySum; }
  void setDimension(int dimension) { if (monomials) { monomials->varCount = dimension; } }
  int getDimension() { if (monomials) { return monomials->varCount; } else { return 0; } }
private:
  monomialSum* monomials;
};

template <class T> class FormLoadConsumer : public FormSumConsumer<T> {
public:
  FormLoadConsumer() {}
  // Take linear form and consume it.
  void ConsumeLinForm(const T& coefficient, int degree, const vec_ZZ& coefs) { insertLinForm(coefficient, degree, coefs, *formSum); }
  void setFormSum(linFormSum& myForms) { formSum = &myForms; }
  void setDimension(int dimension) { if (formSum) { formSum->varCount = dimension; } }
  int getDimension() { if (formSum) { return formSum->varCount; } else { return 0; } }
  ~FormLoadConsumer() {}
private:
  linFormSum* formSum;
};

//block-based

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

//define base consumers to enable flexible string parsing
template <class T> class _MonomialConsumer {
public:
  // Take monomial and consume it.
  virtual void ConsumeMonomial(const T&, int*) = 0;
  virtual void setDimension(int dimension) = 0;
  virtual int getDimension() = 0;
};

template <class T> class _FormSumConsumer {
public:
  // Take linear form and consume it.
  //virtual int ConsumeLinForm(const T&, int, int*) = 0;
  virtual void ConsumeLinForm(const T&, int, const vec_ZZ&) = 0;
  virtual void setDimension(int dimension) = 0;
  virtual int getDimension() = 0;
};

//consumers for loading data structures
template <class T> class _MonomialLoadConsumer : public _MonomialConsumer<T> {
public:
  _MonomialLoadConsumer() {}
  // Take monomial and consume it.
  void ConsumeMonomial(const T& coef, int* exps) { _insertMonomial<ZZ>(coef, exps, *monomials); }
  void setMonomialSum(_monomialSum& mySum) { monomials = &mySum; }
  void setDimension(int dimension) { if (monomials) { monomials->varCount = dimension; } }
  int getDimension() { if (monomials) { return monomials->varCount; } else { return 0; } }
private:
  _monomialSum* monomials;
};

template <class T> class _FormLoadConsumer : public _FormSumConsumer<T> {
public:
  _FormLoadConsumer() {}
  // Take linear form and consume it.
  void ConsumeLinForm(const T& coefficient, int degree, const vec_ZZ& coefs) { _insertLinForm<ZZ>(coefficient, degree, coefs, *formSum); }
  void setFormSum(_linFormSum& myForms) { formSum = &myForms; }
  void setDimension(int dimension) { if (formSum) { formSum->varCount = dimension; } }
  int getDimension() { if (formSum) { return formSum->varCount; } else { return 0; } }
  ~_FormLoadConsumer() {}
private:
  _linFormSum* formSum;
};

#endif
