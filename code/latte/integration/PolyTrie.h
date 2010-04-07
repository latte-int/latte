//will define multivariate polynomial representation, as well as input and output functions
#ifndef POLYTRIE_H
#define POLYTRIE_H

#include "burstTrie.h"
#include <NTL/ZZ.h>
#include <NTL/vec_ZZ.h>

NTL_CLIENT

struct monomialSum
{
	int termCount, varCount;
	BurstTrie<ZZ, int>* myMonomials;
};

//linear forms: sort on degree first, then the form coefficients
struct linFormSum
{
	int termCount, varCount;
	BurstTrie<ZZ, ZZ>* myForms;
};

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

void loadMonomials(monomialSum&, const string&);
//string parsing
void parseMonomials(MonomialConsumer<ZZ>*, const string&);
//data structure operations
void insertMonomial(const ZZ&, int*, monomialSum&);
string printMonomials(const monomialSum&);
void destroyMonomials(monomialSum&);

void loadLinForms(linFormSum&, const string);
//string parsing
void parseLinForms(FormSumConsumer<ZZ>*, const string&);
//data structure operations
void insertLinForm(const ZZ& coef, int degree, const vec_ZZ& coeffs, linFormSum&);
string printLinForms(const linFormSum&);
void destroyLinForms(linFormSum&);

void decompose(monomialSum&, linFormSum&, BurstTerm<ZZ, int>*);
void decompose(monomialSum &myPoly, linFormSum &lForm);

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

/*
class Decomposer: public TrieIterator<ZZ, int> {
    public:
    monomialSum* input;
    linFormSum* output;
     
    void init(monomialSum *myPoly, linFormSum *lForm)
    {
        input = myPoly;
	lForm->varCount = myPoly->varCount;
	lForm->termCount = 0;
	output = lForm;
    }
    
    void consumeTerm(term<ZZ, int>* myTerm)
    {
	//use curprefix, etc.
        decompose(*input, *output, getFullTerm(myTerm));
    }
    
    burstTrie<ZZ, ZZ>* getResults()
    {
	return output->myForms;
    }
    
    void decompose(monomialSum &myPoly, linFormSum &lForm, term<ZZ, int>* myTerm)
    {
	vec_ZZ myExps; myExps.SetLength(lForm.varCount);
	
	if (myTerm->expCount == 0) //constant
	{
		for (int j = 0; j < lForm.varCount; j++) { myExps[j] = 0; }
		insertLinForm(*myTerm->coef, 0, myExps, lForm);
		return;
	}
	
	ZZ formsCount = to_ZZ(myTerm->exponents[0] + 1); //first exponent
	int totalDegree = myTerm->exponents[0];
	for (int i = 1; i < dimension; i++)
	{
		formsCount *=  myTerm->exponents[i] + 1;
		totalDegree +=  myTerm->exponents[i];
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
		for (myIndex = 0; counter[myIndex] > myTerm->exponents[myIndex]; myIndex++)
		{
			counter[myIndex] = 0;
			binomCoeffs[myIndex] = to_ZZ(1);
			counter[myIndex+1] += 1;
		}
		binomCoeffs[myIndex] *= myTerm->exponents[myIndex] - counter[myIndex] + 1; // [n choose k] = [n choose (k - 1) ] * (n - k + 1)/k
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
		temp *= *myTerm->coef;
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
};*/

#endif
