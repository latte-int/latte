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

void decompose(term<ZZ, int>*, linFormSum&);
void decompose(BTrieIterator<ZZ, int>* it, linFormSum&);

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

#endif
