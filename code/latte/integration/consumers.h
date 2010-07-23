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

template <class T>
void _insertLinForm(const T& coef, int degree, const vec_ZZ& coeffs, _linFormSum& formSum);

template <class T>
void _insertMonomial(const T& coefficient, int* exponents, _monomialSum& monomials);

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
  void ConsumeMonomial(const T& coef, int* exps) { _insertMonomial<RationalNTL>(coef, exps, *monomials); }
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
  void ConsumeLinForm(const T& coefficient, int degree, const vec_ZZ& coefs) { _insertLinForm<T>(coefficient, degree, coefs, *formSum); }
  void setFormSum(_linFormSum& myForms) { formSum = &myForms; }
  void setDimension(int dimension) { if (formSum) { formSum->varCount = dimension; } }
  int getDimension() { if (formSum) { return formSum->varCount; } else { return 0; } }
  ~_FormLoadConsumer() {}
private:
  _linFormSum* formSum;
};

#endif
