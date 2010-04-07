//will define multivariate polynomial representation, as well as input and output functions
#ifndef POLYREP_H
#define POLYREP_H

#define BLOCK_SIZE 64
#include <NTL/vec_ZZ.h>
#include <NTL/ZZ.h>

NTL_CLIENT

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

void _loadMonomials(_monomialSum&, const string&);
//string parsing
void _parseMonomials(_MonomialConsumer<ZZ>*, const string&);
//data structure operations
template <class T> void _insertMonomial(const T&, int*, _monomialSum&);
string _printMonomials(const _monomialSum&);
void _destroyMonomials(_monomialSum&);

void _loadLinForms(_linFormSum&, const string);
//string parsing
void _parseLinForms(_FormSumConsumer<ZZ>*, const string&);
//data structure operations
template <class T> void _insertLinForm(const T& coef, int degree, const vec_ZZ& coeffs, _linFormSum&);
string _printLinForms(const _linFormSum&);
void _destroyLinForms(_linFormSum&);

void _decompose(_monomialSum&, _linFormSum&, int);

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
