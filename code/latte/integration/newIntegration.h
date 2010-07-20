#ifndef NEWINTEGRATION_H
#define NEWINTEGRATION_H
#include <NTL/vec_vec_ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/ZZ.h>
#include "PolyTrie.h"
#include "PolyRep.h"
//#include "PolyRep.cpp"
//#include "PolyTrie.cpp"
#include "residue.h"

NTL_CLIENT

struct simplexZZ
{
					// Brandon's notes.
	int d;			// dimension of the space? so s.length = d + 1?
	vec_vec_ZZ s;	// s[i] = vector of the ith vertex
	ZZ v; 			//valuation/integral answer?
};

void update(ZZ &a, ZZ &b, vec_ZZ l, simplexZZ mySimplex,int m, ZZ coe, ZZ de);
void delSpace(string &line);
void convertToSimplex(simplexZZ&, string);
void integrateLinFormSum(ZZ &a, ZZ &b, PolyIterator<ZZ, ZZ>* it, const simplexZZ &mySimplex);
void integrateMonomialSum(ZZ &numerator, ZZ &denominator, monomialSum &monomials, const simplexZZ &mySimplex);
void _integrateMonomialSum(ZZ &numerator, ZZ &denominator, _monomialSum &monomials, const simplexZZ &mySimplex);

template <class T>
class FormIntegrateConsumer : public FormSumConsumer<T> {
public:
  FormIntegrateConsumer() { }
  // Take linear form and consume it.
  void ConsumeLinForm(const ZZ& coefficient, int degree, const vec_ZZ& coefs);
  void setFormSum(const string& myForms) { linForms = myForms; }
  string getFormSum() { return linForms; }
  void setDimension(int dimension) { mySimplex->d = dimension; }//also stored in mySimplex
  int getDimension() { return mySimplex->d; }
  void setSimplex(simplexZZ& simplex) { mySimplex = &simplex; numerator = to_ZZ(0); denominator = to_ZZ(0); }
  void getResults(ZZ& num, ZZ& den);
  ~FormIntegrateConsumer() {}
private:
  string linForms;
  simplexZZ* mySimplex;
  ZZ numerator, denominator;
};

template <class T>
void FormIntegrateConsumer<T>::ConsumeLinForm(const ZZ& coefficient, int degree, const vec_ZZ& coefs)
{
	ZZ de = to_ZZ(1);
	for (int i=1;i<=mySimplex->d+degree;i++)
	{
		de=de*i;
	};
	update(numerator, denominator, coefs, *mySimplex, degree, coefficient, de);
}

template <class T>
void FormIntegrateConsumer<T>::getResults(ZZ& num, ZZ& den)
{
	if (denominator < 0)
	{
		num = to_ZZ(-1) * numerator; den = to_ZZ(-1) * denominator;
	}
	else
	{
		num = numerator; den = denominator;
	}
}

#endif
