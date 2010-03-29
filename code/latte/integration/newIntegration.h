#ifndef NEWINTEGRATION_H
#define NEWINTEGRATION_H
#include <NTL/vec_vec_ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/ZZ.h>
#include "PolyTrie.h"
#include "residue.h"

NTL_CLIENT

struct simplexZZ
{
	int d;
	vec_vec_ZZ s;
	ZZ v;
};

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
  void update(vec_ZZ l,int m, ZZ coe, ZZ de);
};

template <class T>
void FormIntegrateConsumer<T>::ConsumeLinForm(const ZZ& coefficient, int degree, const vec_ZZ& coefs)
{
	ZZ de = to_ZZ(1);
	for (int i=1;i<=mySimplex->d+degree;i++)
	{
		de=de*i;
	};
	update(coefs,degree,coefficient,de);
}

template <class T>
void FormIntegrateConsumer<T>::update(vec_ZZ l,int m, ZZ coe, ZZ de)
{
	ZZ sum,lcm,total,g,tem;
	int i,j;
	vec_ZZ inner_Pro,sum_Nu,sum_De;
	inner_Pro.SetLength(mySimplex->d+1);
	sum_Nu.SetLength(mySimplex->d+1);
	sum_De.SetLength(mySimplex->d+1);
	total=0;
	lcm=1;
	for (i=0;i<=mySimplex->d;i++)
	{
		sum=0; for (j=0;j<mySimplex->d;j++) {sum=sum+l[j]*mySimplex->s[i][j];};
		inner_Pro[i]=sum;
	};//stores inner product for use
	for (i=0;i<=mySimplex->d;i++)
	{
		sum_Nu[i]=1;for (j=0;j<m+mySimplex->d;j++) sum_Nu[i]=sum_Nu[i]*inner_Pro[i];
		sum_De[i]=1;for (j=0;j<=mySimplex->d;j++) if (i!=j) sum_De[i]=sum_De[i]*(inner_Pro[i]-inner_Pro[j]);
		if (sum_De[i]==0) {denominator = to_ZZ(0); return;}; //irregular
		lcm=lcm*sum_De[i]/(GCD(lcm,sum_De[i]));
	};
	for (i=0;i<=mySimplex->d;i++)
	{
		total+=sum_Nu[i]*(lcm/sum_De[i]);
	};
	lcm=lcm*de;
	total=total*mySimplex->v*coe;
	if (numerator==0) {numerator=total;denominator=lcm;}
	else if ((lcm!=0)&&(denominator!=0)) {tem=denominator*lcm/GCD(denominator,lcm);numerator=numerator*tem/denominator+total*tem/lcm;denominator=tem;};	
	g=GCD(numerator,denominator);
	if (g!=0) 
	{
	numerator=numerator/g;
	denominator=denominator/g;};
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

void delSpace(string &line);
void convertToSimplex(simplexZZ &mySimplex, string line);
//the below should be changed to more comprehensible names - for example, integrateForms + integrateStringForms?
void integrateListString(ZZ &a, ZZ &b, string line, const simplexZZ &mySimplex);
void integrateList(ZZ &a, ZZ &b, const linFormSum &forms, const simplexZZ &mySimplex);
//this does decomposition, then calls the above functions. why is this needed?
void integrateFlatVectorString(ZZ &a, ZZ &b, string line, const simplexZZ &mySimplex);
void integrateFlatVector(ZZ &numerator, ZZ &denominator, monomialSum &monomials, const simplexZZ &mySimplex);

class Integrator: public TrieIterator<ZZ, ZZ> {
    public:
    ZZ numerator, denominator;
    const linFormSum *myForm;
    const simplexZZ *mySimplex;
     
    void init(const linFormSum *lForm, const simplexZZ *simplex)
    {
        myForm = lForm;
	mySimplex = simplex;
	numerator = denominator = 0;
    }
    
    void consumeTerm(term<ZZ, ZZ>* myTerm)
    {
	vec_ZZ myCoeffs;
	myCoeffs.SetLength(TrieIterator<ZZ, ZZ>::dimension);
	for (int i = 0; i < myTerm->expCount; i++)
	{
		myCoeffs[i] = myTerm->exponents[i];
	}

	ZZ fact = to_ZZ(1);
	for (int i = 1; i <= mySimplex->d + myTerm->degree;i++)
	{
		fact*=i;
	}
        update(numerator, denominator, myCoeffs, *mySimplex, myTerm->degree, *myTerm->coef, fact);
    }
    
    void getResults(ZZ &num, ZZ &den)
    {
	num = numerator;
	den = denominator;
    }
    
    void update(ZZ &a, ZZ &b, vec_ZZ l, simplexZZ mySimplex,int m, ZZ coe, ZZ de)
    {
	ZZ sum,lcm,total,g,tem;
	int i,j;
	vec_ZZ inner_Pro,sum_Nu,sum_De;
	inner_Pro.SetLength(mySimplex.d+1);
	sum_Nu.SetLength(mySimplex.d+1);
	sum_De.SetLength(mySimplex.d+1);
	total=0;
	lcm=1;
	bool repeat[50];
	for (i=0;i<=mySimplex.d;i++)
	{
		sum=0; for (j=0;j<mySimplex.d;j++) sum=sum+l[j]*mySimplex.s[i][j];
		inner_Pro[i]=sum;
		repeat[i]=0;
		for (j=0;j<i;j++) if (inner_Pro[j]==inner_Pro[i]) {repeat[i]=1;break;};
	};//stores inner product for use
	for (i=0;i<=mySimplex.d;i++)
	if (!repeat[i])
	{
		sum_Nu[i]=1;for (j=0;j<m+mySimplex.d;j++) sum_Nu[i]=sum_Nu[i]*inner_Pro[i];
		sum_De[i]=1;for (j=0;j<=mySimplex.d;j++) if (i!=j) sum_De[i]=sum_De[i]*(inner_Pro[i]-inner_Pro[j]);
		if ((sum_Nu[i]<0)&&(sum_De[i]<0)) {sum_Nu[i]=-sum_Nu[i];sum_De[i]=-sum_De[i];};
		if (sum_De[i]==0)
		{
			//cout<<l<<" is not regular at "<<inner_Pro[i]<<". Residue:";
			vec_ZZ ProDiff;
			ProDiff.SetLength(mySimplex.d+1);
			for (j=0;j<=mySimplex.d;j++) ProDiff[j]=inner_Pro[i]-inner_Pro[j];
			//cout<<"The inner product difference is"<<ProDiff<<" inner product at vertex is "<<inner_Pro[i]<<". Now compute residue:"<<endl;
			computeResidue(mySimplex.d,m,ProDiff,inner_Pro[i],sum_Nu[i],sum_De[i]);
		} //irregular
		//else cout<<l<<" is regular at "<<inner_Pro[i]<<". Sum:";
		//cout<<sum_Nu[i]<<"/"<<sum_De[i]<<endl;
		if (sum_De[i]!=0) {lcm=lcm*sum_De[i]/(GCD(lcm,sum_De[i]));};
	};
	for (i=0;i<=mySimplex.d;i++)
	if (!repeat[i])
	{
		total+=sum_Nu[i]*(lcm/sum_De[i]);
	};
	lcm=lcm*de;
	total=total*mySimplex.v*coe;
	if (a==0) {a=total;b=lcm;}
	else if ((lcm!=0)&&(b!=0)) {tem=b*lcm/GCD(b,lcm);a=a*tem/b+total*tem/lcm;b=tem;};	
	g=GCD(a,b);
	if (g!=0) 
	{
	a=a/g;
	b=b/g;};
      }
};

#endif
