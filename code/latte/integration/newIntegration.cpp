
#include "newIntegration.h"
#include <NTL/vec_ZZ.h>
#include <NTL/ZZ.h>
#include <NTL/mat_ZZ.h>
#include <iostream>

using namespace std;

NTL_CLIENT

//this function deletes space from a given string
void delSpace(string &line)
{
	for (int i=0;i<line.length();i++)
	{
		while ((i<line.length())&&(line.at(i)==32)) {line.erase(i,1);};
	};
};

//this function converts a given string into a simlexZZ mySimplex
//for example, string [[0,0],[1,1],[7,8]] is converted to a two-dimensional vector of ZZs.
void convertToSimplex(simplexZZ &mySimplex, string line)
{
	delSpace(line);
	int index,i,t,j,c;
	string temp,subtemp;
	t=2;mySimplex.d=1;
	t=line.find("[",t)+1;t=line.find("[",t)+1;
	temp=line.substr(t,line.find("]",t)-t);
	for (i=0;i<temp.length();i++) mySimplex.d+=(temp.at(i)==',');
	c=0;
	for (i=0;i<line.length();i++) c+=(line.at(i)==']');
	if (c-2!=mySimplex.d) {cout<<"The d-simplex should have d+1 vertices. Please check."<<endl; exit(1);};
 	(mySimplex.s).SetLength(mySimplex.d+1);
	index=1;
	for (i=0;i<=mySimplex.d;i++)
	{
		temp=line.substr(index,line.find("]",index)-index+1);
		c=0;
		for (j=0;j<temp.length();j++) c+=(temp.at(j)==',');
		if (c!=mySimplex.d-1) {cout<<"Each vertex should have d coordinates. Please check."<<endl; exit(1);};
		(mySimplex.s[i]).SetLength(mySimplex.d);
		t=1;
		for (j=0;j<mySimplex.d-1;j++)
		{
			subtemp=temp.substr(t,temp.find(",",t)-t);
			t=temp.find(",",t)+1;
			mySimplex.s[i][j]=to_ZZ(subtemp.c_str());
		};
		subtemp=temp.substr(t,temp.find(",",t)-t+1);
		t=temp.find(",",t);
		mySimplex.s[i][mySimplex.d-1]=to_ZZ(subtemp.c_str());
		index=line.find("]",index)+2;
	};
	mat_ZZ matt;	
	matt.SetDims(mySimplex.d,mySimplex.d);
	for (i=1;i<=mySimplex.d;i++) matt[i-1]=mySimplex.s[i]-mySimplex.s[0];
	mySimplex.v=determinant(matt);if (mySimplex.v<0) mySimplex.v=-mySimplex.v;
};

//The purpose of this function is to modify a given fraction a/b. 
//Input: l is the exponents of a linear form, mySimplex is the Simplex we are integrating over with d+1 vertices
//m is the power the linear form is raised to
//coe is the coefficient of a linear form
//de is the extra factor in the formulae that we we multiply the result by

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
	bool repeat[1000];
	for (i=0;i<=mySimplex.d;i++)
	{
		sum=0; for (j=0;j<mySimplex.d;j++) sum=sum+l[j]*mySimplex.s[i][j];
		inner_Pro[i]=sum;
		repeat[i]=0;
		for (j=0;j<i;j++) if (inner_Pro[j]==inner_Pro[i]) {repeat[i]=1;break;};//record repetitions
	};//stores inner product for use
	for (i=0;i<=mySimplex.d;i++)
	if (!repeat[i])
	{
		sum_Nu[i]=1;for (j=0;j<m+mySimplex.d;j++) sum_Nu[i]=sum_Nu[i]*inner_Pro[i];
		sum_De[i]=1;for (j=0;j<=mySimplex.d;j++) if (i!=j) sum_De[i]=sum_De[i]*(inner_Pro[i]-inner_Pro[j]);
		if ((sum_Nu[i]<0)&&(sum_De[i]<0)) {sum_Nu[i]=-sum_Nu[i];sum_De[i]=-sum_De[i];};
		if (sum_De[i]==0)
		{
			vec_ZZ ProDiff;
			ProDiff.SetLength(mySimplex.d+1);
			for (j=0;j<=mySimplex.d;j++) ProDiff[j]=inner_Pro[i]-inner_Pro[j];
			computeResidue(mySimplex.d,m,ProDiff,inner_Pro[i],sum_Nu[i],sum_De[i]);
		} 
		if (sum_De[i]!=0) {lcm=lcm*sum_De[i]/(GCD(lcm,sum_De[i]));};
	};
	for (i=0;i<=mySimplex.d;i++)
	if ((!repeat[i])&&(sum_De[i]!=0))
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
};

//This function computes a given fraction a/b, the integral of the linear form forms, over the simplex mySimplex
void integrateLinFormSum(ZZ& numerator, ZZ& denominator, const linFormSum &forms , const simplexZZ &mySimplex)
{
	ZZ v,de,counter,tem,coe;
	int i,j,index,k,m;
	vec_ZZ l;
	if (forms.varCount!=mySimplex.d) {cout<<"The dimensions of the polynomial and simplex don't match. Please check!"<<forms.varCount<<"<>"<<mySimplex.d<<endl;exit(1);};
	l.SetLength(mySimplex.d);
	numerator=0;
	denominator=0;
	BTrieIterator<ZZ, ZZ>* it = new BTrieIterator<ZZ, ZZ>();
	it->setTrie(forms.myForms, forms.varCount);
	it->begin();
	
	term<ZZ, ZZ>* temp;
	while (temp = it->nextTerm())
	{
		coe=temp->coef;m=temp->degree;		//obtain coefficient, power
		l.SetLength(temp->length);		//obtain exponent vector
		for (j=0;j<temp->length;j++)
		{
			l[j]=temp->exps[j];
		}
		de=1;
		for (i=1;i<=mySimplex.d+m;i++)
		{
			de=de*i;			
		};					//de is (d+m)!. Note this is different from the factor in the paper because in our 								storage of a linear form, any coefficient is automatically adjusted by m!

		update(numerator,denominator,l,mySimplex,m,coe,de);//We are ready to compute the integral of one linear form over the simplex
	};
	delete temp;
	if (denominator<0) {denominator *= to_ZZ(-1); numerator *= to_ZZ(-1);};
};

//this is the old block data structure integration
void _integrateLinFormSum(ZZ& numerator, ZZ& denominator, const _linFormSum &forms , const simplexZZ &mySimplex)
{
	ZZ v,de,counter,tem,coe;
	int i,j,index,k,m;
	vec_ZZ l;
	if (forms.varCount!=mySimplex.d) {cout<<"The dimensions of the polynomial and simplex don't match. Please check!"<<forms.varCount<<"<>"<<mySimplex.d<<endl;exit(1);};
	l.SetLength(mySimplex.d);

	lBlock* formTmp = forms.lHead; cBlock<ZZ>* coeffTmp = forms.cHead;
	for (int i = 0; i < forms.termCount; i++)
	{
		if (i > 0 && i % BLOCK_SIZE == 0)			//just a different way of traversing the data structure
		{
			formTmp = formTmp->next;
			coeffTmp = coeffTmp->next;
		}
		coe=coeffTmp->data[i % BLOCK_SIZE];m=formTmp->degree[i % BLOCK_SIZE];
		for (int j = 0; j < forms.varCount; j++)
		{
			l[j]=formTmp->data[i % BLOCK_SIZE][j];
		}
		de=1;
		for (int j=1;j<=mySimplex.d+m;j++)
		{
			de=de*j;
		};
		update(numerator,denominator,l,mySimplex,m,coe,de);
	};
};
void integrateMonomialSum(ZZ &a, ZZ &b, monomialSum &monomials, const simplexZZ &mySimplex)//integrate a polynomial stored as a Burst Trie
{
	linFormSum forms;
	forms.termCount = 0;
	forms.varCount = monomials.varCount;
	BTrieIterator<ZZ, int>* it = new BTrieIterator<ZZ, int>();
	it->setTrie(monomials.myMonomials, monomials.varCount);
	decompose(it, forms);				//decomposition
	delete it;
	integrateLinFormSum(a,b,forms, mySimplex);
};
/*void _integrateMonomialSum(ZZ &a, ZZ &b, _monomialSum &monomials, const simplexZZ &mySimplex)
{
	_linFormSum forms;
	forms.termCount = 0;
	forms.varCount = monomials.varCount;		
	//for (int i=0;i<monomials.termCount;i++) _decompose(monomials, forms, i);
	_integrateLinFormSum(a,b,forms, mySimplex);
};*/
