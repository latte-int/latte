#include <NTL/vec_ZZ.h>
#include <NTL/ZZ.h>

#include "PolyRep.h"
#include "multiply.h"

/*void addToFraction(const int &s[50],int n, const ZZ &t, const int &index[50], const int &counter[50], ZZ &a, ZZ &b)
{
     ZZ de,nu;
     de=1;nu=1;
     int i,j,w;
     j=1;
     for (i=0;i<=k;i++)
     {
         for (w=1;w<=s[i];w++)
         {
             nu*=j;j++;
             de*=w;
         };
     };
     nu/=de;de=1;
     if (((n-s[0])%2)==1) nu=-nu;
     for (i=0;i<s[0];i++) nu*=(t-i);
     for (i=1;i<k;i++)
     {
          for (j=1;j<s[i];j++) nu*=(counter[i]+j);
     };
     nu*=Power(p,t-s[0]);
     for (i=1;i<=k;i++)
     {
	de*=Power(index[i],counter[i]+s[i]);
     };
     ZZ g=GCD(nu,de);
     nu=nu/g;
     de=de/g;
     if (b==0) {a=nu;b=de;return;};
     ZZ lcm=b*de/GCD(b,de);
     a=a*lcm/b+nu*lcm/de;
     b=lcm;
     g=GCD(a,b);
     a=a/g;
     b=b/g;
     if ((a<0)||(b<0)) {a=-a;b=-b;};
}

     
void enumerate(const int &s[50], int i, int k, int n, const ZZ &p, const ZZ &t, const int &index[50], const int &counter[50], ZZ &a, ZZ &b)
{
     if (i==k) {s[i]=n; addToFraction(s,n,t,index,counter,a,b); return;};
     for (int j=0;j<=n;j++)
     {
         s[i]=j;enumerate(s,i+1,k,n-j,p,t,index,counter,a,b);
     };
};*/
ZZ AChooseB(int a,int b)
{
	ZZ t=to_ZZ(1);
	if (b>a) return to_ZZ(0);
	if (2*b>a) b=a-b;
	for (int i=1;i<=b;i++)
	{t=t*(a-i+1)/i;};
	return t;
}		
void computeResidue(int d, int M, const vec_ZZ &innerProDiff, const ZZ &p, ZZ &a, ZZ &b)
{
     if (p==0) {a=0;b=1; return;}; //vertex vanishes, return 0;
     int k,i,j;
     int counter[50];//counter counts number of appearances of each index[i]
     vec_ZZ index;//collecting different terms in the innerProDiff passed in
     int s[50]; //s is the array we use to enumerate all possible weak decompostitions of counter[0]
     bool found;
     ZZ de,nu,c,g;
     int e[1];
     int mindeg[1];
     int maxdeg[1];	
     k=1;
     index.SetLength(d);	
     index[k-1]=0;
     counter[k-1]=0;
     for (i=0;i<=d;i++)
     {
         found=0;
         for (j=0;j<k;j++)
         {
             if (innerProDiff[i]==index[j]) {counter[j]++; found=1; break;}
         };
         if (!found) {k++;index[k-1]=innerProDiff[i];counter[k-1]=1;};
     };
     counter[0]--;   //excluding the vertex itself
			//so far I've been doing book keeping stuff, i.e. how many different differences of each kind are there? index stores the differences and counter stores the multiplicity. 0th entry means 0 itself. It is very important because it's the multiplicity of a vertex itself.
			//actual calculations, I want the appropriate coefficient in a product of one polynomial and some power series. (everything is truncated).
     nu=1;de=1;
     for (i=1;i<=counter[0];i++) nu*=i;
     for (j=1;j<=k-1;j++) de=de*Power_ZZ(index[j],counter[j]+counter[0]);
     monomialSum m1;
     monomialSum sub;					//sub is the substitution for m1, which alternatively stores the product for each other
     m1.varCount=1;m1.termCount=0;
     sub.varCount=1;sub.termCount=0;
     for (i=0;i<=counter[0];i++) 
     {
	c=AChooseB(M+d,i)*Power_ZZ(p,M+d-i);
	e[0]=i;
	insertMonomial<ZZ>(c,e,m1);
     };
     //cout<<printMonomials(m1)<<endl;
     for (i=1;i<k;i++)
     {
	monomialSum m2;
	m2.varCount=1;m2.termCount=0;
	for (j=0;j<=counter[0];j++)
	{
		c=AChooseB(counter[i]+j-1,j)*Power_ZZ(index[i],counter[0]-j);
		if (j % 2==1) c=-c;
		e[0]=j;
		insertMonomial<ZZ>(c,e,m2);
	};
	mindeg[0]=0;
	maxdeg[0]=counter[0];
	if (i % 2==1) {multiply<ZZ>(m1,m2,sub,mindeg,maxdeg);}//cout<<"times "<<printMonomials(m2)<<" gives "<<printMonomials(sub)<<endl;}
	else {multiply<ZZ>(sub,m2,m1,mindeg,maxdeg);}//cout<<"times "<<printMonomials(m2)<<"gives "<<printMonomials(m1)<<endl;};
     };
     ZZ findCoeff;	
     eBlock* myExps; cBlock<ZZ>* myCoeffs;
     if (k % 2==1) 					//choose which one to pick result from
     {myExps = m1.eHead; myCoeffs = m1.cHead;
	     for (i=0;i<m1.termCount;i++)
     	{
		if (i>0 && i % BLOCK_SIZE ==0) 
		{
			myExps = myExps->next; myCoeffs=myCoeffs->next;
		};
		if (myExps->data[i % BLOCK_SIZE]== counter[0]) {findCoeff=myCoeffs->data[i % BLOCK_SIZE];break;};
      	};
      }
      else 	
      {myExps = sub.eHead; myCoeffs = sub.cHead;
	     for (i=0;i<sub.termCount;i++)
     	{
		if (i>0 && i % BLOCK_SIZE ==0) 
		{
			myExps = myExps->next; myCoeffs=myCoeffs->next;
		};
		if (myExps->data[i % BLOCK_SIZE]== counter[0]) {findCoeff=myCoeffs->data[i % BLOCK_SIZE];break;};
      	};
      };
	
	
      a=nu*findCoeff;
      b=de;
      g=GCD(a,b);
      if (g!=0) {a=a/g;b=b/g;};      
      return;
};
