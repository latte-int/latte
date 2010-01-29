#include <NTL/vec_ZZ.h>
#include <NTL/ZZ.h>
#include "PolyRep.h"

void addToFraction(const int &s[50],int n, const ZZ &t, const int &index[50], const int &counter[50], ZZ &a, ZZ &b)
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
};
void computeResidue(int d, ZZ M, const vec_ZZ &innerProDiff, const ZZ &p, ZZ &a, ZZ &b)
{
     if (p==0) {a=0;b=1; return;}; //vertex vanishes, return 0;
     int k,i,j;
     int counter[50];//counter counts number of appearances of each index[i]
     int index[50];//collecting different terms in the innerProDiff passed in
     int s[50]; //s is the array we use to enumerate all possible weak decompostitions of counter[0]
     bool found;
     k=1;
     index[k-1]=0;
     counter[k-1]=0;
     for (i=0;i<=d;i++)
     {
         found=0;
         for (j=0;j<k;j++)
         {
             if (innerProDiff[i]==index[j]) {counter[j]++; found=1; break;}
         };
         if (!found) {k++;index[k-1]=innerProDiff[i];counter[k-1]++;};
     };
     counter[0]--;   //excluding the vertex itself
     for (i=0;i<=k;i++) {s[i]=0;}
     enumerate(s,0,k,counter[0],p,M+d,index,counter,a,b);
};
