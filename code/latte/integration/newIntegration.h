#include "PolyRep.h"

#include <NTL/vec_ZZ.h>
#include <NTL/ZZ.h>
#include <NTL/vec_vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include <iostream>

using namespace std;

NTL_CLIENT

void print_Integrate(int d, const linearPoly &lForm , string &line)
{
  ZZ v,de,sum,lcm,total,g,a,b,counter,m,tem;
	int i,j,index,k,t;
	vec_ZZ inner_Pro,sum_Nu,sum_De,l;
	vec_vec_ZZ s;
	string temp,subtemp;
	//convert the string into a vec_vec_ZZ
 	s.SetLength(d+1);
	index=1;
	cout<<line<<endl;
	for (i=0;i<=d;i++)
	{
		temp=line.substr(index,line.find("]",index)-index+1);
		s[i].SetLength(d);
		t=1;
		for (j=0;j<d-1;j++)
		{
			subtemp=temp.substr(t,temp.find(",",t)-t);
			t=temp.find(",",t)+1;
			s[i][j]=to_ZZ(subtemp.c_str());
		};
		subtemp=temp.substr(t,temp.find(",",t)-t+1);
		t=temp.find(",",t);
		s[i][d-1]=to_ZZ(subtemp.c_str());
		index=line.find("]",index)+2;
	};
	inner_Pro.SetLength(d);
	sum_Nu.SetLength(d);
	sum_De.SetLength(d);
	l.SetLength(d);
	lcm=1; //least common multiple
	total=0; //total numerator
	v=1;
	mat_ZZ matt;	
	matt.SetDims(d,d);
	for (i=1;i<=d;i++) matt[i-1]=s[i]-s[0];
	v=determinant(matt);if (v<0) v=-v;
	lBlock* temForm=lForm.lHead;
	cBlock* temCoef=lForm.cHead;
	k=-1;
	a=0;
	b=0;
	for (counter=0;counter<lForm.termCount;counter++)
	{
		k++;
		if ((k>0)&&(k % BLOCK_SIZE==0)) {temForm=temForm->next;temCoef=temCoef->next;k=0;}
		
		for (i=0;i<d;i++) l[i]=temForm->data[k][i];
		m=temForm->degree[k];
		de=1;
		total=0;
		lcm=1;
		for (i=1;i<=d;i++)
		{
			de=de*(m+j);
		};
		for (i=0;i<=d;i++)
		{
			sum=0; for (j=0;j<d;j++) {sum=sum+l[j]*s[i][j];};
			inner_Pro[i]=sum;
		};//stores inner product for use
		for (i=0;i<=d;i++)
		{
			sum_Nu[i]=1;for (j=0;j<m+d;j++) sum_Nu[i]=sum_Nu[i]*inner_Pro[i];
			sum_De[i]=1;for (j=0;j<=d;j++) if (i!=j) sum_De[i]=sum_De[i]*(inner_Pro[i]-inner_Pro[j]);
			if (sum_De[i]==0) {cout<<"l is not regular!"; return ;}; //irregular
			lcm=lcm*sum_De[i]/(GCD(lcm,sum_De[i]));
		};
		for (i=0;i<=d;i++)
		{
			total+=sum_Nu[i]*(lcm/sum_De[i]);
		};
		lcm=lcm*de;
		total=total*v;
		if (a==0) {a=total;b=lcm;}
		else {tem=b*lcm/GCD(b,lcm);a=a*tem/b+total*tem/lcm;b=tem;};	
		g=GCD(a,b);
		a=a/g;
		b=b/g;
	}
	if (b<0) {b=-b;a=-a;};
	cout<<"The desired integral is equal to:"<<a<<"/"<<b<<endl;
};	
