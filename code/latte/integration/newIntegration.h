#include "PolyRep.h"

#include <NTL/vec_ZZ.h>
#include <NTL/ZZ.h>
#include <NTL/vec_vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include <iostream>

using namespace std;

NTL_CLIENT

void convertToSimplex(int d, vec_vec_ZZ &s, ZZ &v, string line)
{
	int index,i,t,j;
	string temp,subtemp;
 	s.SetLength(d+1);
	index=1;
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
	v=1;
	mat_ZZ matt;	
	matt.SetDims(d,d);
	for (i=1;i<=d;i++) matt[i-1]=s[i]-s[0];
	v=determinant(matt);if (v<0) v=-v;
};

void update(ZZ &a, ZZ &b, vec_ZZ l, vec_vec_ZZ s,ZZ m, ZZ coe, ZZ v, int d, ZZ de)
{
	ZZ sum,lcm,total,g,tem;
	int i,j;
	vec_ZZ inner_Pro,sum_Nu,sum_De;
	inner_Pro.SetLength(d+1);
	sum_Nu.SetLength(d+1);
	sum_De.SetLength(d+1);
	total=0;
	lcm=1;
	for (i=0;i<=d;i++)
	{
		sum=0; for (j=0;j<d;j++) {sum=sum+l[j]*s[i][j];};
		inner_Pro[i]=sum;
	};//stores inner product for use
	for (i=0;i<=d;i++)
	{
		sum_Nu[i]=1;for (j=0;j<m+d;j++) sum_Nu[i]=sum_Nu[i]*inner_Pro[i];
		sum_De[i]=1;for (j=0;j<=d;j++) if (i!=j) sum_De[i]=sum_De[i]*(inner_Pro[i]-inner_Pro[j]);
		if (sum_De[i]==0) {cout<<l<<"is not regular!"; return ;}; //irregular
		lcm=lcm*sum_De[i]/(GCD(lcm,sum_De[i]));
	};
	for (i=0;i<=d;i++)
	{
		total+=sum_Nu[i]*(lcm/sum_De[i]);
	};
	lcm=lcm*de;
	total=total*v*coe;
	if (a==0) {a=total;b=lcm;}
	else {tem=b*lcm/GCD(b,lcm);a=a*tem/b+total*tem/lcm;b=tem;};	
	g=GCD(a,b);
	a=a/g;
	b=b/g;
}

void integrateList(string line, string line2)
{
	ZZ v,de,a,b,counter,m,tem,coe;
	int t,tt,i,j,index,k,d;
	vec_ZZ l;
	vec_vec_ZZ s;
	string temp,subtemp;
	t=2;d=1;
	t=line.find("[",t)+1;t=line.find("[",t)+1;
	temp=line.substr(t,line.find("]",t)-t);
	for (i=0;i<temp.length();i++) d+=(temp.at(i)==',');
	convertToSimplex(d,s,v,line2);//convert the string into a vec_vec_ZZ
	l.SetLength(d);
	index=2;
	a=0;
	b=0;
	while (index<line.length())
	{
		t=line.find("]]]",index);
		temp=line.substr(index,t-index);
		subtemp=temp.substr(0,temp.find(","));
		t=temp.find("[",0)+1;
		coe=to_ZZ(subtemp.c_str());
		subtemp=temp.substr(t,temp.find(",",t)-t);
		m=to_ZZ(subtemp.c_str());
		de=1;
		for (i=1;i<=d;i++)
		{
			de=de*(m+i);
		};
		t=temp.find("[",t)+1;
		temp=temp.substr(t,temp.length()-t+1);
		t=0;
		for (i=0;i<d-1;i++)
		{
			tt=temp.find(",",t);
			subtemp=temp.substr(t,tt-t);
			l[i]=to_ZZ(subtemp.c_str());
			t=tt+1;
		}
		subtemp=temp.substr(t,temp.length()-t+1);
		l[d-1]=to_ZZ(subtemp.c_str());
		update(a,b,l,s,m,coe,v,d,de);
		index=line.find("]]]",index)+5;
	};
	if (b<0) {b=-b;a=-a;};
	cout<<"The desired integral is equal to:"<<a<<"/"<<b<<endl;
};

void integrateFlatVector(int d, const linFormSum &lForm , string line)
{
  ZZ v,a,b,de,counter,m,tem,coe;
	int i,j,index,k;
	vec_ZZ l;
	vec_vec_ZZ s;
	convertToSimplex(d,s,v,line);//convert the string into a vec_vec_ZZ
	l.SetLength(d);
	lBlock* temForm=lForm.lHead;
	cBlock* temCoef=lForm.cHead;
	k=-1;
	a=0;
	b=0;
	cout<<"The polynomial is decomposed to "<<lForm.termCount<<" terms."<<endl;
	for (counter=0;counter<lForm.termCount;counter++)
	{
		k++;
		if ((k>0)&&(k % BLOCK_SIZE==0)) {temForm=temForm->next;temCoef=temCoef->next;k=0;}
		for (i=0;i<d;i++) l[i]=temForm->data[k][i];
		m=temForm->degree[k];
		coe=temCoef->data[k];
		de=1;
		for (i=1;i<=d+m;i++)
		{
			de=de*i;
		};
		update(a,b,l,s,m,coe,v,d,de);
	}
	if (b<0) {b=-b;a=-a;};
	cout<<"The desired integral is equal to:"<<a<<"/"<<b<<endl;
};

void integrateFlatVector(ZZ& numerator, ZZ& denominator, const linFormSum &forms , string line)
{
  ZZ v,de,counter,m,tem,coe;
	int i,j,index,k;
	int d = forms.varCount;
	vec_ZZ l;
	vec_vec_ZZ s;
	convertToSimplex(d,s,v,line);//convert the string into a vec_vec_ZZ
	l.SetLength(d);
	lBlock* temForm=forms.lHead;
	cBlock* temCoef=forms.cHead;
	k=-1;
	numerator=0;
	denominator=0;
	//cout<<"The polynomial is decomposed to "<<forms.termCount<<" terms."<<endl;
	for (counter=0;counter<forms.termCount;counter++)
	{
		k++;
		if ((k>0)&&(k % BLOCK_SIZE==0)) {temForm=temForm->next;temCoef=temCoef->next;k=0;}
		for (i=0;i<d;i++) l[i]=temForm->data[k][i];
		m=temForm->degree[k];
		coe=temCoef->data[k];
		de=1;
		for (i=1;i<=d+m;i++)
		{
			de=de*i;
		};
		update(numerator,denominator,l,s,m,coe,v,d,de);
	}
	if (denominator<0) {denominator *= to_ZZ(-1); numerator *= to_ZZ(-1);};
};	
