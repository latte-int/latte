#include <NTL/vec_ZZ.h>
#include <NTL/ZZ.h>
#include <NTL/vec_vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include <iostream>

using namespace std;

NTL_CLIENT

void print_Integrate(int d, int m, vec_ZZ l, vec_vec_ZZ s)
{
	ZZ v,de,sum,lcm,total,g;
	int i,j;
	vec_ZZ inner_Pro;
	vec_ZZ sum_Nu;
	vec_ZZ sum_De;
	inner_Pro.SetLength(d);
	sum_Nu.SetLength(d);
	sum_De.SetLength(d);
	lcm=1; //least common multiple
	total=0; //total numerator
	v=1;
	mat_ZZ matt;	
	matt.SetDims(d,d);
	for (i=1;i<=d;i++) matt[i-1]=s[i]-s[0];
	v=determinant(matt);if (v<0) v=-v;
	de=1;
	for (i=1;i<=d;i++)
	{
		de=de*(m+i);
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
	total=total*v;
	lcm=lcm*de;
	g=GCD(total,lcm);
	total=total/g;
	lcm=lcm/g;
	if (lcm<0) {lcm=-lcm;total=-total;};
	cout<<"The desired integral is equal to:"<<total<<"/"<<lcm<<endl;
};
	
	
	
	
