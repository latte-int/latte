/*******************************************************************
   Auther: Ruriko Yoshida
   July 24th, 2002
   Update: September 6th, 2002.
   This program computes Barvinok's decomposition of a cone.
   This program is for the project "LattE."

*********************************************************************/

#include <fstream.h>
#include <cstdlib>
#include <cstring>

#include "Cone.h"


int Test_Points(int level, RR (*Min_Max)[2], int dim, mat_RR *A_inverse, vec_RR *point);

ZZ lcm(const ZZ& a, const ZZ& b) { return a * ( b / GCD(a, b)); }

ZZ norm(const vec_ZZ& x, long m){
  ZZ normal;

  for(int i = 1; i <= m; i++)
    if(normal < abs(x(i)))    
       normal = abs(x(i));
  //cout << normal << endl;
  return normal;
}


RR norm2(const vec_RR& x, long m){
  RR  normal;

  for(int i = 1; i <= m; i++)
    if(normal < abs(x(i)))    
       normal = abs(x(i));
  //cout << normal << endl;
  return normal;
}

vec_ZZ ComputeOmega( const mat_ZZ & B, long m, int l, int y)
{
   mat_ZZ U; 
   mat_ZZ W;
   
   ZZ cm;
   set(cm);
 
   U.SetDims(m, m);
   W.SetDims(m, m);

  ZZ D = determinant(B);
  W = B;
  /*   for(int i = 1; i <= m; i++){
     for(int j = 1; j <= m; j++){
       if(W(i, j) != 0)
        cm = lcm( cm,  W(i, j));
        }
      }*/


  ZZ det2;

//     long lc;

   mat_RR R, R5;
   R.SetDims(m, m);
   R5.SetDims(m, m);
 
   for(int i = 1; i <= m; i++){
     for(int j = 1; j <= m; j++){
          conv(R(i, j), W(i, j));
          }
     }
  
   R5 = determinant(R) * inv(R);
   
   for(int i = 1; i <= m; i++){
     for(int j = 1; j <= m; j++){
            R5(i, j) = round(R5(i, j));
     }
   }

 
  mat_ZZ L, Z1;
   vec_ZZ Z;

   L.SetDims(m, m);
   Z1.SetDims(m, m);
   for(int i = 1; i <= m; i++){
     for(int j = 1; j <= m; j++){
            conv(L(i, j), R5(i, j));
     }
   }
   
   mat_ZZ LL;
   LL = L;
  LLL( det2, L, U, 1, 1);
  //cout << U << endl << L << endl;
  int index = 1;
  ZZ p, pp;
//    long tmp = 1000000;
   conv(pp, 1000000);
  for(int i = 1; i <= m; i++){
    // for(int j = 1; j <= m; j++){
         //p = abs(L(i, j)) * abs(L(i, j)) + p;}
          p = norm(L(i), m);
         if((p != 0) && (pp > p))
           { pp = p;
	   index = i;}

           }
  Z = U(index);
     
	if(l != 0)
	{
  		ZZ Array[m];
  		/*for(int i = 0; i < m; i++)
    		Array[i] = 0; */

  		Array[index - 1] = 1;
  		vec_ZZ ShortV;
  		ShortV = L(index);

  		ZZ norm1;
   		conv(norm1, 100000000);

  		ZZ count;
  		conv(count, 1);

   		vec_ZZ dummyV;
    		dummyV.SetLength(m);
  		//dummyV = L(index);
  		ZZ norm2;
 		// norm2 = norm(dummyV, m);
  		ZZ number[m];

  		/*  for(int i = 0; i < m; i++)
      		number[i] = 0; */
  		ZZ tmp3;
  		ZZ tmp, tmp2;
  		if(l == 2)
  		{
     			for(int i = 0; i < pow(10, m) - 1; i++)
			{
       				tmp = count;
       				for(int j = m - 1; j > -1; j--)
				{	
	 				tmp2 = tmp % 10;
         				number[j] = tmp2;
         				tmp = tmp / 10;
      	 			}
       				/*       for(int j = 0; j < m; j++)
         			cout << number[j] << " ";
	 			cout << endl;*/
       				for(int j = 0; j < m; j++)
        			{	 
					number[j] -= 5;
        				conv(dummyV(j + 1),0);
				}
	       			for(int j = 0; j < m; j++)
				{
       		  			conv(tmp3, number[j]);
		 			dummyV = dummyV + tmp3 * L(j + 1); //cout << dummyV << endl;
        			}
       
       				norm2 = norm(dummyV, m); // cout << norm2 << endl;
       				if((norm2 < norm1) && (IsZero(norm2) != 1))
				{
     		    			for(int j = 0; j < m; j++)
       			    			Array[j] = number[j];
         				
					ShortV = dummyV; 
         				norm1 = norm2;
      				}
              			count++;
     			}
     			// cout << ShortV << endl;

  			vec_ZZ INDEX;
  			INDEX.SetLength(m);

    			for(int i = 1; i <= m; i++)
       				INDEX(i) = Array[i-1]; 
    	
			// Z = INDEX * U;
    			LatticeSolve(Z, LL, ShortV);
    			//cout << INDEX << endl;
    			count = 0;
   			INDEX.kill();
  		}
  
   		if(l == 1)
		{
    			if(((index % m) + y) <= m)
       				Z = U((index % m) + y);
    			else if(((index % m) + y) > m)
       				Z = U(((index % m) + y) % m + 1);
   		}
 
  	dummyV.kill();
  
  	}
  
  U.kill();
  W.kill();
  R.kill();
  R5.kill();
  L.kill();
  Z1.kill();
  LL.kill();
  return Z;
 }


RR k_root(RR *a, int k, int iterations)
{
	RR x;


	x = 1;

	for(int i = 0; i < iterations; i++)
	{
		x = ( (k-1) * power(x, k) + (*a) )    /    ( k * power(x, k-1) );
	}

	return x;	

}



vec_ZZ ComputeOmega_2(mat_ZZ & B, long m)
{ 
   mat_ZZ W;
   
   W.SetDims(m, m);

  ZZ D = determinant(B);
  
  
  
  W = B;
  

   mat_RR R, R_inverse;
   R.SetDims(m, m);
   
   vec_RR point;

   point.SetLength(m);
   
 
   RR det;
   conv(det, D);

   // now det will have 1/abs(det^1/m)
   det = abs(det);
   cout << endl << "Newtons method called...";
   det = k_root(&det, m, 200000);
   cout << "done." << endl;

   cout << endl << det << endl;

   // sets R = W (i.e., B)

   for(int i = 1; i <= m; i++){
     for(int j = 1; j <= m; j++){
          conv(R(i, j), W(i, j));
          }
     }
  

   det = 1 / det;
   
    R = R * det;
    
	transpose (R,R);

	RR	Min_Max [m][2];

	RR	Min,Max;

	for (int i = 1; i <= m; i++)
	{
		Min = R(i)(1);
		Max = R(i)(1);
		for ( int q = 1; q <=m; q++)
		{
			if ( R(i)(q) < Min )
				Min = R(i)(q);
			else if ( R(i)(q) > Max)
				Max = R(i)(q);
					
		}

		Min = Min - 0.3;
		Max = Max + 0.3;

		if ( abs (Min) > Max)
			Min = -1 * Max;
		else
			Max = -1 * Min;
		
		trunc (Min,Min);
		trunc (Max,Max);
		Min_Max[i-1][0] = Min;
		Min_Max[i-1][1] = Max;
    	}
	
	
    	R_inverse = inv(R);
     	//det = determinant(R_inverse);
     	//cout << endl << "det of R_inverse:  " << det << endl;

	int result = 1;
     
   	if( Test_Points(0, Min_Max, m, &R_inverse, &point) == -1)   
    	{
		cout << "NOoooooo!  :(" << endl;
		result = -1;
    	}
	   
   
          

 cout << "Deleting stuff" << endl;
   
  W.kill();
  R.kill();
  R_inverse.kill();
 
vec_ZZ integral_point;
  
	integral_point.SetLength(m);
	
if( result == 1)
{
		for(int i = 1; i <= m; i++)
		{
			// order correct?
			//conv(point(i), integral_point(i));
			conv ( integral_point(i), point (i) );
		}
			
	       	
		cout << "yay!! :)";

}

else
	for(int i = 1; i <= m; i++)
		integral_point(i) = 0;
 
	cout << "Returning" << endl;	
  return integral_point;
}


int Test_Points(int level, RR (*Min_Max)[2], int dim, mat_RR *A_inverse, vec_RR *point)
{

	
	int non_zero = 0;

	if(level == dim)
	{

				
		for(int i = 1; i <= dim; i++)
			if((*point)(i) > 0.5 || (*point)(i) < -0.5)
				non_zero = 1;
		
		vec_RR coordinates;
		RR Norm;
		coordinates.SetLength(dim);

		coordinates = (*A_inverse) * (*point);
		
		Norm = norm2(coordinates, dim);

		
		if((Norm < 1.0001) && non_zero == 1)
		{

			//cout << "Test_Points: Success! here is the point: ";
			//for(int i = 1; i <= dim; i++)
			//	cout << (*point)(i) << " ";
			//cout << endl;
			return 1;
		
		}	
		else
			return -1;

	}
	
	else
	{
		RR current;

		//cout << "Test_Points:  Min = " << Min_Max[level][0] << "  Max = " << Min_Max[level][1] << endl;
		current = Min_Max[level][0];

		while(current <= ( Min_Max[level][1] + 0.5))
		{
			(*point)(level+1) = current;
			if(Test_Points(level + 1, Min_Max, dim, A_inverse, point) == 1)
				return 1;

			current = current + 1;

		}

		return -1;

	
	}

}
