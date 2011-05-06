/* SubspaceAvoidingDecomposition.cpp -- Barvinok decomposition that avoids a prescribed subspace
	       
   Copyright 2006 Matthias Koeppe
   Copyright 2009 Robert Hildebrand
   Copyright 2011 Christof Soeger

   This file is part of LattE.
   
   LattE is free software; you can redistribute it and/or modify it
   under the terms of the version 2 of the GNU General Public License
   as published by the Free Software Foundation.

   LattE is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with LattE; if not, write to the Free Software Foundation,
   Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
*/

#include <cassert>
#include <iostream>
#include <math.h>
#include "barvinok/SubspaceAvoidingDecomposition.h"

#include "genFunction/matrix_ops.h"



using namespace std;

/* The maximum norm. */
static ZZ
max_norm(const vec_ZZ& x)
{
  ZZ norm;
  int m = x.length();
  int i;
  for(i = 1; i <= m; i++) {
    ZZ absx = abs(x(i));
    if (norm < absx) norm = absx;
  }
  return norm;
}

//method from Matthias
vec_ZZ
ComputeShortVectorAvoidingSubspace_old(const mat_ZZ & B, const mat_ZZ &Dual)
{
  int dim = B.NumRows();
  assert(dim == B.NumCols());
   
  mat_ZZ L = -transpose(Dual);
  mat_ZZ U; 
  U.SetDims(dim, dim);

  ZZ det2;
  LLL(det2, L, U, 1, 1);

  // Find shortest basis vector that avoids the subspace.

  int i;
  ZZ least_max_norm;
  int least_max_norm_index;
  bool anyone = false;
  for(i = 1; i <= dim; i++) {
    if (U(i, dim) != 0) {
      ZZ n = max_norm(L(i));
      if (!anyone || n < least_max_norm) {
        least_max_norm = n;
        least_max_norm_index = i;
        anyone = true;
      }
    }
  }
  assert(anyone); // Full-dimensional, so we must have one vector that
		  // avoids the subspace.

  vec_ZZ result_U = U(least_max_norm_index);
  vec_ZZ result_L = L(least_max_norm_index);

//cerr<<"U = "<<U<<endl;
//cerr<<"L = "<<L<<endl;
//cerr<<"least max norm = "<<least_max_norm<<" at index "<<least_max_norm_index<<endl;

#if 0
  // Try to greedily reduce the norm of result_L by adding or subtracting basis
  // vectors that lie in the subspace.
  bool change;
  cerr << "L = " << result_L << endl;
  do {
    change = false;
    for(i = 1; i <= dim; i++) {
      if (U(i, dim) == 0) {
	vec_ZZ tL = result_L - L(i);
	ZZ tn = max_norm(tL);
	if (tn < least_max_norm) {
	  least_max_norm = tn;
	  result_L = tL;
	  cerr << "L = " << result_L << endl;
	  result_U -= U(i);
	  change = true;
	}
	else {
	  tL = result_L + L(i);
	  tn = max_norm(tL);
	  if (tn < least_max_norm) {
	    least_max_norm = tn;
	    result_L = tL;
	    cerr << "L = " << result_L << endl;
	    result_U += U(i);
	    change = true;
	  }
	}
      }
    }
  } while (change);
#endif
  return result_U;
}



//now methods from Robert

vec_ZZ ComputeShortVectorAvoidingSubspace(const mat_ZZ & B, const mat_ZZ &Dual,
            long m, int l, int y);
vec_ZZ ComputeShortVectorAvoidingSubspace_3(const mat_ZZ & B, const mat_ZZ &Dual,
            long m, int l, int y);
vec_ZZ ComputeShortVectorAvoidingSubspace_30(const mat_ZZ & B, const mat_ZZ &Dual,
            long m, int l, int y);


vec_ZZ
ComputeShortVectorAvoidingSubspace(const mat_ZZ & B, const mat_ZZ &Dual) {
  int dim = B.NumRows();
  assert(dim == B.NumCols());
  return ComputeShortVectorAvoidingSubspace(B, Dual, dim, 0, 0);
//  return ComputeShortVectorAvoidingSubspace_old(B, Dual);
}







static void
rescaleMat(mat_ZZ &P,ZZ &scaling, int m, int n){
ZZ d;	
d = scaling;
//	cerr << scaling << endl;

	for(int i=1; i<= n; i++){
		for(int j=1; j<=n; j++){
			GCD(d, (P(i))(j), d);
		}
	}
	
	div(scaling, scaling,d);
	for(int i=1; i<= n; i++){
		for(int j=1; j<=n; j++){
			div((P(i))(j), (P(i))(j), d);
		}
	}
}

static void
rescaleVec(vec_ZZ &P,ZZ &scaling, int n){
ZZ d;	
d = scaling;
//	cerr << scaling << endl;

	for(int i=1; i<= n; i++){
		GCD(d, P(i), d);
	}
	
	div(scaling, scaling,d);

	for(int i=1; i<= n; i++){
		div(P(i), P(i), d);
		
	}
}


static void
sub(vec_ZZ &x, ZZ &xDen, vec_ZZ a, ZZ d, vec_ZZ b, int m){
	xDen = d;
	mul(x,b,d);
	sub(x,a, x);
	rescaleVec(x,xDen, m);
		
}



ZZ norm2ZZ(const vec_ZZ& x, int m){
  ZZ normal;

  InnerProduct(normal, x,x);
  return normal;
}







vec_ZZ Blomer(vec_ZZ target, ZZ targetDen, mat_ZZ L,int m, int n, vec_ZZ q) {


	vec_ZZ closest;      // the vector that will be returned
	int fillClosest = 0; // changes to 1 after first iteration of the loop

	vec_ZZ guess;        //current guess to the current CVP problem
	vec_ZZ cbn;          // c*b_n, where c is the int (ZZ) that is being looped upon
	vec_ZZ guessDiff;    // represents the difference between the guess and the target
	vec_ZZ closeDiff;    // represents the difference between the current closest and the target

	ZZ lower;  //Lower bound on the for loop
	ZZ upper;  //Upper bound on the for loop
		

	if (n == 1){		//When n=1, we enumerate c values of q-1, q, and q+1 and find the closest.
				lower = -1; //about 1 in 5000 there is an anomoly that is fixed by this being -1
				upper =  0; //seems to work fine as zero, but might choose it to be 1 just in case

				add(lower, lower, q(1));	
				add(upper, upper, q(1));
		for(ZZ c = lower; c <= upper; add(c,c,1)){
			//Recursively call this part of the algorithm
				mul(guess, L(n), c);
								
			//Compare guess to other close vectors
				if(fillClosest == 0){ 
					closest = guess; 
					mul(closeDiff, closest, targetDen);
					sub(closeDiff, target, closeDiff);
					fillClosest = 1; 
				}
				
				mul(guessDiff, guess, targetDen);
				sub(guessDiff, target, guessDiff);
				
				if(norm2ZZ(closeDiff,m) > norm2ZZ(guessDiff,m)) { 
					closest = guess;
					closeDiff = guessDiff;
				}	
				
		}		
		
			//  ////  Original Code - doesn't seem to always work.					
			// Set closest = L(1)*c
			//cerr << "q(1):  "  << q(1) << endl;
			//mul(closest, L(1), q(1));  					
	} 


	else { 
		
		//Create A = first n-1 vectors of L
			mat_ZZ A;
			A.SetDims(n-1, m);

			for(int j=1; j <= n-1; j++){
				A(j) = L(j);
			}
		
		//Solve for Projection Matrix  P = transpose(A) inv(A transpose(A)) A;
			mat_ZZ P;	
			mul(P, A, transpose(A));
		
			mat_ZZ Pstep;
			ZZ detP;
		
			inv(detP,Pstep,P); 
			mul(Pstep, Pstep, A);   
			mul(P, transpose(A), Pstep);

			//Rescale
			//rescaleMat(P, detP, m,m);
				
		//Declare w
			vec_ZZ w;
			ZZ wDen;
					
		//Solve for boundaries in the for loop
			lower = -(floor(n/2)); //-1;
			upper = floor((n/2)+ 1);

			add(lower, lower, q(n));	
			add(upper, upper, q(n));
		
		//Declare several integers for deciding the proper search directions
			int count = 0; // used to determine when to switch directions (if necessary)
		 			 // ---- this could be done more elegantly, but needed to iterations to run before deciding and this just works 

			ZZ c = q(n);  //start the enumeration at q since the optimal c is likely q or close to it
			
			ZZ searchDirection;  //tells whether to increase or decrese c based on checking values at q, and q+1
			searchDirection = 1;

		while( c <= upper && c >= lower) {
			//Compute t - c*b_n 
				mul(cbn, L(n), c);
				sub(w, wDen, target,targetDen, cbn, m);
		
			//Project onto the space of the other lattice vectors
				mul(w,w,P);
				mul(wDen, detP, wDen);
				rescaleVec(w, wDen, m);

			//Recursively call this part of the algorithm
				
				guess = Blomer(w, wDen, A,m, n-1,q);
				
			//Add back c*b_n
				add(guess, guess, cbn);
				
			//Compare guess to other close vectors
				if(fillClosest == 0){ //gives the result closest an intial value - done this way to assign this within the loop
					closest = guess; 
					mul(closeDiff, closest, targetDen);
					sub(closeDiff, target, closeDiff); //solve for the 'difference' between the closest and the target
					fillClosest = 1; //set to one so that we don't enter this if statement again
				}
				
				mul(guessDiff, guess, targetDen);
				sub(guessDiff, target, guessDiff); //solve for the 'difference' between the guess and the target
			
				
				int isCloser = 0;
				if(norm2ZZ(closeDiff,m) > norm2ZZ(guessDiff,m)) { 
					closest = guess;
					closeDiff = guessDiff;
					isCloser = 1;
				}
				
				if(count ==1 && isCloser == 0){ //if q is closer than q+1, then search in the negative direction
					c = q(n); //set to q so that then we add below, we will be at q-1
					searchDirection = -1;
					
				}
				
			count++;
			add(c,c,searchDirection);
		} 
	}
return(closest);

}



static vec_ZZ CVP(mat_ZZ L,vec_ZZ target, ZZ targetDen, int m, int n){
//  According to Johannes Blomer's algorithm from 2000 

	mat_ZZ B;
	mat_ZZ dualB;
	mat_ZZ dualL;
	mat_ZZ U;


//  Set up lattice to be a dual HKZ basis 

	ZZ d;
	ZZ detB;
	ZZ detDual;

	//Set B = L^T
		transpose(B,L);

	//Set dualB = B(B^T B)^{-1}
		mul(dualB, L, B);
		inv(detB,dualB,dualB);
		mul(dualB, B, dualB);

	//Rescale
		rescaleMat(dualB,detB, m, n);


	//Reduce the dualB
		transpose(dualL, dualB);
		BKZ_FP(dualL,U,.99, 1000);
		transpose(dualB, dualL);


	//Solve for the dual of the dualB and store into B
		mul(B, transpose(dualB), dualB);
		inv(detDual,B,B);
		mul(B, dualB, B);
		

	//Rescale
		mul(B, B, detB);
		rescaleMat(B,detDual, m, n);
	
	//Scaling Check (should be an integer matrix)
		if(detDual != 1){
     			cerr<< "error: " << detDual << endl;  //THIS SHOULD BE 1
		}

	//Store L in original orientation
		transpose(L,B);

	//Re-order the lattice basis such that <b_i,b^*_j> = 1  iff i + j = n+1, and = 0 otherwise
		mat_ZZ tempL;
		tempL = L;
		for(int i=1; i<= n; i++){
			L(i) = tempL(n-i+1);
		}
		

//  Compute q where Bq = target 
		vec_ZZ q;
		ZZ r;

	//solves q*L = d*target where d = det(L);
		solve(d, q, L, target);  
		//cerr << "q   " << d << q << endl;

	//when i>=2, we will take the floor of the value q/d.
		for(int i=2; i <= m; i++){
			div(q(i), q(i), d);
		}
	//when i = 1 we want to round q(1) to the nearest integer
		DivRem(q(1), r, q(1), d);
		
		//void DivRem(ZZ& q, ZZ& r, const ZZ& a, const ZZ& b);
		// q = floor(a/b), r = a - b*q.
		
		mul(r,r,2);
		if(r < 0){mul(r,r,-1);}
		if(d < 0){mul(d,d,-1);}
		
		if(r > d){add(q(1),q(1),1);}

//Call Blomer recursion
   		vec_ZZ closer;	
      		closer = Blomer(target, targetDen,  L, m,  n, q);
return(closer);
}



vec_ZZ
SVP_primeCVP(int m, mat_ZZ L){  //m is the size of L

	mat_ZZ U;
	U.SetDims(m,m);

/* To solve SVP_prime we will run CVP([2^(j+1) b_1, ..., b_n], t_j = 2^j b_1) for j = 1, ..., floor(log_2(||b_1||/alpha))
alpha = dist(b_1, span(b_2, ..., b_n)).  In our special case, alpha =  b_1(n). */

	vec_ZZ shortVec;
	ZZ shortLength;
	shortLength = 10000000;
		
	vec_ZZ temp;
	ZZ tempLength;

	vec_ZZ target;  //taget for CVP
	mat_ZZ reducedL;  //Separate storage of lattive basis for computations

	ZZ bound;
	ZZ det2;
	long bound2;
	InnerProduct(bound, L(m), L(m));
	//divide(bound, bound, (L(m))(size));//this is currently wrong 
		//- should divide by the norm of the projection onto the orthogoal complement of the other vectors
	bound2 = NumBits(bound)/2;  //similar to taking log_2(bound)/2
	
	if(bound2 > 30){
		bound2 = 30;  //computing errors occurr if we try to raise 2^31
	}

	for(int j=0; j <= bound2 + 5; j++){
		mul(target, pow(2.0,j), L(m));
		reducedL = L;
		mul(reducedL(m), pow(2.0,j+1), reducedL(m));
		
				//LLL(det2, reducedL, U,1,1);  //Must reduce L before calling NearVector
				//BKZ_FP(reducedL,U);  //Alternative reduction method that is better and runs faster

				//NearVector(temp, reducedL, target);  //Uses the reduced basis and rounding to approximate a solution
		ZZ targetDen;
		targetDen = 1;
		temp = CVP(reducedL,target, targetDen, m, m);
		sub(temp , temp, target);  //Subtract to put vector back close to the origin

		tempLength = max_norm(temp);  //Tests conclude we should minimize the inf norm

		if(tempLength < shortLength){	
			shortLength = tempLength;
			shortVec = temp;
		}
		
	} 
	
return(shortVec);
}








vec_ZZ
SVP_prime(int m, mat_ZZ L){  //m is the size of L

	mat_ZZ U;
	U.SetDims(m,m);

/* To solve SVP_prime we will run CVP([2^(j+1) b_1, ..., b_n], t_j = 2^j b_1) for j = 1, ..., floor(log_2(||b_1||/alpha))
alpha = dist(b_1, span(b_2, ..., b_n)).  In our special case, alpha =  b_1(n). */

	vec_ZZ shortVec;
	shortVec = L(m);
	ZZ shortLength;
	shortLength = max_norm(shortVec);
	//shortLength = 10000000;
		
	vec_ZZ temp;
	ZZ tempLength;

	vec_ZZ target;  //taget for CVP
	mat_ZZ reducedL;  //Separate storage of lattive basis for computations

	ZZ bound;
	ZZ det2;
	long bound2;
	InnerProduct(bound, L(m), L(m));
	//divide(bound, bound, (L(m))(size));//this is currently wrong 
		//- should divide by the norm of the projection onto the orthogoal complement of the other vectors
	bound2 = NumBits(bound)/2;  //similar to taking log_2(bound)/2
	
	if(bound2 > 30){
		bound2 = 30;  //computing errors occurr if we try to raise 2^31
	}

	for(int j=0; j <= bound2 + 5; j++){
		mul(target, pow(2.0,j), L(m));
		reducedL = L;
		mul(reducedL(m), pow(2.0,j+1), reducedL(m));
		
		//LLL(det2, reducedL, U,1,1);  //Must reduce L before calling NearVector
		BKZ_FP(reducedL,U);  //Alternative reduction method that is better and runs faster

		NearVector(temp, reducedL, target);  //Uses the reduced basis and rounding to approximate a solution
		
		sub(temp , temp, target);  //Subtract to put vector back close to the origin

		tempLength = max_norm(temp);  //Tests conclude we should minimize the inf norm

		if(tempLength < shortLength){	
			shortLength = tempLength;
			shortVec = temp;
		}
		
	} 


 	//vec_ZZ Zvec;
 	//LatticeSolve(Zvec, originalL, shortVec);
 	//cerr << "Winner: " << max_norm(shortVec) << " , "<< Zvec << endl;
	
return(shortVec);
}




vec_ZZ ComputeShortVectorAvoidingSubspace(const mat_ZZ & B, const mat_ZZ &Dual,
		    long m, int l, int y)
{
  mat_ZZ U; 
  U.SetDims(m, m);

  ZZ D = determinant(B);
  if (D<0)  D = -D;

  ZZ det2;

  mat_ZZ L = -transpose(Dual);
  vec_ZZ Z;

  mat_ZZ LL;
  LL = L;
  
  
  //LLL( det2, L, U, 1, 1); // Can be swapped with BKZ(L,U,30);
  BKZ_FP(L,U);  
  int index = 1;
  ZZ p, pp;

   conv(pp, 1000000);
  for(int i = 1; i <= m; i++){ //Determine shortest vector in the reduced lattice basis
         p = max_norm(L(i));
         if((p != 0) && (pp > p))
           { pp = p;
	   index = i;}

           }
 Z = U(index);


ZZ lastTerm;
vec_ZZ Zother;

lastTerm = (U(index))(m);
//cerr << lastTerm << endl;

/* Decide to use SVP_prime  */
if (lastTerm == 0 && l==0) { // && max_norm(L(index)) != 1) {


	
	vec_ZZ shortVec = SVP_prime(m,LL);	
	//cerr << "lattice solving..." << endl;
	//cerr << shortVec << endl;
	LatticeSolve(Zother, LL, shortVec);
	//cerr << "done." << endl;
	if( D - max_norm(shortVec) <= 0){
	//cerr << "diff, norms:  " << D  << " , " << max_norm(shortVec) << " , " << max_norm(L(index))  << endl;
	}//cerr << "Lattice Solve: " << Zother << endl;
	 //cerr << "U(index): " << Z << endl;


	//need some criterion to use the avoidance vector
	//if(max_norm(shortVec) < max_norm(L(index)) + 10){
	//cerr << "Using SVP_prime vector" << endl;
	Z = Zother;
	
} 

if (lastTerm == 0 && l==1) { // && max_norm(L(index)) != 1) {


	
	vec_ZZ shortVec = SVP_primeCVP(m,LL);	
	LatticeSolve(Zother, LL, shortVec);
	if( D - max_norm(shortVec) <= 0){
	//cerr << "diff, norms:  " << D  << " , " << max_norm(shortVec) << " , " << max_norm(L(index))  << endl;
	}//cerr << "Lattice Solve: " << Zother << endl;
	 //cerr << "U(index): " << Z << endl;


	//need some criterion to use the avoidance vector
	//if(max_norm(shortVec) < max_norm(L(index)) + 10){
	//cerr << "Using CVP vector" << endl;
	Z = Zother;
	
}



  U.kill();
  L.kill();
  LL.kill();
  return Z;
 }

vec_ZZ ComputeShortVectorAvoidingSubspace_3(const mat_ZZ & B, const mat_ZZ &Dual,
		    long m, int l, int y)
{
  mat_ZZ U; 
   
  U.SetDims(m, m);

  ZZ D = determinant(B);

  ZZ det2;

  mat_ZZ L = -transpose(Dual);
  vec_ZZ Z;

  mat_ZZ LL;
  LL = L;
  
  
  //Tests done on hickerson-10  //num unimodular cones // run time (sec)
   //LLL( det2, L, U, 1, 1); // Can be swapped with BKZ(L,U,30);
	BKZ_FP(L,U);  
// LLL( det2, L, U, 2, 7); // 4410  35
  // LLL( det2, L, U, 3,7);// 4106  1.8
  // LLL( det2, L, U, 4,7);// 3291 .9
  // LLL( det2, L, U, 5, 7);// 3202 .9
  // LLL(  det2, L,U,6,7);  //3128  .85
  // LLL( det2, L,U,7,7);
  // BKZ_FP(L,U);  // 3108 .85
  // cerr << "BKZ computing" << endl;
   
  //BKZ(L,U,30);
  //BKZ(L,U,40);
  //BKZ(L,U,50);
  //cerr << U << endl << L << endl;
  int index = 1;
  ZZ p, pp;
//    long tmp = 1000000;
   conv(pp, 1000000);
  for(int i = 1; i <= m; i++){ //Determine shortest vector in the reduced lattice basis
    // for(int j = 1; j <= m; j++){
         //p = abs(L(i, j)) * abs(L(i, j)) + p;}
          p = max_norm(L(i));
         if((p != 0) && (pp > p))
           { pp = p;
	   index = i;}

           }
 Z = U(index);


	if(l != 0)
	{
	  vec_ZZ Array;
	  Array.SetLength(m);
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
 		// norm2 = max_norm(dummyV);
  		vec_ZZ number;
		number.SetLength(m);

  		/*  for(int i = 0; i < m; i++)
      		number[i] = 0; */
  		ZZ tmp3;
  		ZZ tmp, tmp2;
  		if(l == 2)
  		{
		  for(int i = 0; i < pow(10.0, (int)m) - 1; i++)
			{
       				tmp = count;
       				for(int j = m - 1; j > -1; j--)
				{	
	 				tmp2 = tmp % 10;
         				number[j] = tmp2;
         				tmp = tmp / 10;
      	 			}
       				/*       for(int j = 0; j < m; j++)
         			cerr << number[j] << " ";
	 			cerr << endl;*/
       				for(int j = 0; j < m; j++)
        			{	 
					number[j] -= 5;
        				conv(dummyV(j + 1),0);
				}
	       			for(int j = 0; j < m; j++)
				{
       		  			conv(tmp3, number[j]);
		 			dummyV = dummyV + tmp3 * L(j + 1); //cerr << dummyV << endl;
        			}
       
       				norm2 = max_norm(dummyV); // cerr << norm2 << endl;
       				if((norm2 < norm1) && (IsZero(norm2) != 1))
				{
     		    			for(int j = 0; j < m; j++)
       			    			Array[j] = number[j];
         				
					ShortV = dummyV; 
         				norm1 = norm2;
      				}
              			count++;
     			}
     			// cerr << ShortV << endl;

  			vec_ZZ INDEX;
  			INDEX.SetLength(m);

    			for(int i = 1; i <= m; i++)
       				INDEX(i) = Array[i-1]; 
    	
			// Z = INDEX * U;

			cerr<< "that's strange"<< endl;

			//cerr << INDEX << endl;
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
  L.kill();
  LL.kill();
  return Z;
 }






vec_ZZ ComputeShortVectorAvoidingSubspace_30(const mat_ZZ & B, const mat_ZZ &Dual,
		    long m, int l, int y)
{
  mat_ZZ U; 
   
  U.SetDims(m, m);

  ZZ D = determinant(B);

  ZZ det2;

  mat_ZZ L = -transpose(Dual);
  vec_ZZ Z;

  mat_ZZ LL;
  LL = L;
  
  
 
  LLL( det2, L, U, 1, 1); // Can be swapped with BKZ_FP(L,U,30);
  int index = 1;
  ZZ p, pp;
//    long tmp = 1000000;
  conv(pp, 1000000);
  for(int i = 1; i <= m; i++){ //Determine shortest vector in the reduced lattice basis

         
          p = max_norm(L(i));
         if((p != 0) && (pp > p))
           { pp = p;
	   index = i;}
           }

  ZZ lastTerm;
  lastTerm = (L(index))(m);
  vec_ZZ ShortV = L(index);

  if (lastTerm == 0) {
	ShortV = SVP_prime(m,L);	
  }


	
 		
  LatticeSolve(Z, LL, ShortV);
  U.kill();
  L.kill();
  LL.kill();
  return Z;
 }
