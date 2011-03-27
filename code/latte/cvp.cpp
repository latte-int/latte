/* SAP.cpp -- Shortest Avoiding Vector Program 

   Copyright 2009 Robert Hildebrand
   

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

// #include <string.h>
// #include <stdio.h>
// #include <cassert>

// #include "config.h"
// #include "latte_ntl_integer.h"
// #include "barvinok/dec.h"
// #include "barvinok/barvinok.h"
// #include "barvinok/Triangulation.h"
// #include "vertices/cdd.h"
// #include "genFunction/maple.h"
// #include "genFunction/piped.h"
// #include "cone.h"

// #include "dual.h"
// #include "RudyResNTL.h"
// #include "Residue.h"
// #include "Grobner.h"
// #include "preprocess.h"
// #include "print.h"
// #include "ramon.h"
// #include "rational.h"
// #include "timing.h"
// #include "flags.h"


// #include "banner.h"
// #include "convert.h"
// #include "latte_system.h"
// #include "Polyhedron.h"
// #include "ReadPolyhedron.h"
// #include "ProjectUp.h"

#include <stdio.h>
#include <stdlib.h>
#include "genFunction/matrix_ops.h"
#include <math.h>
#include <timing.h>
#include "latte_ntl.h"




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


ZZ norm(const vec_ZZ& x, int m){
  ZZ normal;

  for(int i = 1; i <= m; i++)
    if(normal < abs(x(i)))    
       normal = abs(x(i));
    return normal;
}

ZZ norm2(const vec_ZZ& x, int m){
  ZZ normal;

	InnerProduct(normal, x,x);
    return normal;
}







vec_ZZ Blomer(vec_ZZ target, ZZ targetDen, mat_ZZ L,int m, int n, vec_ZZ q) {
//Blomer is the recurion called within the CVP algorithm


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
				
				if(norm2(closeDiff,m) > norm2(guessDiff,m)) { 
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
				if(norm2(closeDiff,m) > norm2(guessDiff,m)) { 
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



static vec_ZZ   CVP(mat_ZZ L,vec_ZZ target, ZZ targetDen, int m, int n){
//  According to Johannes Blomer's algorithm from 2000 



// cerr << " numerators  ---  " << target.numerators() << endl;
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

	//account for targetDen
		//will want q*L/targetDen = d*target/targetDen, fix this by changing d
		mul(d,d,targetDen);

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
		
		if(r > d){add(q(1),q(1),1);
			//cerr << "yes " << endl;
		}

//Call Blomer recursion
   		vec_ZZ closer;	
      		closer = Blomer(target, targetDen,  L, m,  n, q);
return(closer);
}




static
void Random_Matrix(int size, int max, mat_ZZ &outL, int seed)
{
	mat_ZZ L;
	L.SetDims(size, size);
	double r;		// random value in range [0,1) 
	srand(seed);		// initialize random number generator
   
   	for(int i =1; i<= size; i++){
		for(int j=1; j<=size; j++){
       	r = (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
      	L(i,j) = (int)(r * max);// random integer in the range [0,max) 
			//{including 0, not including max}. If M is an integer then the range is [0,max-1]  
		}
	}
	
	outL = L;
}


/* ----------------------------------------------------------------- */


static ZZ  practiceSVP(mat_ZZ L, int size){

		//Timer::Timer(const std::string &a_name, bool start_timer) :
		//  name(a_name), ticks_elapsed(0), started(false)

	mat_ZZ  V , U;
		//L.SetDims(size, size);
		V.SetDims(size, size);
		U.SetDims(size, size);

	

	//cerr << "deterninant(L):  " << determinant(L) << endl;


	mat_ZZ originalL = L;
	mat_ZZ L2 = L;

	ZZ det3;
	LLL(det3,L2,U,1,1);

	vec_ZZ shorty;
	shorty = L2(1);

	ZZ shortyNorm = norm(shorty, size);
	for(int i = 2; i <= size; i++){
		if (shortyNorm > norm(L2(i), size)){
			shorty  = L2(i);
			shortyNorm = norm(shorty, size);
		}
	}




/* Polytime reduction to CVP */
	mat_ZZ B = transpose(L);

	mat_ZZ HB; /* SpaceAvoid = {x : Hx = 0}= {x: x_n = 0}  HB = H*B, 
		which is just B with the last row set to zero */

	HB.SetDims(size,size);
	HB(size) = B(size);
 


	SmithNormalForm(HB, V,U);
	V.kill();
	U.kill();
	HB.kill();

	L = transpose(B);
	B.kill();


//cerr <<" HB " <<endl << HB << endl;
//cerr <<" Unimodular Transformation " << endl << U << endl;
//cerr <<" New Lattice Basis " << endl << B << endl;
//cerr <<" first vector" << endl << L(1) << endl;


/* To solve SAP we will run CVP([2^(j+1) b_1, ..., b_n], t_j = 2^j b_1) for j = 1, ..., floor(log_2(||b_1||/alpha))
alpha = dist(b_1, span(b_2, ..., b_n)).  In our special case, alpha =  b_1(n). */

	int m = size;
	vec_ZZ shortVec = L(m);
	ZZ shortLength;

	InnerProduct(shortLength, shortVec, shortVec); //Solve for shortLength
	vec_ZZ temp;
	ZZ tempLength;
	vec_ZZ target;  //taget for CVP
	mat_ZZ reducedL;  //Separate storage of lattive basis for computations



	ZZ bound;
	ZZ det2;
	long bound2;
	InnerProduct(bound, L(m), L(m));
	divide(bound, bound, (L(m))(size));
	bound2 = NumBits(bound)/2;
			if(bound2 > 3){
				bound2 = 3;
			}
//Temporary
	bound2 = 0;

	vec_ZZ temp2;
	ZZ tempLength2;
	vec_ZZ shortVec2 = L(m);

	ZZ shortLength2;
	InnerProduct(shortLength2, shortVec2, shortVec2); //Solve for shortLength


		ZZ oneZZ;
		oneZZ = 1;
		rationalVector targetVec;
		mat_ZZ Lagain;

	for(int j=0; j <= bound2; j++){
		mul(target, pow(2.0,j), L(m));
		reducedL = L;
		
				
		mul(reducedL(m), pow(2.0,j+1), reducedL(m));
			
		Lagain = reducedL;
		
		//LLL(det2, reducedL, U,1,1);  //Must reduce L before calling NearVector
		BKZ_FP(reducedL,U);


		NearVector(temp, reducedL, target);  //Uses the reduced basis and rounding to approximate a solution
		

		ZZ targetDen;
		targetDen = 1;
		targetVec = rationalVector(target, oneZZ);
		temp2 = CVP(Lagain,target, targetDen, m, m);

		
	//	cerr << " temp"  << temp <<endl;
	//	cerr << " temp2" << temp2 <<endl;
		sub(temp , temp, target);  //Subtract to put vector back close to the origin
		 

		InnerProduct(tempLength, temp, temp);
		
		if(tempLength < shortLength){	
			shortLength = tempLength;
			shortVec = temp;
		}

		
				
		//targetVec = rationalVector(target, oneZZ);
		
		sub(temp2 , temp2, target);
		
		
		InnerProduct(tempLength2, temp2, temp2);
		
		if(tempLength2 < shortLength2){	
			shortLength2 = tempLength2;
			shortVec2 = temp2;
		}
				
		
		
	} 


//cerr << "Near Avoiding Vector" << endl << shortVec << endl;



//Double check the feasibility of the solution
//	vec_ZZ coeff;
//	solve(d, coeff, originalL, shortVec);
//	cerr << "Coefficients: " << coeff << d << endl;




 		vec_ZZ Zvec;
 		vec_ZZ Zvec2;

 		LatticeSolve(Zvec, originalL, shortVec);
 		LatticeSolve(Zvec2, originalL, shortVec2);

		ZZ z1;
		ZZ z2;
		ZZ z3;
		ZZ z4;
		ZZ z5;
		z1 = norm(shortVec,m);
		z2 = norm(shortVec2,m);
		z4 = norm2(shortVec,m);
		z5 = norm2(shortVec2,m);
		
		//if(z1 < z2){z3 = 1;}
		//if (z2< z1){z3 = 2;}
		//if(z1 == z2){z3 = 0;}

		if(z4 < z5){z3 = 1;}
		if (z5 < z4){z3 = 2;}
		if(z4 == z5){z3 = 0;}

		if(z4 < z5){
		
		cerr << "LLL_inf = " << z1 << "  CVP_inf =" << z2 << "  LLL_2 = " << z4 << "  CVP_2 ="<< z5 << endl;
		}
	//	cerr << "comparison:  " << z3 << endl;
		
	//	cerr << "New iteration" << endl;
	//	cerr << "LLL shortest vector" << shortVec << endl;
	//	cerr << "CVP shortest vector" << shortVec2 << endl;

	//	cerr << "LLL coefficients" << Zvec << endl;
	//	cerr << "CVP coefficients" << Zvec2 << endl;

return(z3);
}
/* ----------------------------------------------------------------- */




int main(){


cerr << "Time: " << GetTime() << " sec\n\n";

//  Regular Testing
		ZZ checks;
		ZZ addChecks;
		ZZ LLLshorter;
		ZZ CVPshorter;
		ZZ Bothequal;
		LLLshorter = 0;
		CVPshorter = 0;
		Bothequal = 0;
		ZZ one;
		one = 1;
		checks = 0;
		for(int i = 0; i <= 1000; i++){
		mat_ZZ L;
		Random_Matrix(5, 10, L, 1000 +i);
		if(determinant(L) !=0){
			addChecks = practiceSVP(L, 5);
			if(addChecks == 0){add(Bothequal, Bothequal, one);}
			if(addChecks != 0){add(checks, checks,one);}
			if(addChecks == 1){ cerr <<"i=" << i <<  " LLL was shorter " << endl;
				add(LLLshorter, LLLshorter,one);}
			if(addChecks == 2){ add(CVPshorter, CVPshorter, one); } //cerr <<"i="<< i << "CVP was shorter" << endl;}
		}
		
		//targetVec = rationalVector(temp, oneZZ);
		}
		cerr << "LLLshorter = "<<LLLshorter << "  CVPshorter = "<<CVPshorter <<"  SameLength= "<< Bothequal << "   checks= " <<checks << endl << endl << endl<< endl << endl;
cerr << "Time: " << GetTime() << " sec\n\n";



return(0);

}




