with(combinat):with(LinearAlgebra):with(linalg):with(numtheory):
# Let A=[A_1,..,A_{N+1}] be a list of  positive integers and t a variable.
# We compute the knapsack for A: the number of  integral solutions for A_1x_1+..+A_{N+1}x_{N+1}=t.
# The output of this as a step polynomial function  of t.
#
# This program is aimed at computing efficiently the highest k coefficients
#  of this step polynomial.
# 
# 
# 
# Main Functions
# 1. complete_knapsack(A,t). Computes the Ehrhart polynomial for the knapsack A1x1+..A_{N+1}x_{N+1}=t, t a variable or an integer. The output is a quasi polynomial function of t, expressed in function of t^k and the function frac(T)=T-floor(T);
# 
# 2. coeff_Nminusk_knapsack(A,t,k). Computes the N-k-term in the Ehrhart polynomial. The output is a periodic function of t;
#
#
#
# Example usage
#  L:=[1,2,3,4,5];
#  coeff4minus0:=coeff_Nminusk_knapsack(L, t, 0); #find t^4 term
#  coeff4minus3:=coeff_Nminusk_knapsack(L, t, 3); #find t^1 term
#  eval(subs({T=10, t=10, Frac=fractionalpart}, coeff4minus3)); #evaluate the t^1 term at t=T=10.
#
# License
#   This Maple script is part of the LattE integrale 1.7 package made available
#   under the GNU General Public License at http://www.math.ucdavis.edu/~latte/.
#	




############################
#Start of Library functions#
############################


#Global values:
USE_DUAL:=true:						#uses dual cone decomposition if true
MAX_DET_part_f_of_knapsack:=0:	#global variable that collects a statistic on the performance of the part_f_of_knapsack function



#  Basic programs: primitive_vector, ortho_basis, ComplementList, random_vector
# PRIMITIVE VECTOR
# Input: A :a vector with rational coordinates.
# Output: A vector with integral coordinates:
# Math: the primitive vector on the half line R^+A;
# Example: #primitive_vector([0,-1/2])->[0,-1];
primitive_vector:=proc(A) local d,n,g;
	d:=nops(A);
	n:=ilcm(seq(denom(A[i]),i=1..d));
	g:=igcd(seq(n*A[i],i=1..d));
	if g<>0 then
		[seq(n*A[i]/g,i=1..d)];
     else [seq(n*A[i],i=1..d)];
	fi;
	end:
#COMPLEMENT LIST 
#The output is the Complement  List, within the list [1,..,d]
ComplementList:=proc(K,d);
	RETURN([seq (`if` (member(i,K)=false, i, op({})),i=1..d)]);end:

random_vector:=proc(N,d) local R;
	R:=rand(N);
	[seq(R()+1,i=1..d)]:
end:

ortho_basis:=proc(d) local i,v;
	for i from 1 to d do
		v[i]:=[seq(0,j=1..i-1),1,seq(0,j=i+1..d)]
	od;
	[seq(v[j],j=1..d)];
end:
# Cone decomposition

#  Signed decomposition into unimodular cones
# A "simplicial cone" is a list of  d linearly independent  vectors in Z^d, sometimes assumed primitive. 
# short_vector(A)
# # Input:   A is a list of d linearly independent vectors.
# # Output: sho is a vector of dimension d.
short_vector:=proc(A) local n,base,i,sho;  
	n:=nops(A);
	base:=IntegerRelations[LLL](A);
	sho:=base[1];
	i:=1; 
	while i<=n-1 do
	    if max(seq(abs(sho[j]),j=1..n))<=max(seq(abs(base[i+1][j]),j=1..n))
	             then sho:=sho; else sho:=base[i+1];
	    fi;
	    i:=i+1;
	od;
	sho;
end:


# # sign_entries_vector(V)
# #  Input : vector V of dimension d.
# # Output:  L=[ Lplus,Lminus,Lzero] is a partition of [1..d] into three sublists,
# #               according to the signs of the entries of the vector V.
sign_entries_vector:=proc(V) local d,i,Lplus,Lminus,Lzero; 
       d:=nops(V); Lplus:=[]; Lminus:=[];Lzero:=[];

       for i from 1 to d do 
          if type(V[i],positive)      then Lplus:=[op(Lplus),i];
          elif type(V[i],negative) then Lminus:=[op(Lminus),i];
                                else Lzero:=[op(Lzero),i];
          fi;
       od;
	[Lplus,Lminus,Lzero];
end:

# # good_vector(G)
# # Input   G  is a  "simplicial cone"
# # Output consists of 2 elements: 
# #              V is a vector in Z^d. 
# #               L=[ Lplus,Lminus,Lzero] is a partition of [1..d] into three sublists,
# #               according to the signs of the entries of the vector V. in the basis G. 
good_vector:=proc(G) local n,A,Ainverse,B,sho,V,L;
       n:=nops(G);  
       #A:=Transpose(Matrix(G));      
       Ainverse:=inverse(matrix(G));
       B:=convert(Ainverse, listlist);
       #B:=[seq(convert(Ainverse[i,1..n],list),i=1..n)]; 
       sho:=short_vector(B); 
       V :=[seq(add(G[j][i]*sho[j],j=1..n),i=1..n)];
       L:= sign_entries_vector(sho);
	[V,L];
end:

# # signed_decomp(eps,G,v,L)
# # Input :  eps = 1 or -1
# #             G  is a  "simplicial cone"
# #              V is a vector of dim d
# #              L= [ Lplus,Lminus,Lzero] is a partition of [1..d] into three sublists,
# # Output : [Nonuni,Uni] 
# #              Nonuni and Uni are  lists of terms  [eps,detG,G],  where
# #               eps=1 or -1, 
# #               detG is an integer,  
# #               G  is a  list of  d linearly independant primitive  vectors in Z^d. 
signed_decomp:=proc(eps,G,v,L) local Nonuni,Uni,Lplus,Lminus,Lzero,kplus,kminus,kzero,i,j, C,M, detC, Csigned ; 
	Nonuni:=[]; Uni:=[];
	Lplus:=L[1]; Lminus:=L[2]; Lzero:=L[3];
	kplus:=nops(Lplus); kminus:=nops(Lminus); kzero:=nops(Lzero);
	if kplus>0 then
		for i from 1 to kplus do
        		C:=[seq(G[Lplus[j]],j=1..i-1),seq(-G[Lplus[j]],j=i+1..kplus),v,seq(G[Lminus[j]],j=1..kminus),seq(G[Lzero[j]],j=1..kzero)];

		        detC := det(matrix(C));        
		        Csigned:=[eps*(-1)^(i+kplus),detC,C];       
 
		       if abs(detC)>1 then
       			Nonuni:=[op(Nonuni),Csigned] else Uni:=[op(Uni),Csigned];
        		fi;
   		od;
 	fi;

	if kminus>0 then
		for i from 1 to kminus do
       		C:=[seq(G[Lplus[j]],j=1..kplus),-v,seq(-G[Lminus[j]],j=1..i-1),seq(G[Lminus[j]],j=i+1..kminus),seq(G[Lzero[j]],j=1..kzero)]; 
        
       		detC := det(matrix(C));
       		Csigned:=[eps*(-1)^(i+1),detC,C];      
       
      			if abs(detC)>1 then
            			Nonuni:=[op(Nonuni),Csigned] else Uni:=[op(Uni), Csigned];
         		fi;
    		od;
 	end if;
 	[Nonuni,Uni];
 end:
 
# # good_cone_dec(eps,G)
# #  Input: eps = 1 or -1
# #             G  is a  simplicial cone
# #  Output:  two lists [Nonuni,Uni] as in procedure signed_decomp: 
good_cone_dec:=proc(eps,G) local n,A,R,Output, det_A;
	n:=nops(G);  A:=matrix([seq(G[i],i=1..n)]);
	det_A:=det(A);   
	if abs(det_A)=1 then
		Output:=[[],[[eps,det_A,G]]];
     	else R:=good_vector(G);
          	Output:=signed_decomp(eps,G,R[1],R[2]);
   	fi;
end:

# # more_decomposition_in_cones(cones)
# # Input:  cones =[cones[1],cones[2]] as in procedure signed_decomp
# # Output: [Newnonuni,Newuni] as in procedure signed_decomp
more_decomposition_in_cones:=proc(cones) local i,Newuni,Newnonuni,newcones:
	Newnonuni:=[]; 
	Newuni:=cones[2];
   	for i from 1 to nops(cones[1]) do
    		newcones:=good_cone_dec(cones[1][i][1],cones[1][i][3]);
   		Newnonuni:=[op(Newnonuni),op(newcones[1])];
   		Newuni:=[op(Newuni),op(newcones[2])];
 	od;
	[Newnonuni,Newuni];
end:  
     ###############################

# # cone_dec(G)
# # Input:  G is a "simplicial cone"
# # Output: A list of  terms [eps,detG,G] where
# #               eps =1 or -1, 
# #               detG is an integer ( hopefully 1 or -1),  
# #               G  is a  "simplicial cone", (hopefully unimodular)
cone_dec:=proc(G) local seed, i,ok;
	if G=[] then 
		RETURN([[1,1,[]]]);
	fi:
	seed:=good_cone_dec(1,G);
 	ok:=0;
	i:=1; 
	while ok=0  do
		seed:=more_decomposition_in_cones(seed); 
		if seed[1]=[] then
       		ok:=1;else ok:=0;i:=i+1;
     		fi;
	od;
    	RETURN(seed[2]);
end:

# new Cone dec via polar

#  Signed decomposition into unimodular cones via polar dec
# A "simplicial cone" is a list of  d linearly independent  vectors in Z^d, sometimes assumed primitive. 
# short_vector(A)
# # Input:   A is a list of d linearly independent vectors.
# # Output: sho is a vector of dimension d.
new_short_vector:=proc(A) local n,base,i,sho;  
	n:=nops(A);
	base:=IntegerRelations[LLL](A);
	sho:=base[1];
	i:=1; 
	while i<=n-1 do
	    if max(seq(abs(sho[j]),j=1..n))<=max(seq(abs(base[i+1][j]),j=1..n))
	             then sho:=sho; else sho:=base[i+1];
	    fi;
	    i:=i+1;
	od;
	sho;
end:



# # good_vector(G)
# # Input   G  is a  "simplicial cone"
# # Output is  a vector V in Z^d. 
# #                

polar_good_vector:=proc(G) local n,A,Ainverse,B,sho,V,L,newV,VV,ok,i,out;
          n:=nops(G);  
          A:=Transpose(Matrix(G));      
          Ainverse:=MatrixInverse(A);
          B:=[seq(convert(Ainverse[1..n,i],list),i=1..n)]; #print("B",B);
          sho:=new_short_vector(B); 
          V :=[seq(add(G[j][i]*sho[j],j=1..n),i=1..n)];
 newV:=LinearSolve(A,Vector(V));
VV:=convert(newV,list);
ok:=0;i:=1; while i<= nops(VV) and ok=0 do  
 if VV[i]<=0 then i:=i+1;
else ok:=1; i:=nops(VV)+1;
fi; od:   
if ok=1 then out:=V; else out:=[seq(-V[i],i=1..nops(V))];
fi;
out;
end:

# # signed_decomp(G,v)
# # Input :  signum = 1 or -1 
# #             G  is a  "simplicial cone"
# #              v is a vector of dim d            
# # Output : [Nonuni,Uni] 
# #              Nonuni and Uni are  lists of terms  [eps,G],  where
# #               eps=1 or -1,  
# #               G  is a  list of  d linearly independant primitive  vectors in Z^d. 
 

polar_signed_decomp:=proc(signum,G,v) local eps,Nonuni,Uni,Lplus,Lminus,Lzero,kplus,kminus,kzero,i,j, C,M, detC, Csigned ; 
	Nonuni:=[]; Uni:=[];
	
  eps:=sign(Determinant(Matrix(G)));
		for i from 1 to nops(G) do
        		C:=[seq(G[j],j=1..i-1),v,seq(G[j],j=i+1..nops(G))];

		        detC := Determinant(Matrix(C));        
		        Csigned:=[signum*eps*sign(detC),C];       
  if detC<>0then
		       if abs(detC)>1 then
       			Nonuni:=[op(Nonuni),Csigned] 
          
else Uni:=[op(Uni),Csigned];
        		fi;else fi;
   		od;
 	[Nonuni,Uni];
 end:
# # polar_good_cone_dec(signum,G)
# #  Input: signum = 1 or -1
# #             G  is a  simplicial cone
# #  Output:  two lists [Nonuni,Uni] as in procedure polar_signed_decomp: 
polar_good_cone_dec:=proc(signum,G) local eps,n,A,R,Output, det_A;
	n:=nops(G);  A:=Matrix([seq(G[i],i=1..n)]);
	det_A:=Determinant(A); eps:=sign(det_A);
       #print("G",det_A);
	if abs(det_A)=1 then #print(itisuni);
		Output:=[[],[[1,G]]];
      
     else R:=polar_good_vector(G);#print("goodvector",R,G);
          	Output:=polar_signed_decomp(signum,G,R);
   	fi;
end:

# # polar_more_decomposition_in_cones(cones)
# # Input:  cones =[cones[1],cones[2]] as in procedure polar_signed_decomp
# # Output: [Newnonuni,Newuni] as in procedure signed_decomp
polar_more_decomposition_in_cones:=proc(cones) local i,Newuni,Newnonuni,newcones:
	Newnonuni:=[]; #print("cones",cones);
	Newuni:=cones[2];
   	for i from 1 to nops(cones[1]) do
    		newcones:=polar_good_cone_dec(cones[1][i][1],cones[1][i][2]);

   		Newnonuni:=[op(Newnonuni),op(newcones[1])];
   		Newuni:=[op(Newuni),op(newcones[2])];
 	od;
	[Newnonuni,Newuni];
end:  

     ###############################

# # polar_cone_dec(G)
# # Input:  G is a "simplicial cone"
# # Output: A list of  terms [eps,G] where
# #               eps =1 or -1, 
# #               detG is an integer ( hopefully 1 or -1),  
# #               G  is a  "simplicial cone", (hopefully unimodular)
polar_cone_dec:=proc(G) local A,eps,seed, i,ok;
	if G=[] then 
		RETURN([[1,[]]]);
	fi:  A:=Matrix([seq(G[i],i=1..nops(G))]);
	eps:=sign(Determinant(A));#print("DET",eps);
	if abs(Determinant(A))=1 then 
		#seed:=polar_good_cone_dec(G);#print(absol=1,seed);
		RETURN([[1,G]]);
	else
	seed:=polar_good_cone_dec(1,G);fi;
 	#print("start",seed);
	ok:=0;
	i:=1; 
	while ok=0  do
		seed:=polar_more_decomposition_in_cones(seed); #print("i,SEED",i,seed);
		if seed[1]=[] then 
			ok:=1;else ok:=0;i:=i+1;
		fi;
	od;
	RETURN(seed[2]);
end:


# Cone decomposition via polar cone
primitive_vector:=proc(A) local d,n,g;
	d:=nops(A);
	n:=ilcm(seq(denom(A[i]),i=1..d));
	g:=igcd(seq(n*A[i],i=1..d));if g<>0 then
	[seq(n*A[i]/g,i=1..d)];else [seq(n*A[i],i=1..d)];fi;
end:

#The polar cone.
polarcone:=proc(A) local n,B,Ainverse,out;
	n:=nops(A);
	#MatrixInverse(Transpose(Matrix(A)));
	B:=Matrix(A); #print("B",Matrix(A),B);     
   Ainverse:=MatrixInverse(B);#print(Ainverse);
	out:=[seq(convert(Transpose(Column(Ainverse,i)),list),i=1..n)];
end:

#polarcone([[1,0],[1,1]]);
#Compute the cone decomposition going to polar cone and back
new_dual_dec_cone:=proc(G) 
	local out,dual_cone,L,dec,i;
	out:=[];
	dual_cone:=polarcone(G);#print("dual",dual_cone);
	L:=[seq(primitive_vector(dual_cone[m]),m=1..nops(dual_cone))];#print("L",L);
	dec:=polar_cone_dec(L);#print("decpolar",dec,L);
	#print("enteringpolarcone",dec);
	for i from 1 to nops(dec) do
		out:=[op(out),[dec[i][1],polarcone(dec[i][2])]];
	od:
	out:
end:

#  HERE STARTS BACK TO KNAPSACK;
# Let  W be a cone basis w_i in R^k. Lambda a lattice with base alpha_i in R^k, s=[s_1..s_k] symbolic in R^k.
# Compute coordinates and change of coordinates: cone_in_basis_Lambda; Basis_Lattice, the_vertex.
# INPUT : W:=a list of   p-vectors w_i of lenght r (with rational entries), Lambda a list of p-vectors alpha_i of lenght r (with rational entries). 
# OUTPUT: a list of r numvectors with integral coordinates.
# MATH: we look for the primitive vectors of the cone W generated by w_i in the basis of alpha_i
# EXAMPLE: Lambda:=[[2/3,-5/2],[3/7,-8/5]];W:=[[1,0],[1,1]];cone_in_basis_Lambda(W,Lambda):=[[-16, 25], [-426, 665]];Then -426*Lambda[1]+665*Lambda[2]=[1,1]=w[2]
cone_in_basis_Lambda:=proc(W,Lambda) 
	local M,MM,newW,i;
	newW:=[];
	M:=Transpose(Matrix([seq(Lambda[i],i=1..nops(Lambda))])):
	#print("M", M);
	MM:=MatrixInverse(M);
	#print("MM", MM);
	for i from 1 to nops(W) do
		newW:=[op(newW),primitive_vector(convert(Transpose(Multiply(MM,Vector(W[i]))),list))];
	od;
	newW;
end
:
# INPUT: Small :=[a1,a2,...,ar] a list of  positive integers, f an integer.
# OUTPUT: a list of r lists. Each of the list is a list of integers.
# MATH: we compute a basis of the lattice Lambda, consisting of the elements 
# [k1,k2,...,kr] such that k1*a1+k2*a2+...kr*ar is a multiple of f
# Example: Basis_Lattice([1,2,2],2):=[[2, 0, 0], [-2, 1, 0], [-2, 0, 1]];
Basis_Lattice:=proc(Small,f) local P,H,U,Lambda,hh,newf;
	P:=matrix([[seq(Small[i],i=1..nops(Small))]]);
	H:=ihermite(transpose(P),u); 
	hh:=row(H,1)[1]; 
	U:=eval(u):
	#print("U", U);
	newf:=f/igcd(hh,f);
	Lambda:=[[seq(newf*convert(row(U,1),list)[i],i=1..nops(Small))],seq([seq(convert(row(U,i),list)[k],k=1..nops(Small))],i=2..nops(Small))];
	Lambda; 
end:


# INPUT: Small :=[a1,a2,...,ar] a list of  positive integers, f an integer.
# OUTPUT: a list of r integers.
# MATH: We compute s1,....sr,s0 such that a1s1+a2s2+..arsr+s0f=1. we output only [s1..sr].
# Example: the_vertex([1,2,2],2):=[1, 0, 0];
the_vertex:=proc(Small,f) local P,H,out,u; 
	P:=matrix([[seq(Small[i],i=1..nops(Small)),f]]);
	H:=ihermite(transpose(P),u);  
	out:=[seq(convert(row(u,1),list)[i],i=1..nops(Small))]; 
	out;
end:

# Some tools: B_and_S; F_k;
# INPUT : A:=a list of  positive integers, f a positive integer.
# OUTPUT: a list of two lists: [B,S].
# MATH: B is the list of elements in A consising of the elements A[i] such that f divides A[i]. S is the complement list
# Example: B_and_S([2,2,2,3,4,5],2):=[[2, 2, 2, 4], [3, 5]];
# 
B_and_S:=proc(A,f) local j,B,S;
	B:=[];
	S:=[];
	for j from 1 to nops(A) do
		if modp(A[j],f)=0 then 
			 B:=[op(B),A[j]];
			 S:=S;
		 else
		 	B:=B; S:=[op(S),A[j]];
		 fi;
	od;
	[B,S];
end:


# find all k subsets of a list of objects
# @parm inputList: a list. ex [x, y, z, w]
# @parm k: the subset size
myChoose:=proc(inputList, k)
	local n, L;
	local i;
	local counter;
	local limit;
	local outputList:=[];
	local newSubset;
	local numDiff;

	#use the maple choose funciton if there are not many unique numbers.
	sortedInputList = sort(inputList);
	numDiff := 0;
	for i from 1 to nops(sortedInputList)-1 do
		if sortedInputList[i] <> sortedInputList[i+1] then
			numDiff := numDiff + 1;
		end;
	end;

	if numDiff = 2 then
		return choose(inputList, k);
	end;


	#L is going to be a list of index values to add to the output list.
	#L starts off as [n-k, ..., 3, 2, 1].
	#We keep adding values to L[1] and if L[i] > n - i +1 we add 1 to the next element
	#Hence the last element added is [n, n-1, n-2, ..., n-k+1]

	n:=nops(inputList);
	L:=[seq(0, i=1..k)];
	for i from 1 to k do
		L[i]:=k-i+1;
	end;
	
	limit:= binomial(n,k);
	counter:=1;
	while counter < limit do
	
		#add a subset of inputList to outputList indexed by L
		#print(counter, L);
		newSubset:=[seq(0,i=1..k)];
		for i from 1 to k do
			newSubset[i]:=inputList[L[k-i+1]];
		end;
		outputList:=[op(outputList), newSubset];
	
		L[1]:=L[1]+1;
		
		i:=1;
		while L[i] > (n - i + 1) do
			L[i+1]:=L[i+1]+1;
			i:=i+1;
		end;
		
		while i > 1 do
			L[i-1]:=L[i]+1;
			i:=i-1;	
		end;
		counter:=counter+1;
	end;
	
	#add the last element to the output list.
	#print(counter, L);
	newSubset:=[seq(0,i=1..k)];
	for i from 1 to k do
		newSubset[i]:=inputList[L[k-i+1]];
	end;
	outputList:=[op(outputList), newSubset];
end:





# Input:= a list A of positive integers, an integer k, 0<=k<=nops(A)
# Output:= a list
# Math: We compute S_k (Poles>=N+1-k); in the application  k is big
# F_k([2,3,3,55],2):=[1,3]; F_k([2,3,3,55],3):=[1, 2, 3, 55];
F_k:=proc(A,k) local N,g,C,i,out,a;
	out:={};
	N:=nops(A)-1;
	for a from 0 to k do
		C:=myChoose(A,N+1-a);
		for i from 1 to nops(C) do 
			g:= igcd(op(C[i]));out:={op(out),g};
		od:
	od:
	out:=[op(out)];
end:

# The coefficient of t^h in the function E(A,f)(t):
# 
#smallstep:=proc(a,n) local out;
#if  a=0 then out:=0;
#elif type(simplify(a),integer)=true  then out:=0;
#else out:=Smallstep(a*n);
#fi;
#end:
fract:=proc(a,n) local out,p,q,newa,t; 
	if  a=0 then out:=0; 
	elif type(simplify(a),integer)=true  then out:=0;
	else 
		p:=numer(a); q:=denom(a); 
		newa:=modp(p,q);
		t:=Frac(newa/q*n);
	fi;
end:

smallstep:=proc(a,n); fract(-a,n); end:

Todd:=proc(x,order):
-add(1/j!*bernoulli(j,0)*x^j,j=0..order):
end:

# Input: U is a  cone, Lambda is a lattice,
# Here: U is given in coordinates with respect to the  lattice basis,
# s is a pvector (with symbolic entries);
# xi is  a dvector (with symbolic entries), order a positive integer.
# Output:  a list of a nonnegative integer and a function of s,xi, with entries formal functions rest().
# Math: we compute the function e^<(smallstep(-T*s)),xi> prod_{i=1}^r 1/(1-e^(xi,g_j)), taking the Laurent series up to order and an nonnegative integer. 
# Here the g_j are the p-vectors generators of the cone U, it will be unimodular in the application.
# Later  we evaluate it at xi=[(a1+epsilon reg_1)*x,(a2+epsilon reg_2)x]...compute the Laurent series in x up to some length (later to be N+2-h);
# and we take the coefficient in epsilon=0. In view of this applicatiuon we evaluate  the the nonnegative integer which is the order of the pole in the variable epsilon
# s is written as sum_i x_i g_i and rest(s) is the p vector sum_i (ceil(x_i)-x_i)
# We also will use the formal functions rest(s)=ceil(s)-s;
# Example:=Function_M_uni_expanded([T*s1,T*s2],[[1,2],[1,3]],[[1,0],[0,1]],[2,-1],[xi1,xi2],3):=[1, (1+xi1*(Smallstep((-3*s1+s2)*T)+Smallstep((2*s1-s2)*T))+xi2*(2*Smallstep((-3*s1+s2)*T)+3*Smallstep((2*s1-s2)*T))+(1/2)*(xi1*(Smallstep((-3*s1+s2)*T)+Smallstep((2*s1-s2)*T))+xi2*(2*Smallstep((-3*s1+s2)*T)+3*Smallstep((2*s1-s2)*T)))^2+(1/6)*(xi1*(Smallstep((-3*s1+s2)*T)+Smallstep((2*s1-s2)*T))+xi2*(2*Smallstep((-3*s1+s2)*T)+3*Smallstep((2*s1-s2)*T)))^3)*(-1+(1/2)*xi1+xi2-(1/12)*(xi1+2*xi2)^2)*(-1+(1/2)*xi1+(3/2)*xi2-(1/12)*(xi1+3*xi2)^2)/((xi1+2*xi2)*(xi1+3*xi2))];
Function_M_uni_expanded:=proc(T,s,U,Lambda,Small,xi,order) local r, M1,M2,M,pro,QQ,Rests,i,a,news,b,FM,Eb,pole,coex:
	pro:=1; 
	r:=nops(U);
	M1:=Transpose(Matrix([seq(Lambda[i],i=1..nops(Lambda))]));
	M2:=Transpose(Matrix([seq(U[i],i=1..nops(Lambda))]));
	M:=Multiply(M1,M2);
	QQ:=[seq(convert(col(M,k),list),k=1..r)];##print("QQ",QQ);
	news:=convert(Multiply(MatrixInverse(M),Vector(s)),list); #small_step:=[seq(smallmove(-news[j],T),j=1..nops(news))]
	#Rests:=[seq(smallstep(news[j]),j=1..nops(news))];pole:=0;
	Rests:=[seq(smallstep(-news[j],T),j=1..nops(news))];pole:=0; 
	for i from 1 to nops(U) do
		coex:=add(Small[j]*QQ[i][j],j=1..nops(QQ[i])); 
		if coex=0 then 
			pole:=pole+1; 
		fi;
		a:=add(xi[j]*QQ[i][j],j=1..nops(QQ[i]));
		pro:=pro*1/a*Todd(a,order): ##print(i,pro);
	od;
	b:=add(xi[i]*(add(Rests[k]*QQ[k][i],k=1..nops(QQ))),i=1..r);
	Eb:=add(b^m/m!,m=0..order); ##print('Eb',Eb);
	[pole,Eb*pro]; 
end:



#  Input: AA list of positive integers ; f a positive integer,
# x a variable; T a symbolic variable, f a nonnegative integer, f<>1;
# Output; a function of x,T;
# Math: the coefficient of the function H(AA,f)(T,x) in x^(-h-1), that is wecompute the coefficient of t^h of E(A,f)(t).  
# Example:=part_f_of_knapsack([2,2,3,3],T,3,1):=17/48+(1/6)*Smallstep((1/3)*T)-(1/2)*Smallstep((1/3)*T)^2;

part_f_of_knapsack:=proc(AA,T,f,h) 
	local N,BB,Small,LAMBDA,fB,i,WW,newW,uni_cones,Tvertex,out,xxi,reg,FF,FM,FMfB,pole,newWmatrix,veryNewWmatrix,tempLambda,t1,j;
	global USE_DUAL;
	global MAX_DET_part_f_of_knapsack;
	N:=nops(AA)-1;
	BB:=B_and_S(AA,f)[1]; 
	Small:=B_and_S(AA,f)[2];
	reg:=random_vector(500,nops(Small));
	fB:=1: 
	for i from 1 to nops(BB) do #print(i,nops    (BB),Todd(BB[i]*x,N-h));
		fB:=fB*1/(BB[i]*x)*Todd(BB[i]*x,N-h);
	od; ##print(aqui);
	LAMBDA:=Basis_Lattice(Small,f); 
	WW:=ortho_basis(nops(Small));
	#print("WW=", WW);
	newW:=cone_in_basis_Lambda(WW,LAMBDA); #print("thecone",newW);

	#start  of new ideas for the dual cone computation
	newWmatrix:=Matrix(newW); #convert the list to a matrix.
	veryNewWmatrix:=MatrixInverse(newWmatrix); #this matrix will hold the dual cone rays.

	#find scaling for the new matrix
	for i from 1 to nops(newW) do
		tempLambda:=ilcm(seq(denom(veryNewWmatrix[j, i]),j=1..nops(newW)));

	  for j from 1 to nops(newW) do 
   	 veryNewWmatrix[j,i]:=veryNewWmatrix[j,i]*tempLambda;
	  od;# for jth row.
	od; #for ith col

	#print("org. matrix", Transpose(newWmatrix));
	#print("org. det= ", Determinant(newWmatrix));
	#print("new matrix", veryNewWmatrix);
	#print("new det=", Determinant(veryNewWmatrix));
	#print("f=",f);
	if USE_DUAL = false then
		MAX_DET_part_f_of_knapsack:=max(	MAX_DET_part_f_of_knapsack, abs(Determinant(newWmatrix)));
	else
		MAX_DET_part_f_of_knapsack:=max(	MAX_DET_part_f_of_knapsack, abs(Determinant(veryNewWmatrix)));
	fi;
	

	#end of new ideas

	#Tvertex:=[seq(the_vertex(-T*Small,f)[i],i=1..nops(Small))];
	Tvertex:=[seq(the_vertex(Small,f)[i],i=1..nops(Small))];
	 ##print(aqui,Tvertex);
	
	t1:=time();

	#here, branch on if we are using the dual cones or not.
	if USE_DUAL = false then 
	  #print("using primal cones");
	  uni_cones:=cone_dec(newW); #
	else
	  	#print("using dual cones");
		#print("dual cone computation is not correct");
		#todo: we have to dualize back.
		#uni_cones:=cone_dec(convert(Transpose(veryNewWmatrix), listlist)); 
		uni_cones:=new_dual_dec_cone(newW);
	fi;
	
	#print("total time in uni_cones ", time() -t1);
	 #print(nopsunicones,nops(uni_cones),uni_cones);
	out:=0;
	xxi:=[seq((Small[z]+epsilon*reg[z])*x,z=1..nops(Small))];
	for i from 1 to  nops(uni_cones) do
		#####print(datas,Tvertex,uni_cones[i][3],LAMBDA,xxi);
		if USE_DUAL = false then 
			FF:=Function_M_uni_expanded(T,Tvertex,uni_cones[i][3],LAMBDA,Small,xxi,N-h);
		else 
			FF:=Function_M_uni_expanded(T,Tvertex,uni_cones[i][2],LAMBDA,Small,xxi,N-h);
		fi;
		FM:=FF[2]; 
		FMfB:=FM*fB;;
		pole:=FF[1]; ##print(polex,pole);
		
		if pole=0 then 
			FM:=subs(epsilon=0,FM);
			FMfB:=coeff(FMfB,x,-h-1);
		else 
			FMfB:=simplify(coeff(FMfB,x,-h-1));
		fi;
		FMfB:=coeff(convert(series(FMfB,epsilon=0,nops(Small)),polynom),epsilon,0);
		out:=out+uni_cones[i][1]*FMfB;
	od;
	#print("total time in uni_cones+the rest",time()-t1);

	simplify(-f*out*(-1)^h/h!);
end: 

# Input: AA list of positive integers ; h a nonnegative integer
# x a variable; 
# Output; a function of x,
# Math: the coefficient of the function H(AA,f)(T,x) in x^(-h-1), when f=1, that is we compute the coefficient of t^h of E(A,1)(t).  
# part_one_of_knapsack([1,4,9,10],3):=1/2160;
part_one_of_knapsack:=proc(AA,h) local fB,i,N,HH,HHH;

	fB:=1: N:=nops(AA)-1;
	for i from 1 to nops(AA) do ;
		fB:=fB*1/(AA[i]*x)*Todd(AA[i]*x,N-h);
	od; ##print(aqui);
	 HHH:=-coeff(fB,x,-h-1)*(-1)^h/h!;
end:

# The knapsack the straight way: to test some values
XXX:=proc(n);[seq(seq(seq([i,j,k],i=0..n),j=0..n),k=0..n)];end:
XXXX:=proc(n);[seq(seq(seq(seq([i,j,k,ell],i=0..n),j=0..n),k=0..n),ell=0..n)];end:

XXXXX:=proc(n);[seq(seq(seq(seq(seq([i,j,k,ell,p],i=0..n),j=0..n),k=0..n),ell=0..n),p=0..n)];end:
valueknapsackXXX:=proc(AAA,t) local test,u,R,s; test:=XXX(t);s:=0;
for u from 1 to nops(test) do R:=add(AAA[i]*test[u][i],i=1..3); 
if R=t then  s:=s+1;fi; od; s;end:
  
valueknapsackXXXX:=proc(AAAA,t) local test,u,R,s; test:=XXXX(t);s:=0;
for u from 1 to nops(test) do R:=add(AAAA[i]*test[u][i],i=1..4); 
if R=t then   s:=s+1;fi; od; s;end:

#printlevel:=0;valueknapsackXXX(A123,3);

valueknapsackXXXXX:=proc(AAAAA,t) local test,u,R,s; test:=XXXXX(t);s:=0;
for u from 1 to nops(test) do R:=add(AAAAA[i]*test[u][i],i=1..5); 
if R=t then   s:=s+1;fi; od; s;end: 
# Patch Function
# 
# INPUT: a list of positive integers L and  a positive integer i 
# OUTPUT: a list 
# Math: we compute the elements in L that are divided by i. If i is not a member of L then we return the empty list
# Example:is_divisor([2,4,3,7,12],5):=[]; is_divisor([2,4,3,7,12],4):=[4,12]
# 
is_divisor:=proc(L,i) local t,LL,newi,j,newL;
	LL:=sort(L);newL:=[];
	if member(i,LL,'t')=true then
		newi:=LL[t];
		for j from t to nops(LL) do
			if modp(LL[j],newi)=0 then newL:=[op(newL),LL[j]];
			else newL:=newL;
			fi:
		od:
	fi;
	newL;
end:

# INPUT: a list of  positive integers L and  a positive integer i 
# OUTPUT: a number
# Math: we compute the elements in L that divide  i. If i is not a member of L then we return the empty list
# Example:is_divided([2,4,3,7,12],12):=[2, 3, 4, 12];; 
is_divided:=proc(L,i) local t,LL,newL,newi,j;
	LL:=sort(L);newL:=[];
	if member(i,LL,'t')=true then 
		newi:=LL[t];
		for j from 1 to t do
			if modp(newi,LL[j])=0 then newL:=[op(newL),LL[j]];
			else newL:=newL;
			fi:
		od:
	fi;
	newL;
end:

# Input: a list of  positive integers L, i and j two elements of L with i<=j
# Output: an integer
# Math: we compute mu(i,j) following the patch formula 
# Example: mmu([2,4,3,7,12],3,12);=-1; mmu([2,4,3,7,12],3,7):=0;mmu([2,4,3,7,12],12,12):=1
mmu:=proc(L,i,j) local out,L1,L2,indexi,newL,indexj;#L is a list of common divisor
	if i=j then out:=1;
	elif modp(j,i)<>0 #that is L[i] doesn't divide L[j]
	then out:=0;
	else 
		L1:=is_divisor(L,i);
		newL:=is_divided(L1,j);
		newL:=subsop(1=NULL,newL);
		if member(j,newL,'s')=true then indexj:=s;fi;
		out:=-add(mmu(newL,newL[k],newL[indexj]),k=1..nops(newL));
	fi;
	out;
end:


# Input: a list of  positive integers L, i  an element of L
# Output: an integer
# Math: we compute ϱ(i) following the patch formula 
# Example: rho_i([2,3,4,6,8,7],2):= -1;
rho_i:=proc(L,i) local L1,k,rho;
	L1:=is_divisor(L,i):
	rho:=add(mmu(L,i,L1[k]),k=1..nops(L1));###print(seq(mmu1(L,i,k),k=1..nops(L1));
	rho;
end:

# Input: a list of  positive integers L
# Output:  a linear combination of G(i), i in L
# Math: we compute ϱ following the patch formula.
# Example :rho_knap([1,2,3]):=-G(1)+G(2)+G(3);
rho_knap:=proc(L) local out2,co,i,L1,rho_i,x ;
	out2:=0;co:=0;
	for i from 1 to nops(L) do 
		L1:=is_divisor(L,L[i]):
		rho_i:=add(mmu(L,L[i],L1[k]),k=1..nops(L1));
		out2:=out2+rho_i*G(L[i]);co:=co+rho_i;
	od;
	x:=co-1;
	out2:=out2-x*G(1);
end:

# The Complete knapsack and one coeff.
# Input: t a variable, A a list of nonnegative integers
# Math: We compute the knapsack for A and t: that is the whole Ehrhart polynomial
# Example:=complete_knapsack([1,2,3],17)=403/12+Smallstep(-17/2)^2-Smallstep(-17/2)+(3/2)*Smallstep(-17/3)^2-(3/2)*Smallstep(-17/3);
complete_knapsack:=proc(A,t) local h,out,F,p,i,f,c,N;
	out:=0;
	N:=nops(A)-1;
	F:=F_k(A,N):##print("F",F);
	p:=rho_knap(F): #print(p);
	for i from 1 to nops(F) do
		f:=F[i];
		c:=coeff(p,G(f));
		if f=1 then 
			out:=out+c*add(part_one_of_knapsack(A,h)*t^h,h=0..N);
		else 
		out:=out+c*add(part_f_of_knapsack(A,T,f,h)*t^h,h=0..N);fi ;
	od;
	eval(subs(T=t,simplify(out)));
end:

#printlevel:=0;complete_knapsack([1,2,3],t);
# Math: compute the Ehrhart coefficient of t^(N-k)
coeff_Nminusk_knapsack:=proc(A,t,k) local out,N,F,i,c,p,f,noError,Tseries;
	out:=0;
	N:=nops(A)-1;
	F:=F_k(A,k):##print(rootsk,F);
	#print("F=", F);
	p:=rho_knap(F):
	#print("p=", p);
	for i from 1 to nops(F) do
		f:=F[i];
		c:=coeff(p,G(f));
		#print("i=",i,"c=", c);
		if c = 0 then
			next;
		fi;
		
		noError := 0;
		Tseries:=0;
		while noError = 0 do
			try
				if f=1 then 
					Tseries:=part_one_of_knapsack(A,N-k)*t^(N-k);
				else 
					Tseries:=part_f_of_knapsack(A,T,f,N-k)*t^(N-k);
				fi;
				noError := 1;
			catch:
				fprintf("Error. trying again\n");
			end;
		od;
		
		out:=out+c*Tseries;
			
	od;
	#eval(subs(T=t,simplify(out)));
	eval(simplify(out));
end:

SMALLSTEP:=proc(T);
ceil(T)-T;
end:

numeric_knapsack:=proc(A,value);
	eval(subs(t=value,subs(Smallstep=SMALLSTEP,complete_knapsack(A,t))));
end:

# SOME EXPERIMENTS;
ASou:=proc(dimension) local n; n:=dimension; 
	[1,seq(n!/i,i=1..n)]; 
end:

random_list:=proc(bound,dimension) local R;
	R:=rand(bound);
	[1,seq(R()+1,i=1..dimension)]:
end:

testTopknapsackASou:=proc(dimension,k) local start1,finish; start1:=time();
	print(coeff_Nminusk_knapsack(ASou(dimension),t,k));finish:=time()-start1;
end:

testTopknapsackrandom:=proc(dimension,k) local start1,AA,n;
	n:=dimension;
	AA:=random_list(100,n); print(AA); 
	start1:=time();print(coeff_Nminusk_knapsack(AA,t,k));time()-start1;
end:

#testTopknapsackASou(18,2); #22 seconds
#testTopknapsackrandom(18,2);# 114 seconds;
#testTopknapsackrandom(11,3);

# First we have a few random computations (but without correctness verification) 
# Next we compare the calculations with 10 well-known knapsack problems for which
# we know everything (inside the file testedknapsacks) thanks to LattE. We provide
# some 
# 
#read("testedknapsacks");

#for i from 1 to 10 do
#informa[i]:=[]:
#od:
#for i from 6 to 6 do 
# st:=time();
# if numeric_knapsack(-L[i],m[i]) <> truenum[i] then
# print(`MISTAKE`);
# fi:
# informa[i]:=time()-st;
#od;




smallstep3 := proc (s) eval(ceil(s)-s) end proc:


#This function computes the top k+1 coefficients 
#and returns the numerical answer. So this 
#function is best called with a numeric t 
#(but a symbolic t would work too)

topk_numeric_polynomial:=proc(A,t,k) local poly, i, onlyOneTerm;
	poly:=0;
	for i from 0 to k do
	  onlyOneTerm:=coeff_Nminusk_knapsack(A, t, i);
	  poly:=poly + onlyOneTerm;
	  print("sum of first ", i, "=", evalf(eval(subs(Smallstep=smallstep3, poly))));
	  print("only the term ",i,evalf(eval(subs(Smallstep=smallstep3,onlyOneTerm))));

	end;
	return evalf(subs(Smallstep=smallstep3, poly));
end:


# 

#makes a list [1, 2, 3, ..., k]
makePartitionEquation:=proc(k)
	local i, L;
	L:=[];
	for i from 1 to k do
	  L:=[op(L), i];
	end; #for i
	return L;
end:


#print the function myceil(x) in latex ceil
	`latex/myceil`:=proc(x)
	return cat("\\ceil{",convert(x, string),"}");
end:

#print the function Smallstep as {{-x}} in latex
`latex/Smallstep`:=proc(x)
 	local y;
	y:=-1*x;
 	return cat("\\{", convert(y,string),"\\}");
end:





#make a table of the e.coefficients for the 
#partition knapsacks.

runPartitionTable:=proc()
local fptr, m, k, c;
local fileName;

for m from 5 to 5 do
  fileName:="ToughTable"||".m"||m||".tex";

  for k from ceil(m/2) to ceil(m/2) do
     printf("Starting partition %d, coefficient k=%d\n",m, k);
    # writeto("/dev/null"); 
     #don't print anything from coeff_Nminusk_knapsack
     c:=coeff_Nminusk_knapsack(makePartitionEquation(m),t,k);
     #c:=eval(subs(Smallstep=smallstep3, c));
    #c:=simplify(eval(subs(ceil=myceil,c)));
     print("c=",c);
     interface(screenwidth=9999);
     #print("latex");
     appendto(fileName);
     printf("& $");
     latex(c); printf("$");
     interface(screenwidth=79);
     writeto(terminal);
  od; #for k
od; #for m
end: writeto(terminal):



#will find the Ehrhart polynomial for some knapsacks and print the resutls to a file.

printPartitionPolynomial:=proc()
	local fptr, m, k, c;
	local fileName;
	fileName:="partitionTable3.mpl";
	fptr:=fopen(fileName, WRITE, TEXT);
	for m from 2 to 5 do
	  fprintf(fptr, "polynomial%d:= \\;", m);
	  for k from 0 to (m-1) do
	  #for k from 0 to (5-1) do
	     printf("Starting partition %d, coefficient k=%d\n",m, k);
	     writeto("/dev/null"); 
	     c:=coeff_Nminusk_knapsack(makePartitionEquation(m),t,k);
	     #c:=coeff_Nminusk_knapsack([1, 20, 7, 70, 400],t,k);
	     writeto(terminal);
	     fprintf(fptr, "\n+ %a \\", c);
	  od; #for k
	  fprintf(fptr, ";\n\n");
	od; #for m
	fclose(fptr);
end: 

writeto(terminal);

fractionalpart:=proc(s) local our,T;
    if  type(s,rational) then our:=s-floor(s);
    else our:={-s};
    fi;
    RETURN(our);
end:


