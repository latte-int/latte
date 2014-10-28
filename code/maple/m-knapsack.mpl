with(combinat):with(LinearAlgebra):with(linalg):with(numtheory):
# Outline
# Let A=[A_1,..,A_{N+1}] be a list of  positive integers and t a variable.
# We compute the number of  integral solutions for A_1x_1+..+A_{N+1}x_{N+1}=t.
# The output of this as a step polynomial function  of t.
#
# This program is aimed at computing efficiently the highest k coefficients
#  of this step polynomial.
# 
# 
# 
# Main Functions
# 1. complete_knapsack(A,t,T). Computes the Ehrhart polynomial for the knapsack A1x1+..A_{N+1}x_{N+1}=t, t a variable or an integer. The output is a quasi polynomial function of t and T, expressed as a function of T^k and MOD(a*t, b). MOD(a*t,b)  is the number in the half-open interval [0,b) that is equal to  a*t mod n.
# 
# 2. coeff_Nminusk_knapsack(A,t,T,k). Computes the N-k-degree term in the Ehrhart polynomial. The output is a periodic function of t and polynomial in T.
#
# 3. latteMod(a,b). Computes the number in the half-open interval [0,b) that is equal to  a mod b.
#
# Note that in the Ehrhart polynomial, the periodic coefficients are functions of t, while the monomial terms use T. This is done so that you can use Maple's coeff() function if you wish. T and t should always be evaluated at the same point.  
#
# Example usage
#  L:=[1,2,3,4,5];
#  coeff4minus0:=coeff_Nminusk_knapsack(L, t, T, 0); #find the T^4 term
#  coeff4minus3:=coeff_Nminusk_knapsack(L, t, T, 3); #find the T^1 term
#  eval(subs({T=10, t=10, MOD=latteMod}, coeff4minus3)); #evaluate the t^1 term at t=T=10.
#  
#  ehrhartPoly:=complete_knapsack(L,t,T);
#  coeff(ehrhartPoly, T^4); #get the leading term's coefficient
#  eval(subs({T=10, t=10, MOD=latteMod}, ehrhartPoly)); #number of solutions to x_1 + 2x_2 + 3x+3 + 4x_4 +5x_5 = 10, x_i >= 0...the answer is 30
#
# License
#   This Maple script is part of the LattE integrale 1.7 package made available
#   under the GNU General Public License at http://www.math.ucdavis.edu/~latte/.
#	
# Testing
#   Define this variable to have the code check self-tests (some are expensive)
#  (this is for "make check", not for production code).
#CHECK_EXAMPLES := true:

#   Enable checking ASSERTions
kernelopts(assertlevel=1):       

##################################
#Start example checking functinos#
##################################

# Like ASSERT, but with better reporting (for automatic tests).
# Example usage: TEST_EQUAL("sqrt(4)", "2", "sqrt test #1");
TEST_EQUAL:=proc(a_expr, b_expr, message) local a, b;
    printf("Checking %s\n", message);
    a := eval(parse(a_expr));
    b := eval(parse(b_expr));
    if not (a = b) then
        printf("FAIL: Results differ:\n  Got:       %a\n  Should be: %a\n", a, b);
    fi:
end:

# Find out whether we are to check examples.
# Uses the global CHECK_EXAMPLES
# Example usage 
#   if check_examples() then
#     TEST_EQUAL("sqrt(4)", "2", "sqrt test #1");
#     TEST_EQUAL("sqrt(0)", "0", "sqrt test #2");
#     TEST_EQUAL("2^3", "8", "power test #1");
#   fi;
check_examples:=proc()
    type(CHECK_EXAMPLES, boolean) and CHECK_EXAMPLES;
end:



############################
#Start of library functions#
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
if check_examples() then
	TEST_EQUAL("primitive_vector([2,2,2,8,16,100])", "[1, 1, 1, 4, 8, 50]", "test primitive_vector #1");
	TEST_EQUAL("primitive_vector([1234/1345, 1234/122, 2457/3568, 4, -24/2457, -6])", "[219965080608, 2425024864080, 165097758735, 959004970560, -2341892480, -1438507455840]", "test primitive_vector #2");
fi;


#COMPLEMENT LIST 
#The output is the Complement  List, within the list [1,..,d]
ComplementList:=proc(K,d);
	RETURN([seq (`if` (member(i,K)=false, i, op({})),i=1..d)]);
end:

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
if check_examples() then
	TEST_EQUAL("short_vector([[2200,2200],[2,100]])", "[2, 100]", "test short_vector #1");
	TEST_EQUAL("short_vector([[-2200,-2200,-2346],[2,100,-565], [457,568,-2457]])", "[449, 168, -197]", "test short_vector #2");		
fi;

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
if check_examples() then
	TEST_EQUAL("good_vector([[2200,2129],[101,100]])", "[[31, 30], [[1, 2], [], []]]", "test good_vector #1");	
fi;


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
if check_examples() then
	TEST_EQUAL("good_cone_dec(-1, [[2200,2129],[101,100]])", "[[[1, 70, [[-101, -100], [31, 30]]]], [[-1, 1, [[2200, 2129], [31, 30]]]]]", "test good_cone_dec #1");
	TEST_EQUAL("good_cone_dec(-1, [[2200,2129,5675],[101,100,-4545], [45613,345,-135]])", " [[[-1, -10243805, [[1, 0, 0], [2200, 2129, 5675], [101, 100, -4545]]], [-1, -1554525, [[45613, 345, -135], [-1, 0, 0], [101, 100, -4545]]], [1, -2245290, [[45613, 345, -135], [-1, 0, 0], [-2200, -2129, -5675]]]], []]", "test good_cone_dec #2");	
fi;


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
if check_examples() then
	TEST_EQUAL("cone_dec([[2200,2129],[101,100]])", "[[1, 1, [[2200, 2129], [31, 30]]], [-1, 1, [[-1, -1], [31, 30]]], [1, 1, [[-1, -1], [101, 100]]]]", "test cone_dec #1");
	TEST_EQUAL("cone_dec([[7,12,3],[101,0,0], [7,5,-5]])", "[[-1, 1, [[3, 4, 0], [-7, -12, -3], [7, 5, -5]]], [1, 1, [[3, 4, 0], [-7, -12, -3], [0, 1, 1]]], [1, 1, [[0, 1, 1], [1, 1, 0], [-1, 0, 0]]], [1, 1, [[-2, -2, 1], [7, 5, -5], [0, 1, 1]]], [1, 1, [[-1, -1, 1], [2, 2, -1], [0, 1, 1]]], [-1, 1, [[0, -1, -1], [0, 1, 0], [1, 0, 0]]], [1, 1, [[1, 1, -1], [0, 1, 0], [1, 0, 0]]], [1, 1, [[1, 1, -1], [0, 1, 1], [0, -1, 0]]], [1, 1, [[5, 4, -3], [0, -1, -1], [7, 5, -5]]], [-1, 1, [[5, 4, -3], [-3, -4, 0], [7, 5, -5]]], [1, 1, [[5, 4, -3], [-3, -4, 0], [0, 1, 1]]], [1, -1, [[1, 0, 0], [2, 4, 1], [1, 1, 0]]], [-1, -1, [[7, 12, 3], [2, 4, 1], [1, 1, 0]]], [-1, 1, [[0, 1, 1], [1, 1, 0], [3, 5, 1]]], [1, 1, [[7, 12, 3], [1, 1, 0], [3, 5, 1]]], [-1, 1, [[7, 12, 3], [0, -1, -1], [3, 5, 1]]], [-1, -1, [[-4, -3, 3], [-7, -5, 5], [2, 2, -1]]], [1, -1, [[-4, -3, 3], [1, 1, -1], [2, 2, -1]]]]", "test cone_dec #2");	
fi;

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
if check_examples() then
	TEST_EQUAL("new_short_vector([[2200,2200],[2,100]])", "[2, 100]", "test new_short_vector #1");
	TEST_EQUAL("new_short_vector([[-2200,-2200,-2346],[2,100,-565], [457,568,-2457]])", "[449, 168, -197]", "test new_short_vector #2");		
fi;


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
	ok:=0;i:=1; 
	while i<= nops(VV) and ok=0 do  
		if VV[i]<=0 then i:=i+1;
		else ok:=1; i:=nops(VV)+1;
		fi; 
	od:   
	if ok=1 then 
		out:=V; 
	else 
		out:=[seq(-V[i],i=1..nops(V))];
	fi;
	out;
end:
if check_examples() then
	TEST_EQUAL("polar_good_vector([[2200,2129],[101,100]])", "[31, 30]", "test polar_good_vector #1");		
fi;

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
       		else 
       			Uni:=[op(Uni),Csigned];
        	fi;
        else fi;
   	od;
 	[Nonuni,Uni];
 end:
if check_examples() then
	TEST_EQUAL("polar_signed_decomp(1, [[101,2],[7,5]], [7,6])", "[[[-1, [[7, 6], [7, 5]]], [1, [[101, 2], [7, 6]]]], []]", "test polar_signed_decomp #1");		
fi; 
 
 
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
if check_examples() then
	TEST_EQUAL("polar_good_cone_dec(-1, [[2200,2129],[101,100]])", "[[[-1, [[31, 30], [101, 100]]]], [[-1, [[2200, 2129], [31, 30]]]]]", "test polar_good_cone_dec #1");	
fi;



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
if check_examples() then
	TEST_EQUAL("polar_cone_dec([[2200,2129],[101,100]])", "[[1, [[2200, 2129], [31, 30]]], [-1, [[1, 1], [101, 100]]], [1, [[31, 30], [1, 1]]]]", "test polar_cone_dec #1");
	TEST_EQUAL("polar_cone_dec([[7,12,3],[11,0,0], [7,5,-5]])", "[[1, [[7, 12, 3], [3, 4, 0], [7, 5, -5]]], [-1, [[1, 1, 0], [1, 0, 0], [-1, -1, 1]]], [-1, [[2, 1, -2], [1, 1, 0], [7, 5, -5]]], [-1, [[3, 4, 0], [2, 1, -2], [7, 5, -5]]], [1, [[3, 5, 1], [1, 1, 0], [3, 4, 0]]], [-1, [[7, 12, 3], [3, 5, 1], [3, 4, 0]]], [1, [[7, 12, 3], [1, 1, 0], [3, 5, 1]]], [1, [[2, 4, 1], [1, 0, 0], [1, 1, 0]]], [-1, [[7, 12, 3], [2, 4, 1], [1, 1, 0]]], [-1, [[1, 1, 0], [4, 3, -3], [7, 5, -5]]], [1, [[1, 1, 0], [-1, -1, 1], [4, 3, -3]]], [1, [[3, 3, -1], [1, 1, 0], [2, 1, -2]]], [1, [[3, 4, 0], [3, 3, -1], [2, 1, -2]]], [1, [[3, 4, 0], [1, 1, 0], [3, 3, -1]]]]", "test polar_cone_dec #2");
fi;



# Cone decomposition via polar cone
primitive_vector:=proc(A) local d,n,g;
	d:=nops(A);
	n:=ilcm(seq(denom(A[i]),i=1..d));
	g:=igcd(seq(n*A[i],i=1..d));if g<>0 then
	[seq(n*A[i]/g,i=1..d)];else [seq(n*A[i],i=1..d)];fi;
end:
if check_examples() then
	TEST_EQUAL("primitive_vector([2200,2120,2340,240])", "[110, 106, 117, 12]", "test primitive_vector #1");
fi;


#The polar cone.
polarcone:=proc(A) local n,B,Ainverse,out;
	n:=nops(A);
	#MatrixInverse(Transpose(Matrix(A)));
	B:=Matrix(A); #print("B",Matrix(A),B);     
   Ainverse:=MatrixInverse(B);#print(Ainverse);
	out:=[seq(convert(Transpose(Column(Ainverse,i)),list),i=1..n)];
end:
if check_examples() then
	TEST_EQUAL("polarcone([[2200,2129],[101,100]])", " [[100/4971, -101/4971], [-2129/4971, 2200/4971]]", "test polarcone #1");
	TEST_EQUAL("polarcone([[1,0],[0,-1]])", "[[1,0],[0,-1]]", "test polarcone #2");
fi;


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
if check_examples() then
	TEST_EQUAL("new_dual_dec_cone([[1,0],[1,10]])", "[[1, [[1, 0], [0, 1]]], [1, [[0, -1], [1, 10]]]]", "test new_dual_dec_cone #1");
fi;


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
end:
if check_examples() then
	TEST_EQUAL("cone_in_basis_Lambda([[1,0],[1,10]], [[1,0],[1,-10]])", "[[1, 0], [2, -1]]", "test cone_in_basis_Lambda #1");
fi;

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
if check_examples() then
	TEST_EQUAL("Basis_Lattice([23,43,21,43,3,1,3], 5234)", "[[0, 0, 0, 0, 0, 5234, 0], [1, 0, 0, 0, 0, -23, 0], [0, 1, 0, 0, 0, -43, 0], [0, 0, 1, 0, 0, -21, 0], [0, 0, 0, 1, 0, -43, 0], [0, 0, 0, 0, 1, -3, 0], [0, 0, 0, 0, 0, -3, 1]]", "test Basis_Lattice #1");
fi;

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
if check_examples() then
	TEST_EQUAL("B_and_S([2,2,2,3,4,5],2)", "[[2, 2, 2, 4], [3, 5]]", "test B_and_S #1");
fi;

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
if check_examples() then
	TEST_EQUAL("myChoose([2,2,2,3,4,5],2)", "[[2, 2], [2, 2], [2, 3], [2, 4], [2, 5], [2, 2], [2, 3], [2, 4], [2, 5], [2, 3], [2, 4], [2, 5], [3, 4], [3, 5], [4, 5]]", "test myChoose #1");
fi;




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
if check_examples() then
	TEST_EQUAL("F_k([2,3,3,55],3)", "[1, 2, 3, 55]", "test F_k #1");
	TEST_EQUAL("F_k([2,3,3,55],2)", "[1, 3]", "test F_k #2");	
fi;


# The coefficient of t^h in the function E(A,f)(t):
fract:=proc(a,n) local out,p,q,newa,t; 
	if  a=0 then out:=0; 
	elif type(simplify(a),integer)=true  then out:=0;
	else 
		p:=numer(a); q:=denom(a); 
		newa:=modp(p,q);
		t:=MOD(newa/q*n,1);
	fi;
end:
if check_examples() then
	TEST_EQUAL("fract(2/3,-1)", "MOD(-2/3,1)", "test fract #1");
	TEST_EQUAL("fract(3/4,2)", "MOD(3/2,1)", "test fract #2");	
fi;


smallstep:=proc(a,n); fract(-a,n); end:

Todd:=proc(x,order):
-add(1/j!*bernoulli(j,0)*x^j,j=0..order):
end:
if check_examples() then
	TEST_EQUAL("Todd(t,7)", "-1+1/2*t-1/12*t^2+1/720*t^4-1/30240*t^6", "test Todd #1");
fi;


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
if check_examples() then
	TEST_EQUAL("is_divisor([2,4,3,7,-12],5)", "[]", "test is_divisor #1");
	TEST_EQUAL("is_divisor([2,4,3,7,12],4)", "[4, 12]", "test is_divisor #2");	
fi;


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
if check_examples() then
	TEST_EQUAL("is_divided([2,4,3,7,12],12)", "[2, 3, 4, 12]", "test is_divided #1");
	TEST_EQUAL("is_divided([2,4,3,7,-12],17)", "[]", "test is_divided #2");	
fi;


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
if check_examples() then
	TEST_EQUAL("mmu([2,4,3,7,12],2,12)", "0", "test mmu #1");
	TEST_EQUAL("mmu([2,4,3,7,12],12,12)", "1", "test mmu #2");
	TEST_EQUAL("mmu([2,4,3,7,12],4,7)", "0", "test mmu #3");	
fi;



# Input: a list of  positive integers L, i  an element of L
# Output: an integer
# Math: we compute rho(i) following the patch formula 
# Example: rho_i([2,3,4,6,8,7],2):= -1;
rho_i:=proc(L,i) local L1,k,rho;
	L1:=is_divisor(L,i):
	rho:=add(mmu(L,i,L1[k]),k=1..nops(L1));###print(seq(mmu1(L,i,k),k=1..nops(L1));
	rho;
end:
if check_examples() then
	TEST_EQUAL("rho_i([2,4,3,7,12],2)", "0", "test rho_i #1");
	TEST_EQUAL("rho_i([2,4,3,7,12],12)", "1", "test rho_i #2");
	TEST_EQUAL("rho_i([2,4,3,7,12],2)", "0", "test rho_i #3");	
fi;


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
if check_examples() then
	TEST_EQUAL("rho_knap([2,4,3,7,12])", "G(7)+G(12)-G(1)", "test rho_knap #1");
	TEST_EQUAL("rho_knap([345,234,6245,3457,4526,3126,2457,235,1,1,2,3,45,456,12,43,67,2,46,2])", "G(345)+G(234)+G(6245)+G(3457)+G(4526)+G(3126)+G(2457)+G(235)-G(1)-6*G(2)-5*G(3)+G(45)+G(456)+G(43)+G(67)+G(46)", "test rho_knap #2");	
fi;

#####################################
#Next 3 are main interface functions#
#####################################

# The Complete knapsack and one coeff.
# Input: t a variable, A a list of nonnegative integers
# Math: We compute the knapsack for A and t: that is the whole Ehrhart polynomial
# Example:=complete_knapsack([1,2,3],17)=403/12+Smallstep(-17/2)^2-Smallstep(-17/2)+(3/2)*Smallstep(-17/3)^2-(3/2)*Smallstep(-17/3);
complete_knapsack:=proc(A,t,T) local h,out,F,p,i,f,c,N;
	out:=0;
	N:=nops(A)-1;
	F:=F_k(A,N):##print("F",F);
	p:=rho_knap(F): #print(p);
	for i from 1 to nops(F) do
		f:=F[i];
		c:=coeff(p,G(f));
		if f=1 then 
			out:=out+c*add(part_one_of_knapsack(A,h)*T^h,h=0..N);
		else 
		out:=out+c*add(part_f_of_knapsack(A,t,f,h)*T^h,h=0..N);fi ;
	od;
	#eval(subs(T=t,simplify(out)));
	eval(simplify(out));
end:
if check_examples() then
	TEST_EQUAL("complete_knapsack([2,4,3],t,T)", "1-3/2*MOD(1/3*t,1)+3/2*MOD(1/3*t,1)^2-3/2*MOD(1/2*t,1)+1/32*(4*MOD(3/4*t,1)+4*MOD(1/2*t,1))^2-3/2*MOD(3/4*t,1)^2+(1/4-1/4*MOD(1/2*t,1))*T+1/48*T^2", "test complete_knapsack #1");
	TEST_EQUAL("complete_knapsack([1,6,6,6,6],t,T)", "1-25/12*MOD(1/6*t,1)+35/24*MOD(1/6*t,1)^2-5/12*MOD(1/6*t,1)^3+1/24*MOD(1/6*t,1)^4+(25/72-35/72*MOD(1/6*t,1)+5/24*MOD(1/6*t,1)^2-1/36*MOD(1/6*t,1)^3)*T+(35/864-5/144*MOD(1/6*t,1)+1/144*MOD(1/6*t,1)^2)*T^2+(5/2592-1/1296*MOD(1/6*t,1))*T^3+1/31104*T^4", "test complete_knapsack #2");	
fi;


#Compute the Ehrhart coefficient of t^(N-k)
coeff_Nminusk_knapsack:=proc(A,t,T, k) local out,N,F,i,c,p,f,noError,Tseries;
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
					Tseries:=part_one_of_knapsack(A,N-k)*T^(N-k);
				else 
					Tseries:=part_f_of_knapsack(A,t,f,N-k)*T^(N-k);
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
if check_examples() then
	TEST_EQUAL("coeff_Nminusk_knapsack([2,4,3],t,T,0)", "1/48*T^2", "test coeff_Nminusk_knapsack #1");
	TEST_EQUAL("coeff_Nminusk_knapsack([2,4,3],t,T,1)", "-1/4*(-1+MOD(1/2*t,1))*T", "test coeff_Nminusk_knapsack #2");
	TEST_EQUAL("coeff_Nminusk_knapsack([2,4,3],t,T,2)", "1-3/2*MOD(1/3*t,1)+3/2*MOD(1/3*t,1)^2-3/2*MOD(1/2*t,1)+1/32*(4*MOD(3/4*t,1)+4*MOD(1/2*t,1))^2-3/2*MOD(3/4*t,1)^2", "test coeff_Nminusk_knapsack #3");
	TEST_EQUAL("coeff_Nminusk_knapsack([1,6,6,6,6],t,T,0)", "1/31104*T^4", "test coeff_Nminusk_knapsack #4");
	TEST_EQUAL("coeff_Nminusk_knapsack([1,6,6,6,6],t,T,1)", "-1/2592*(-5+2*MOD(1/6*t,1))*T^3", "test coeff_Nminusk_knapsack #5");
	TEST_EQUAL("coeff_Nminusk_knapsack([1,6,6,6,6],t,T,2)", "1/864*(35-30*MOD(1/6*t,1)+6*MOD(1/6*t,1)^2)*T^2", "test coeff_Nminusk_knapsack #6");
	TEST_EQUAL("coeff_Nminusk_knapsack([1,6,6,6,6],t,T,3)", " -1/72*(2*MOD(1/6*t,1)^3-15*MOD(1/6*t,1)^2+35*MOD(1/6*t,1)-25)*T", "test coeff_Nminusk_knapsack #7");
	TEST_EQUAL("coeff_Nminusk_knapsack([1,6,6,6,6],t,T,4)", "1-25/12*MOD(1/6*t,1)+35/24*MOD(1/6*t,1)^2-5/12*MOD(1/6*t,1)^3+1/24*MOD(1/6*t,1)^4", "test coeff_Nminusk_knapsack #8");
fi;


#Computes the top k+1 terms: T^N, T^{N-1}, ..., T_{N-k}
knapsackKTerms:=proc(A, t, T, k)
	local ans, i;
	ans:=0;
	for i from 0 to k do
		ans:=ans + coeff_Nminusk_knapsack(A, t, T, i);
	end;
	
	return ans;
end:
if check_examples() then
	TEST_EQUAL("knapsackKTerms([2,4,3],t,T,1)", "1/48*T^2-1/4*(-1+MOD(1/2*t,1))*T", "test knapsackKTerms #1");
	TEST_EQUAL("knapsackKTerms([2,4,3],t,T,2)", "1/48*T^2-1/4*(-1+MOD(1/2*t,1))*T+1-3/2*MOD(1/3*t,1)+3/2*MOD(1/3*t,1)^2-3/2*MOD(1/2*t,1)+1/32*(4*MOD(3/4*t,1)+4*MOD(1/2*t,1))^2-3/2*MOD(3/4*t,1)^2", "test knapsackKTerms #2");
fi;



SMALLSTEP:=proc(T);
ceil(T)-T;
end:
if check_examples() then
	TEST_EQUAL("SMALLSTEP(2/3)", "1/3", "test SMALLSTEP #1");
	TEST_EQUAL("SMALLSTEP(3+2/3)", "1/3", "test SMALLSTEP #2");
	TEST_EQUAL("SMALLSTEP(-3-1/3)", "1/3", "test SMALLSTEP #3");
	TEST_EQUAL("SMALLSTEP(-4/5)", "4/5", "test SMALLSTEP #4");
	TEST_EQUAL("SMALLSTEP(5/4)", "3/4", "test SMALLSTEP #5");
	TEST_EQUAL("SMALLSTEP(5)", "0", "test SMALLSTEP #6");
fi;


numeric_knapsack:=proc(A,value);
	eval(subs(T=value,t=value,subs(Smallstep=SMALLSTEP,complete_knapsack(A,t,T))));
end:
if check_examples() then
	TEST_EQUAL("numeric_knapsack([1,1,1],30)", "496", "test numeric_knapsack #1");
fi;

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
if check_examples() then
	TEST_EQUAL("fractionalpart(2/3)", "2/3", "test fractionalpart #1");
	TEST_EQUAL("fractionalpart(3+2/3)", "2/3", "test fractionalpart #2");
	TEST_EQUAL("fractionalpart(-3-1/3)", "2/3", "test fractionalpart #3");
	TEST_EQUAL("fractionalpart(-4/5)", "1/5", "test fractionalpart #4");
	TEST_EQUAL("fractionalpart(5/4)", "1/4", "test fractionalpart #5");
	TEST_EQUAL("fractionalpart(5)", "0", "test fractionalpart #6");
fi;

#### LATTE INTERFACE FUNCTION:
#
# Function to be substituted for the formal MOD function to get
# results.
# (can as well just use value().)
#
# input:
#	a: any rational number or symbolic expression.
#	n: any rational number.
#	return the number in the half-open interval [0,n) that is equal to  a mod n
latteMod:=proc(x, n)
	ASSERT(n > 0);
    ### Note `floor' works fine for symbolics, if `Digits'
    ### is large enough. --mkoeppe
    x - floor(x/n)*n;
end:
if check_examples() then
	TEST_EQUAL("latteMod(5+1/3, 2)", "4/3", "test latteMod #1");
	TEST_EQUAL("latteMod(-1/3, 1)", "2/3", "test latteMod #2");
	TEST_EQUAL("latteMod(sqrt(2), 1)", "sqrt(2) - 1", "test latteMod #3"); # Note the difference to what we output in fractionalpart!
    
fi;

