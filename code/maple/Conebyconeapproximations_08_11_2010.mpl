with(linalg):with(LinearAlgebra):with(combinat):
kernelopts(assertlevel=1):       ### Enable checking ASSERTions

## Define this variable to have the code check self-tests (some expensive)
## (this is for "make check", not for production code; consider errorbreak (-e 2)).
#
# CHECK_EXAMPLES := true:

# PEDAGOCICAL PROGRAM FOR COMPUTING EXAMPLES FOR ARTICLE:
# HIGHEST EHRHART  COEFFICIENTS;
# Version of November 10 -2010 : I have added the Ehrhart polynomial over the reals.
#  .... with changes for LattE.
#
# I PUT SOME EXAMPLES AFTER THE PROCEDURE.
#
#
#
# HERE IS THE LIST OF WHAT THE PROGRAM DOES.
#
#
#
# PROGRAM FOR  FORMULAE A and b for S^{L^I}(s+c) .
# This function  is expressed in terms of the "black box functions"
# TODD, EXP, and either CEIL or MOD.
#
# TODD(s,x) representing the function e^(sx) x/(1-e^x). 

# IF we want to evaluate S_C_k(x) at a regular element reg; the command eval and
# subs (TODD=Todd, EXP=exp, CEIL=ceil, MOD=latteMod) should be used.
#
# The commands are EITHER:
# S_Ispace_Coneformulaa:=proc(vertex,cone, ISpace,xi);
#  (which gives expressions using the CEIL function)
# OR:
# S_Ispace_Coneformulab:=proc(vertex,cone, ISpace,xi);
#  (which gives expressions using the fractional part function MOD).
#
# Here  vertex should be entered  as a symbolic variable [s_1,s_2,s_3,....,s_d] ;
#  cone is our cone  with generators v_i; I space is given by  a subset of [1,2,...,d];
#
# Thus if  I is [];  we get S(s+c);
#  if S is [1,2,...,d] we get I(s+c) ;
# There are some examples  in another sheet:
#S_Ispace_Coneformulaa([s[1],s[2]],[[1,0],[1,2]],[1],xi);
#S_Ispace_Coneformulab([s[1],s[2]],[[1,0],[1,2]],[1],xi);
#
#
# It computes the coefficient in n^m of the weighted Ehrhart quasi polynomial  of a rational simplex;  Eh(Simplex,ell^M)(n)=sum _m E_m(n)*n^m;
# it return  E_m(n) a  periodic function of n; However, it might return an error message: say division by zero;
# In this case we have to use the next procedure.
#
#TopEhrhartweightedluckyell(n,Simplex,ell,M,m);

#
# EXAMPLE
#Delta:=[[0,0],[5/28,0],[5/28,5/14]];
#TopEhrhartweightedluckyell(n,Delta,2,[1,1],1);
#
#
# If luckyell  return an error message: that is ell is singular for one of the thousands cone occuring in the procedure, uses
# the next procedure. it might also return an error message:;
# In this case we have to rerun (there is a random choice there).
#
#TopEhrhartweighted:=proc(n,Simplex,ell,M,m)
#
# EXAMPLE
#Delta:=[[0,0],[5/28,0],[5/28,5/14]];
#TopEhrhartweighted(n,Delta,2,[1,0],1);
#
#
#
#

## FIXME:  Above examples regarding TopEhrhart are outdated; the
# argument list has changed.  See
# Conebyconeapproximations_08_11_2010_examples.mpl instead.

# Find out whether we are to check examples.
check_examples:=proc():
    type(CHECK_EXAMPLES, boolean) and CHECK_EXAMPLES;
end:

#
# Programs on lists: addition on lists, complement of a list, sublist,etc...
#
#The output is the Complement  List, within the list [1,..,d]
ComplementList:=proc(K,d);
    RETURN([seq (`if` (member(i,K)=false, i, op({})),i=1..d)]);
end:
special_lincomb_v:=proc(a,v,n) local out;
    ASSERT(nops(a)=nops(v)," the number of coefficients and vectors do not match");
    if v=[]   then out:=[seq(0,i=1..n)];else
        out:=[seq(add(a[i]*v[i][j],i=1..nops(v)),j=1..nops(v[1]))];
    fi;
    out;
end:

# Miscellanea
#
# Input: A :a vector with rational coordinates.
# Output: A vector with integral coordinates:
# Math: the primitive vector on the half line R^+A;
# Example: #

#
primitive_vector:=proc(A) local d,n,g;
    d:=nops(A);
    n:=ilcm(seq(denom(A[i]),i=1..d));
    g:=igcd(seq(n*A[i],i=1..d));
    if g<>0 then
        [seq(n*A[i]/g,i=1..d)];
    else
        [seq(n*A[i],i=1..d)];
    fi;
end:
if check_examples() then
    ASSERT(primitive_vector([0,-1/2]) = [0,-1], "primitive_vector test #1");
fi;

#
ortho_basis:=proc(d) local i,v;
    for i from 1 to d do
        v[i]:=[seq(0,j=1..i-1),1,seq(0,j=i+1..d)]
    od;
    [seq(v[j],j=1..d)];
end:

#####################################
# Computing Barvinok Signed decomposition into unimodular cones
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
#    printf("good_cone_dec: %A, %A", eps, G);
    n:=nops(G);  A:=matrix([seq(G[i],i=1..n)]);
    det_A:=det(A);
    if abs(det_A)=1 then
        Output:=[[],[[eps,det_A,G]]];
    else R:=good_vector(G);
        Output:=signed_decomp(eps,G,R[1],R[2]);
    fi;
    # print(Output);
    # Output;
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

#
# Projections:
# Input: W is a list of vectors  of V , [v[1],..v[d]], of lenght d.
# Cspace =[i[1]..,i[s]],  a list of integers.
#  b is a vector of lenght d.
# Output: a vector of lenght d.
#
# Math:
# We decompose the space V in lin(Cspace)+lin(ISpace) where lin(Cspace) is the linear span  of the vectors v[i], i in Cspace,
# and lin(ISpace) of the vectors in the complement indices. We project a vector b on lin(Cspace)
# Thus we write b=b_Cspace+b_ISpace;
# Our output is b_Cspace;
projectedvector:=proc(W,Cspace,b) local M,S,j,v,V,m;
    M:=transpose(matrix([seq(W[i],i=1..nops(W))]));
    S:=linsolve(M,b);
    m:=det(M);
    for j from 1 to nops(W) do
        v[j]:=add(S[Cspace[i]]*W[Cspace[i]][j],i=1..nops(Cspace));
    od:
    V:=[seq(v[j],j=1..nops(W))];
end:
if check_examples() then
    ASSERT(projectedvector([[1,0,0],[0,1,2],[0,1,0]],[3],[0,0,1]) = [0,-1/2,0], "projectedvector test #1");
end;

## Same, with precomputed inverse (faster)
projectedvector_with_inverse:=proc(M_inverse, W,Cspace,b) local S,j,v,V;
    S:=multiply(M_inverse,b);
    #m:=det(M);
    for j from 1 to nops(W) do
        v[j]:=add(S[Cspace[i]]*W[Cspace[i]][j],i=1..nops(Cspace));
    od:
    V:=[seq(v[j],j=1..nops(W))];
end:

# Projected lattice
# # Input:  W=[v1,v2,.., vd];  a "Cone"  in  R^d;
# BE CAREFUl: The vectors in W must have integral coordinates.
#  Cspace a subset of [1,2..d] of cardinal k;

# # Output a list [H1,H2,...,Hk] of  k vectors in R^d .
#
#
# projectedlattice:
# Math: we
# decompose V in lin(Cspace)+lin(ISpace);
#  we project the standard lattice (that is Ze[1]+..+Ze[d], that is  Z[1,0,0..0]+... Z.[0,0,0..,1]])
# on lin[Cspace] which is a  subspace of dimension k  of a space of dim d.
# output: (using ihermite) a basis of k elements (of lenght d) of the projected lattice  on lin(Cspace).
# We will use over and over again this list H1,H2,..., Hk, so that we will work in Z^k  (embedded in R^d via H1,H2,..Hk).
#
projectedlattice:=proc(W,Cspace) local m,B, d,k,i,r,S,IS,List,M_inverse, temp_projectedVectors;
    d:=nops(W);
    B:=ortho_basis(d);
    k:=nops(Cspace);
    m:=abs(det(transpose(matrix([seq(W[i],i=1..nops(W))]))));

    M_inverse:=inverse(transpose(matrix([seq(W[i],i=1..nops(W))])));
    for i from 1 to d do

        temp_projectedVectors:=m*projectedvector_with_inverse(M_inverse, W,Cspace,B[i]);
        r[i]:=[seq(temp_projectedVectors[j],j=1..nops(W))];
    od;
    S:=matrix([seq(r[i],i=1..d)]);;
    IS:=ihermite(S);
    List:=[seq(1/m*convert(row(IS,j),list),j=1..k)];
    List;
end:
if check_examples() then
    ASSERT(projectedlattice([[1,3,0],[0,1,0],[0,0,2]],[1,3]) = [[1, 3, 0], [0, 0, 1]], "projectedlattice test #1");
fi;

# Projected cone and projected vertex (expressed in the lattice basis)
# Input: W is a Cone in Z^d and Cspace is a subset of [1,..,d] of
#        cardinality k; ProjLattice = projectedlattice(W,Cspace);
# Output: A "Cone" in Z^k;

# Be careful: our input must have integral coordinates.
# The ouput then will have integral coordinates.
#
#
# Here W is the cone and we are projecting W over lin(Cspace) and
# expressing it in term of the basis H_1,H_2,...,H_k of
# projectedlattice(W,Cspace).
projectedconeinbasislattice:=proc(W,Cspace,ProjLattice) local P,M,output,i,F;
    P:=ProjLattice;
    M:=transpose(matrix([seq(P[i],i=1..nops(P))]));
    output:=[];
    for i from 1 to nops(Cspace) do
        F:=convert(linsolve(M,Vector(W[Cspace[i]])),list);
        output:=[op(output),primitive_vector(F)];
    od;
    output;
end:
if check_examples() then
    ASSERT(projectedconeinbasislattice([[1,1,0],[0,1,0],[0,0,2]],[1,3], projectedlattice([[1,1,0],[0,1,0],[0,0,2]], [1,3])) = [[1,0],[0,1]]);
fi;


# #Input: W a Cone in Z^d;
#         Cspace a subset of [1,2,..d] of cardinal k;
#         ProjLattice := projectedlattice(W,Cspace);
#         s a vector in R^d with rational coordinates (or symbolic coordinates);
# #Ouput: a vector in R^k with rational coordinates (or symbolic coordinates);

# Math: Here W is the cone and we are projecting V over lin( Cspace)  using  V:=lin(Cspace) oplus
#  lin(ISpace). We express the projection of s
# with respect to the basis of the projected lattice. If the ouput is [a1,a2], this means that our
# projected vertex is s_Cspace=a1*H1+a2*H2 where H1,H2 is the basis of the projected lattice computed before.
#
projectedvertexinbasislattice:=proc(W,Cspace,ProjLattice,s) local m,P,M,output,i,F;
    P:=ProjLattice;
    if Cspace=[] then RETURN([]);fi;
    M:=Transpose(Matrix([seq(P[i],i=1..nops(P))]));
    F:=convert(LinearSolve(M,Vector(projectedvector(W,Cspace,s))),list);
    output:=F;
end:
if check_examples() then
    ASSERT(projectedvertexinbasislattice([[1,0,0],[0,2,1],[0,1,1]],[1,3],projectedlattice([[1,0,0],[0,2,1],[0,1,1]], [1,3]), [s1,s2,s3]) = [s1, 2*s3-s2], "projectedvertexinbasislattice test #1");
fi;


# Input: s a vector in R^d with rational coordinates (or symbolic).
# W a cone in Z^d;
# Ispace a subset of [1,2,...,d];
# Output:  a vector in R^d
#
# Math: We decompose V in lin(CSpace) oplus lin (ISpace), (with CSpace spanned by the v|i] in the complementary indices of Ispace );
#  and here we write s=s_Cspace+s_(ISpace): Here the output is
# s_(ISpace);
s_ISpace:=proc(s,W,ISpace) local M,s_in_cone_coord,s_ISpace;
    M:=Matrix([seq(Vector([W[i]]),i=1..nops(W))]);
    s_in_cone_coord:=convert(LinearSolve(M,Vector(s)),list);
    s_ISpace:=[seq(s_in_cone_coord[ISpace[k]],k=1..nops(ISpace))];
    special_lincomb_v(s_ISpace,[seq(W[ISpace[k]],k=1..nops(ISpace))],nops(W));
end:
if check_examples() then
    ASSERT(s_ISpace([s1,s2],[[1,0],[0,1]],[1]) = [s1,0], "s_ISpace test #1");
fi;

# Basic functions
#
#
# Todd(z,t):  the function (e^(zt)*t/(1-exp(t)));
Todd:=proc(z,t);
    exp(z*t)*t/(1-exp(t));
end:
#Todd(z,t);

# Input: a symbolic expression or a number; 
# Output: This gives the formal ceil function, written CEIL(t), or the
# ceil of a number.
#
# To evaluate after substituting a symbolic t by a number, substitute CEIL by ceil.
#
# Examples: see below.
ourceil:=proc(t) local our:
    our := ceil(t);
    if type(our, numeric) then our;
    else CEIL(t);
    fi;
end:
if check_examples() then
    ASSERT(ourceil(1/2) = 1, "ourceil test #1");
    ASSERT(ourceil(-1/2) = 0, "ourceil test #2");
    ASSERT(ourceil(xyzzy) = CEIL(xyzzy), "ourceil test #3");
    ASSERT(ourceil(sqrt(2)) = 2, "ourceil test #4");
    ASSERT(eval(subs({CEIL=ceil, xyzzy=-1/2}, ourceil(xyzzy))) = 0, "ourceil test #4");
fi;

#### LATTE INTERFACE FUNCTION:
#
# Function to be substituted for the formal MOD function to get results.
#
# input:
#	a: any rational number or symbolic expression.
#	n: any rational number.
#	return the number in the half-open interval [0,n) that is equal to  a mod n
latteMod:=proc(a, n)
	local x;
	ASSERT(n > 0);
	x:=a;
	#while ( x >= n or x < 0) do
    ### I don't understand what the intention of this while loop
    ### was.  It fails in the case that `a' is symbolic, such as
    ### sqrt(2) with the error message "cannot determine if this
    ### expression is true or false". 
    ### On the other hand, `floor' works fine for symbolics, if `Digits'
    ### is large enough. --mkoeppe
	x:= x - floor(x/n)*n;
	#end;
	return x;
end:
if check_examples() then
    ASSERT(latteMod(1/3, 1) = 1/3, "latteMod test #1");
    ASSERT(latteMod(-1/3, 1) = 2/3, "latteMod test #2");
    ASSERT(latteMod(sqrt(2), 1) = sqrt(2) - 1, "latteMod test #3"); # Note the difference to what we output in fractionalpart and nfractionalpartreal!  
fi;

# Input: a symbolic expression or a number; 
# Output: This gives the formal fractional part of a function, written
#         MOD(t, 1), or the fractional part in the half-open interval [0, 1) of a number.
#
# To evaluate after substituting a symbolic s by a number, substitute
# MOD by latteMod.
#
# Examples: See below.
fractionalpart:=proc(s) local our;
    our := s - floor(s);
    if type(our, numeric) then our;
    else MOD(s, 1);
    fi;
end:
if check_examples() then
    ASSERT(fractionalpart(1/3) = 1/3, "fractionalpart test #1");
    ASSERT(fractionalpart(-1/3) = 2/3, "fractionalpart test #2");
    ASSERT(fractionalpart(xyzzy) = MOD(xyzzy, 1), "fractionalpart test #3");
    ASSERT(fractionalpart(sqrt(2)) = MOD(sqrt(2), 1), "fractionalpart test #4"); # Note that we do NOT replace it by sqrt(2) - 1; we want to keep MOD as the primitive expression for readibility.
    ASSERT(eval(subs({MOD=latteMod, xyzzy=-1/3}, fractionalpart(-1/3))) = 2/3, "fractionalpart test #5");
fi;

## helper function for nfractionalpart.
ourmod:=proc(p,q,t) local our;
    if q=1 or modp(p,q)=0 then our:=0;
    elif type(t,integer) then our:=modp(t*p,q);
    else our:=MOD(modp(p,q)*t,q);
    fi;
    RETURN(our);
end:

# Input: n (INTEGER or symbolic expression which stands for an integer), integers p, q.
# Output: A number or expression equivalent to fractional part of
#         p*n/q.
#
#         ONLY VALID FOR INTEGERS n!
#
#         If the output is a symbolic expression, it is stylized in
#         the same way of our papers "Computation of the Highest
#         Coefficients..." and "Intermediate Sums On Polyhedra:
#         Computation And Real Ehrhart Theory".  See examples below.
#
# To evaluate after substituting a symbolic n by a number, substitute
# MOD by latteMod.
#
nfractionalpart:=proc(n,p,q)
    if type(n, numeric) and not type(n, integer) then
        error "nfractionalpart may only be called with first argument symbolic or integer";
    fi;
    1/q*ourmod(p,q,n);
end:
if check_examples() then
    ASSERT(nfractionalpart(xyzzy, 13, 5) = 1/5 * MOD(3*xyzzy, 5), "nfractionalpart test #1");
    ASSERT(nfractionalpart(xyzzy, 10, 5) = 0, "nfractionalpart test #2");
    ASSERT(nfractionalpart(1,0,1) = 0, "nfractionalpart test #3");
    ASSERT(eval(subs({MOD=latteMod, xyzzy=2}, nfractionalpart(xyzzy, 13, 5))) = 1/5, "nfractionalpart test #4");
fi;

## helper function for nfractionalpartreal.
ourmodreal:=proc(p,q,t) local our;
    if type(t,integer) then our:=modp(t*p,q);fi;
    if t=0 or p=0 then our:=0;
    else our:=MOD(p*t,q);
    fi;
    our;
end:

# Input: n (number or symbolic expression), integers p, q.
# Output: A number or expression equivalent to fractional part of
#         p*n/q.
#
#         Valid for all real numbers n. 
#
#         If the output is a symbolic expression, it is stylized in
#         the same way of our papers "Computation of the Highest
#         Coefficients..." and "Intermediate Sums On Polyhedra:
#         Computation And Real Ehrhart Theory".  See examples below.
#
# To evaluate after substituting a symbolic n by a number, substitute
# MOD by latteMod.
#
nfractionalpartreal:=proc(n,p,q) local our;
    if  type(n,rational) then our:=fractionalpart(p*n/q)
    else our:=1/q*ourmodreal(p,q,n);
    fi;
    our;
end:
if check_examples() then
    ASSERT(nfractionalpartreal(xyzzy, 13, 5) = 1/5 * MOD(13*xyzzy, 5), "nfractionalpartreal test #1");
    ASSERT(nfractionalpartreal(sqrt(2), 1, 1) = MOD(sqrt(2), 1), "nfractionalpartreal test #2");
    ASSERT(nfractionalpartreal(xyzzy, 0, 17) = 0, "nfractionalpartreal test #3");
    ASSERT(nfractionalpartreal(sqrt(3), 0, 11) = 0, "nfractionalpartreal test #4");
    ASSERT(nfractionalpartreal(1,0,1) = 0, "nfractionalpartreal test #5");
    ASSERT(nfractionalpartreal(7/6, 3, 1) = 1/2, "nfractionalpartreal test #6");
    ASSERT(eval(subs({MOD=latteMod, xyzzy=2/13}, nfractionalpartreal(xyzzy, 13, 5))) = 2/5, "nfractionalpartreal test #7");
fi;

# Relative volume
#
#
# Input: W is a Cone in R^d and Cspace is a subset of [1,..,d] of cardinal k;
# Ouput: a number;
#
# Math: the volume of the Box(v[i], i not in Cspace), with respect to the intersected lattice.
#
volume_ISpace:=proc(W,ISpace) local P,M,H,MM,output;
    if ISpace=[] then output:=1;
    else P:=matrix([seq(W[ISpace[i]],i=1..nops(ISpace))]);
        M:=transpose(matrix(P));
        H:=ihermite(M);
        MM:=matrix([seq(row(H,i),i=1..nops(ISpace))]);
        output:=det(MM);
    fi;
    output;
end:
if check_examples() then
    ASSERT(volume_ISpace([[1,0],[0,1]],[1]) = 1, "volume_ISpace test #1");
fi;

# Necessary  functions to compute S_L
# Input: s a vector in R^d;  W a "Cone" in R^d; Ispace a subset of [1, 2,...,d];
# xi a variable (for a list of  d symbolic variables):

# Output: a list of two functions of xi;
# Math: #We compute integral over the affine cone s+c; sliced by subspaces parallel to ISpace of
# exp^(xi,x) ; the answer is given as [vol*exp (<q,xi>, product of linear forms]
# Representing separately the numerator and the denominator.
# Furthermore, we enter exp as a "black box" EXP(.); later on we might want to replace it.
#
functionIa:=proc(s,W,ISpace,xi)
local s_on_ISpace,d,T,i,y,r,out;
    d:=nops(W);
    s_on_ISpace:=s_ISpace(s,W,ISpace);
    if nops(ISpace)=0
    then out:=[1,1];
    else
        r:=volume_ISpace(W,ISpace);
        T:=1;
        for i from 1 to nops(ISpace) do
            y:=add(W[ISpace[i]][j]*xi[j],j=1..d);
            T:=T*y;
        od;
        T:=(-1)^(nops(ISpace))*T;
        out:=[r*EXP(add(s_on_ISpace[m]*xi[m],m=1..d)),T];
    fi;
    out;
end:
if check_examples() then
    ASSERT(functionIa([0,s], [[1,0],[1,2]],[1,2],xi) = [2*EXP(s*xi[2]), xi[1]*(xi[1]+2*xi[2])], "functionIa test #1");
fi;

# For formula b we do not enter theexponential inside.

# Input: s a vector in R^d;  W a "Cone" in R^d; Ispace a subset of [1, 2,...,d];
# xi a variable (for a list of  d symbolic variables):

# Output: a list of two functions of xi;
# Math: #We compute integral over the cone with vertex at zero  sliced by subspaces parallel to ISpace of
# exp^(xi,x) ; the answer is given as [vol, product of linear forms]
# Representing separately the numerator and the denominator.
# Furthermore, we enter exp as a "black box" EXP(); later on we might want to replace it.
#
functionIb:=proc(s,W,ISpace,xi)
local d,T,i,y,r,out;
    d:=nops(W);
    if nops(ISpace)=0
    then out:=[1,1];
    else
        r:=volume_ISpace(W,ISpace);
        T:=1;
        for i from 1 to nops(ISpace) do
            y:=add(W[ISpace[i]][j]*xi[j],j=1..d);
            T:=T*y;
        od;
        T:=(-1)^(nops(ISpace))*T;
        out:=[r,T];
    fi;
    out;
end:
if check_examples() then
    ASSERT(functionIb([0,s], [[1,0],[1,2]],[1,2],xi) = [2, xi[1]*(xi[1]+2*xi[2])], "functionIb test #1");
fi;

# Input: z =[z1,...,zd], x=[x1,x2,..,xd];  two lists of symbolic expressions (or just z,x), W a cone in R^d.
# Output: a symbolic expression, using the formal function TODD. (substitute by Todd to evaluate).
# Math: Our cone has generator w1,w2,...,wd.
# We replace x by <x,w_i> and we compute  the product of Todd(z_i,<x,w_i>);
prod_Todd:=proc(z,W,xi) local d,E,i,T,y;
    d:=nops(W);
#ASSERT(d = nops(z) and d = nops(xi),
#       "z, x, W need to be of the same length");
    T:=1;
    for i from 1 to d do
#ASSERT(nops(W[i])=nops(xi),"W[i], xi need to be of the same length");
        y:=add(W[i][j]*xi[j],j=1..nops(W[i]));
        T:=T*TODD(z[i],y);
    od;
    T;
end:
if check_examples() then
    ASSERT(prod_Todd(z,[[1,0,0],[1,2,1]],x) = TODD(z[1], x[1])*TODD(z[2], x[1]+2*x[2]+x[3]), "prod_Todd test #1");
fi;

#
#
# Input: z =[z1,...,zd], xi=[xi1,xi2,..,xid];  two lists of symbolic expression, or letters (z,xi); W a cone in R^d.
# Output: a list of two symbolic expressions [P1,Q1].
# Math: P1 is the   product of Todd(z_i,<xi,w_i>), while Q1 is  the product of the (<xi,wi>)
#
functionS:=proc(z,W,xi) local P,Q,y,i;
    P:=prod_Todd(z,W,xi);
    Q:=1;
    for i from 1 to nops(W) do
#ASSERT(nops(W[i])=nops(x),"W[i], x need to be of the same length");
        y:=add(W[i][j]*xi[j],j=1..nops(W[i]));
        Q:=Q*y;
    od;
    [P,Q];
end:
if check_examples() then
    ASSERT(functionS(z,[[1,0,0],[1,1,2],[0,5,1]],xi) 
           = 
           [TODD(z[1],xi[1]) * TODD(z[2],xi[1]+xi[2]+2*xi[3]) * TODD(z[3],5*xi[2]+xi[3]), xi[1] * (xi[1] + xi[2] + 2*xi[3]) * (5*xi[2]+xi[3])],
           "functionS test #1");
fi;

# Input: a Cone W;
#        Cspace a subset of [1..d] of cardinality k;
#        xi a letter;
#        ProjLattice := projectedlattice(W,Cspace);
# Ouput: a list of  k linear forms in variables xi[1],...xi[d].
#
#
#  Math:
# We write R^d=lin(Cspace)+lin(ISpace). We computed a basis H1,H2...H_k of the projection of the lattice Z^d in lin(Cspace).
# Thus the output is the list i <xi,H_i> where H_i are the basis vectors of the projected lattice.
#

changeofcoordinates:=proc(W,Cspace,ProjLattice,xi) local H,newxi,i,d;
    H:=ProjLattice;
    d:=nops(W[1]);
    newxi:=[];
    for i from 1 to nops(H) do
        newxi:=[op(newxi),add(xi[j]*H[i][j],j=1..d)];
    od;
    newxi;
end:
if check_examples() then
    ASSERT(changeofcoordinates([[1,0,0],[1,1,2],[0,5,1]],[1,2], projectedlattice([[1,0,0],[1,1,2],[0,5,1]], [1,2]) ,xi)
           =
           [xi[1], (1/9)*xi[2]+(2/9)*xi[3]],
           "changeofcoordinates test #1");
fi;

# THE FUNCTION (S^Ispace) for a cone. Here we sum  the integrals of e^{xi,x}
#  on slices of the cone
# parallel to  L generated by w_i with i in Ispace.
#
# WE GIVE THE TWO FORMULAE A) and B)
# THESE ARE  THE MAIN  TECHNICAL PROCEDURES.
#
#
# Input: s a vector in Q^d,  or a symbolic variable ; BUT THEN IT HAS TO BE ENTERED AS A LIST OF  d SYMBOLIC VARIABLES
# s:=[s1,s2,...,sd]; W a cone in Z^d, Ispace a subset of [1,...,d]
# xi a list of lenght d  of variable (or xi);

# The output is a function of xi[i].
#
# The subspace $L$ where we integrate is the following face of W: L is the linear span of
# <w[j]>, with j running of Ispace. (thus Ispace should be "big")
#
# Here we take out a function of s, the ceiling, the formal version of
# which is written CEIL.
#
S_Ispace_Coneformulaa:=proc(s,W,ISpace,xi) local i,ss,uni_cones,function_on_Cspace,function_on_ISpace,W_projected,WW,WWW,signuni,signL,j,Cspace,out1,out2,s_in_cone_coord,s_Cspace_in_cone_coord,s_prime_Cspace,M,newxi,dimL,g,testrank,newP,
    s_Cspace_in_lattice_coord,news,
    ProjLattice;
    Cspace:=ComplementList(ISpace,nops(W));
    ProjLattice := projectedlattice(W,Cspace);
    s_Cspace_in_lattice_coord:=projectedvertexinbasislattice(W,Cspace,ProjLattice,s);
    function_on_ISpace:=functionIa(s,W,ISpace,xi);
#from here express in terms of the basis lattice for projected cone.
    W_projected:=projectedconeinbasislattice(W,Cspace,ProjLattice):
    if W_projected=[] then
        out1:=function_on_ISpace[1]/function_on_ISpace[2];
    else
        newxi:=changeofcoordinates(W,Cspace,ProjLattice,xi);
        uni_cones:=cone_dec(W_projected);
        out1:=0;
        for j from 1 to nops(uni_cones) do
            WWW:=uni_cones[j][3];
            signuni:=uni_cones[j][1];
            ASSERT(abs(uni_cones[j][2])=1, "decomposition not unimodular");
            newP:=MatrixInverse(Transpose(Matrix(WWW)));
            news:=convert(Multiply(newP,Vector(s_Cspace_in_lattice_coord)),list); ##print(news);
            s_prime_Cspace:=[seq(ourceil(news[f]),f=1..nops(news))];
            function_on_Cspace:=functionS(s_prime_Cspace,WWW,newxi);
            out1:=out1+signuni*function_on_Cspace[1]/function_on_Cspace[2]*function_on_ISpace[1]/function_on_ISpace[2];
        od:
    fi;
    out1;
end:
if check_examples() then
    ASSERT(S_Ispace_Coneformulaa([s1,s2],[[1,0],[1,2]],[1],xi)
           =
           -TODD(CEIL(s2), (1/2)*xi[1]+xi[2])*EXP((s1-(1/2)*s2)*xi[1])/(((1/2)*xi[1]+xi[2])*xi[1]),
           "S_Ispace_Coneformulaa test #1");
fi;

# Input: s a vector in Q^d,  or a symbolic variable ; BUT THEN IT HAS TO BE ENTERED AS A LIST OF  d SYMBOLIC VARIABLES
# s:=[s1,s2,...,sd]; W a cone in Z^d, Ispace a subset of [1,...,d]
# xi a list of lenght d  of variable (or xi);

# The output is a function of xi[i].
#
# The subspace $L$ where we integrate is the following face of W: L is the linear span of
# <w[j]>, with j running of Ispace. (thus Ispace should be "big")
#
# Here we take out a function of s, the fractional part, the formal
# version of which is written as MOD( . , 1).
#
S_Ispace_Coneformulab:=proc(s,W,ISpace,xi) local i,ss,uni_cones,function_on_Cspace,function_on_ISpace,W_projected,WW,WWW,signuni,signL,j,Cspace,out1,out2,s_in_cone_coord,s_Cspace_in_cone_coord,s_small_move,M,newxi,dimL,g,testrank,newP,
    s_Cspace_in_lattice_coord,news,
    ProjLattice;
    Cspace:=ComplementList(ISpace,nops(W));
    ProjLattice := projectedlattice(W,Cspace);
    s_Cspace_in_lattice_coord:=projectedvertexinbasislattice(W,Cspace,ProjLattice,s);
    function_on_ISpace:=functionIb(s,W,ISpace,xi);
#from here express in terms of the basis lattice for projected cone.
    W_projected:=projectedconeinbasislattice(W,Cspace,ProjLattice):
    if W_projected=[] then
        out1:=function_on_ISpace[1]/function_on_ISpace[2];
    else
        newxi:=changeofcoordinates(W,Cspace,ProjLattice,xi);
        uni_cones:=cone_dec(W_projected);
        out1:=0;
        for j from 1 to nops(uni_cones) do
            WWW:=uni_cones[j][3];
            signuni:=uni_cones[j][1];
            ASSERT(abs(uni_cones[j][2])=1, "decomposition not unimodular");
            newP:=MatrixInverse(Transpose(Matrix(WWW)));
            news:=convert(Multiply(newP,Vector(s_Cspace_in_lattice_coord)),list); ##print(news);
            s_small_move:=[seq(fractionalpart(-news[f]),f=1..nops(news))];
            function_on_Cspace:=functionS(s_small_move,WWW,newxi);
            out1:=out1+signuni*function_on_Cspace[1]/function_on_Cspace[2]*function_on_ISpace[1]/function_on_ISpace[2];
        od:
    fi;
    EXP(add(s[i]*xi[i],i=1..nops(W)))*out1;
end:
if check_examples() then
    ASSERT(S_Ispace_Coneformulab([s1,s2],[[1,0],[1,2]],[1],xi)
           = 
           -EXP(s1*xi[1]+s2*xi[2])*TODD(MOD(-s2, 1), (1/2)*xi[1]+xi[2])/(((1/2)*xi[1]+xi[2])*xi[1]),
           "S_Ispace_Coneformulab test #1");
fi;

#
# WE WILL NOT USE THE FOLLOWING  PROCEDURE, AS OUR CHOICE OF REGULAR VECTOR WILL BE DONE WITH A RANDOM PROCEDURE.
# BUT WE SHOULD WHEN DETERMINING DETERMINISTICALLY A REGULAR VECTOR.
#
#
# Input: W a cone in Z^d, Cspace  a subset of [1,2,..,d]  ,x a variable.
# Output: a list of linear forms.
#  Math: this is  the forms in denominator of the   function S_Ispace_Cone(W,Cspace,x). In practice we will not use this procedure.
# This is useful to determine a "deterministic regular vector", but we will plug a random regular vector;

linindenom:=proc(W,Cspace) local YY,i,ISpace,g,WW,newx,d,a,z,cc,
    WW_projected,uni_cones,t,cleanYY,r,ProjLattice;
    d:=nops(W);
    YY:={};
    ProjLattice := projectedlattice(W,Cspace);
    WW_projected:=projectedconeinbasislattice(W,Cspace,ProjLattice):
###print(WW_projected);
    newx:=changeofcoordinates(W,Cspace,ProjLattice,x);
    uni_cones:=cone_dec(WW_projected);
    for z from 1 to nops(uni_cones) do
        cc:=uni_cones[z][3];
        for t from 1 to nops(cc) do
            YY:={op(YY), add(cc[t][s]*newx[s],s=1..nops(newx))}
        od;

    od;
    cleanYY:={};
    for r from 1 to nops(YY) do
        if member(-YY[r],cleanYY)=false then
            cleanYY:={op(cleanYY),YY[r]}
        fi;
    od;
    cleanYY;
end:
if check_examples() then
    ASSERT(linindenom([[1,0],[1,2]],[1,2])
           = {x[1], x[2], x[1]+2*x[2]},
           "linindenom test #1");
fi;

#  Approximation for a cone;
#
# Input:  s a vector in Q^d,  or a symbolic variable (but has to be entered as a list of d symbolic variables,
# W a cone, order  an integer;
# xi a list of variables.
# Output a function f(xi);
# Math:  let C=s+W;  The output is a function f(xi) such that
# the beginning of the Laurent series of f(t*xi) under dilation should  coinciding with S_C(t*xi) for order+1
#  terms:  given by formula a (with ceil functions)

# EXAMPLES ARE GIVEN AFTER;

approx_Cone_formulaa:=proc(s,W,order,xi) local output,d,j,C,a,K,KK,cc,P;
    output:=0; P:=order;
    d:=nops(W);
    if P=d then
        output:=S_Ispace_Coneformulaa(s,W,[],xi) ;
    else
        for j from 0 to P do
            C:=choose(d,j);
            cc[j]:=(-1)^(P-j)*binomial(d-j-1,d-P-1);
            for a from 1 to nops(C) do
                K:=C[a]; KK:=ComplementList(K,nops(W));
                output:=output+cc[j]*S_Ispace_Coneformulaa(s,W,KK,xi) ;
            od;
        od:
    fi;
    output;
end:
if check_examples() then
    ASSERT(approx_Cone_formulaa([s,1/2], [[1,0],[1,2]],2,xi)
           = 
           TODD(1,xi[2])*TODD(CEIL(s),xi[1])/xi[2]/xi[1]-TODD(CEIL(s),xi[1]+2*xi[2])*TODD(CEIL(2*s-1/2),-xi[2])/(xi[1]+2*xi[2])/xi[2],  ## result has not been checked
           "approx_Cone_formulaa test #1");
    ASSERT(approx_Cone_formulaa([s1,s2], [[1,0],[1,2]],1,xi)
           = -2*EXP(s1*xi[1]+s2*xi[2])/xi[1]/(xi[1]+2*xi[2])+2*TODD(CEIL(2*s1-s2),1/2*xi[1])/xi[1]*EXP(1/2*s2*xi[1]+s2*xi[2])/(-xi[1]-2*xi[2])-TODD(CEIL(s2),1/2*xi[1]+xi[2])/(1/2*xi[1]+xi[2])*EXP((s1-1/2*s2)*xi[1])/xi[1], ## result has not been checked
           "approx_Cone_formulaa test #2");
fi;           


# Input:  s a vector in Q^d,  or a symbolic variable (but has to be entered as a list of d symbolic variables,
# W a cone, order  an integer;
# xi a list of variables.
# Output a function f(xi);
# Math:  let C=s+W; A function f(xi) such that
# the beginning of the Laurent series of f(t*xi) under dilation should  coinciding with S_C(t*xi) for order+1
#  terms:  given by formula b (with fractionalparts functions)

# EXAMPLES ARE GIVEN AFTER;

approx_Cone_formulab:=proc(s,W,order,xi) local output,d,j,C,a,K,KK,cc,P;
    output:=0; P:=order;
    d:=nops(W);
    if P=d then
        output:=S_Ispace_Coneformulab(s,W,[],xi) ;
    else
        for j from 0 to P do
            C:=choose(d,j);
            cc[j]:=(-1)^(P-j)*binomial(d-j-1,d-P-1);
            for a from 1 to nops(C) do
                K:=C[a]; KK:=ComplementList(K,nops(W));
                output:=output+cc[j]*S_Ispace_Coneformulab(s,W,KK,xi) ;
            od;
        od:
    fi;
    output;
end:
if check_examples() then
    ASSERT(approx_Cone_formulab([1/2,1/2], [[1,0],[1,2]],1,xi)
           = 
           -2*EXP(1/2*xi[1]+1/2*xi[2])/xi[1]/(xi[1]+2*xi[2])+2*EXP(1/2*xi[1]+1/2*xi[2])*TODD(1/2,1/2*xi[1])/xi[1]/(-xi[1]-2*xi[2])-EXP(1/2*xi[1]+1/2*xi[2])*TODD(1/2,1/2*xi[1]+xi[2])/(1/2*xi[1]+xi[2])/xi[1], ## result has not been checked
           "approx_Cone_formulab test #1");
    ASSERT(approx_Cone_formulab([s1,s2], [[1,0],[1,2]],1,xi)
           = 
           -2*EXP(s1*xi[1]+s2*xi[2])/xi[1]/(xi[1]+2*xi[2])+2*EXP(s1*xi[1]+s2*xi[2])*TODD(MOD(-2*s1+s2,1),1/2*xi[1])/xi[1]/(-xi[1]-2*xi[2])-EXP(s1*xi[1]+s2*xi[2])*TODD(MOD(-s2,1),1/2*xi[1]+xi[2])/(1/2*xi[1]+xi[2])/xi[1], ## result has not been checked
            "approx_Cone_formulab test #2");
fi;

# ADDING THE CONES APPROXIMATIONS FOR A RATIONAL SIMPLEX;
# Input: A SIMPLEX entered as a list of d+1 rational vectors in R^d; order is an integer, xi is a variable.
# xi can also be entered as a numeric list of lenght d, but there can be then an error message (division by zero).
# OUTPUT: a function of xi;

#
# EXAMPLES ARE GIVEN AFTER;
#
#

cone_by_cone_approxi_simplex_formulaa:=proc(Simplex,order,xi) local F,W,i,st,d,S,y,P;
    F:=0; P:=order; S:=Simplex;
    d:=nops(S)-1;
    for i from 1 to nops(S) do
        W:=[seq(primitive_vector(S[j]-S[i]),j=1..i-1),seq(primitive_vector(S[j]-S[i]),j=i+1..nops(S))];
        #print(datas,S[i],W,P,xi);
        F:=F+approx_Cone_formulaa(S[i],W,P,xi);
    od:
    F:=eval(subs({TODD=Todd,EXP=exp},F));
end:
cone_by_cone_approxi_simplex_formulab:=proc(Simplex,order,xi) local F,W,i,st,d,S,y,P;
    F:=0; P:=order; S:=Simplex;
    d:=nops(S)-1;
    for i from 1 to nops(S) do
        W:=[seq(primitive_vector(S[j]-S[i]),j=1..i-1),seq(primitive_vector(S[j]-S[i]),j=i+1..nops(S))];
        F:=F+approx_Cone_formulab(S[i],W,P,xi);
    od:
    F:=eval(subs({TODD=Todd,EXP=exp},F));
end:
if check_examples() then
    ASSERT(cone_by_cone_approxi_simplex_formulab([[0,0],[1,0],[0,1]], 1,xi)
           = -1/xi[2]/xi[1]-1/(1-exp(xi[1]))/xi[2]-1/(1-exp(xi[2]))/xi[1]+exp(xi[1])/xi[1]/(-xi[1]+xi[2])+exp(xi[1])/(1-exp(-xi[1]))/(xi[1]-xi[2])+exp(xi[1])/(1-exp(-xi[1]+xi[2]))/xi[1]+exp(xi[2])/xi[2]/(xi[1]-xi[2])+exp(xi[2])/(1-exp(-xi[2]))/(-xi[1]+xi[2])+exp(xi[2])/(1-exp(xi[1]-xi[2]))/xi[2], ## result has not been checked
           "cone_by_cone_approxi_simplex_formulab test #1");
fi;

# Approximate  functions  S^L  for a  dilated  cone ns+Cone; HERE n is an integer.
# Input: n a variable,  s a numeric vector in Q^d,
# W a cone, Ispace a subset of [1,2,...d];
# xi a list of variables.
# Output a function f(n,xi);
# This is the function S^{Isplace}(ns+W)(xi), where we emphasize the dependance in n;
# We use formulab;
#

# EXAMPLE ARE GIVEN AFTER THE PROCEDURE:

dilatedS_Ispace_Cone:=proc(n,s,W,ISpace,xi) local i,ss,uni_cones,function_on_Cspace,function_on_ISpace,W_projected,WW,WWW,signuni,signL,ts,j,Cspace,out1,out2,s_in_cone_coord,s_Cspace_in_cone_coord,s_small_move,M,newxi,dimL,g,testrank,newP,dilateds,
    s_Cspace_in_lattice_coord,news,
    ProjLattice;
    #printf("### dilatedS_Ispace_Cone: W = %a, ISpace = %a\n", W, ISpace);
    Cspace:=ComplementList(ISpace,nops(W));
    ProjLattice := projectedlattice(W,Cspace);
    dilateds:=[seq(n*s[i],i=1..nops(W))];
    s_Cspace_in_lattice_coord:=projectedvertexinbasislattice(W,Cspace,ProjLattice,s);# I keep n outside;
    function_on_ISpace:=functionIb(dilateds,W,ISpace,xi);
#from here express in terms of the basis lattice for projected cone.
    W_projected:=projectedconeinbasislattice(W,Cspace,ProjLattice):
    if W_projected=[] then
        out1:=function_on_ISpace[1]/function_on_ISpace[2];
    else
        newxi:=changeofcoordinates(W,Cspace,ProjLattice,xi);
        uni_cones:=cone_dec(W_projected);
        out1:=0;
        for j from 1 to nops(uni_cones) do
            WWW:=uni_cones[j][3];
            signuni:=uni_cones[j][1];
            ASSERT(abs(uni_cones[j][2])=1, "decomposition not unimodular");
            newP:=MatrixInverse(Transpose(Matrix(WWW)));
            news:=convert(Multiply(newP,Vector(s_Cspace_in_lattice_coord)),list); #print("news",news);

            s_small_move:=[seq(nfractionalpart(N,-numer(news[f]),denom(news[f])),f=1..nops(news))];  #print("smallmove",s_small_move);
            function_on_Cspace:=functionS(s_small_move,WWW,newxi);
            out1:=out1+signuni*function_on_Cspace[1]/function_on_Cspace[2]*function_on_ISpace[1]/function_on_ISpace[2];
        od:
    fi;
    EXP(add(n*s[i]*xi[i],i=1..nops(W)))*out1;
end:

if check_examples() then
    ASSERT(dilatedS_Ispace_Cone(n,[1/2,1/2],[[1,0],[1,2]],[1],xi)
           =
           -EXP((1/2)*n*xi[1]+(1/2)*n*xi[2])*TODD((1/2)*MOD(N, 2), (1/2)*xi[1]+xi[2])/(((1/2)*xi[1]+xi[2])*xi[1]),
           "dilatedS_Ispace_Cone test #1");
fi;


random_vector:=proc(N,d) local R;
    R:=rand(N);
    [seq(R()+1,i=1..d)]:
end:




#### LATTE INTERFACE FUNCTION:
##Used by latte to find just the top k coefficients of the ehrhart polynomial.
##Because we do this by finding the polynomial per linear form, we do not print
##  out the coefficients incrementally, and hence do not have to store 
##  partial results which gives this function a smaller memory footprint than 
##  the printIncrementalEhrhartPolynomial function.
##But this function could find the entire ehrhart polynomial if asked to.
#input
#	n: symbolic variable. the coefficients are functions of n. example: 3mod(n,2)^3
#	nn: symbolic variable. The coefficients are graded by nn. example (3mod(n,2)^3 + 2)*nn^3
#	simpleCones: the polytope.
#	linearForms: list of powers of linear forms
#	d: dimension of the polytope
#	useRealDilations: ture=the polynomial can be evaluaded at rational dilations.
#	topK: find the top topK coefficients (not the topK +1) or all of them if topk=-1
#	filename: if -1, the polynomial is not saved to a file. Else, the polynomial is saved to fileName.
findEhrhartPolynomial:=proc(n,nn,simpleCones,linearForms, d,useRealDilations, topK, fileName) 
	local coef, M, ell, ehrhartPoly, mapleLinForm;
	local fPtr;
	

	ASSERT(topK > 0 or topK = -1);
	
	ehrhartPoly:=0;
	
	#loop over every linear form and collect the polynomial.
	for mapleLinForm in linearForms do
		coef:=mapleLinForm[1];
		M   :=mapleLinForm[2][1];
		ell :=mapleLinForm[2][2];

		if (topK > 0) then
			ehrhartPoly:= ehrhartPoly + coef*findEhrhartPolynomial_linearForm(n,nn,simpleCones,ell,M,d,useRealDilations, topK-1);
		else
			ehrhartPoly:= ehrhartPoly + coef*findEhrhartPolynomial_linearForm(n,nn,simpleCones,ell,M,d,useRealDilations, M+d);
		fi;
	end;
	
	if fileName <> -1 then
		fPtr:=fopen(fileName, WRITE, TEXT);
    	fprintf(fPtr, "epoly:= %a;", ehrhartPoly);
		fprintf(fPtr, "\n");
    	fclose(fPtr);
	fi;
	
	return ehrhartPoly;
end:



#Computes the top weighted ehrhart polynomial's coefficients with one power of a linear form weight
#input
#	ell: the linear form
#	M: the power of the linear form
#	d: dimension of the polytope
#	simpleCones: the polytope.
#	n: symbolic variable. the coefficients are functions of n. example: 3mod(n,2)^3
#	nn: symbolic variable. The coefficients are graded by N. example (3mod(n,2)^3 + 2)*nn^3
#	topK: compute the top topK+1 coefficients.
#return: the polynomial
findEhrhartPolynomial_linearForm:=proc(n,nn,simpleCones,ell,M,d, useRealDilations, topK) 
 local order, newOrder, xi;
 local partialF, partialSeries; 
 local totalSeries, term;
 local l, j, i, a, s, W, C, K,KK, reg, output;
 local cone, rays;
 local ehrhartPoly;
 local new_n;
 
  	new_n := `tools/gensym`('myn');
  	
 	ehrhartPoly:=0;    
    order:=min(M+d, topK);
    #order:=M+d;
    
    reg:=random_vector(5000,d);
    xi:=[seq(t*(ell[i]+epsilon*reg[i]),i=1..d)];
          
	partialF:=Array([seq(0,ll=0..order+1)]);
	partialSeries:=Array([seq(0,ll=0..order+1)]);	
    
    for j from 0 to order do
    	C:=choose(d,j);
    	output:=0;
    	#cc[j]:=(-1)^(order-j)*binomial(d-j-1,d-order-1);
    	
    	#compute the valuation over each cone.
    	for cone in simpleCones do

    		s:=cone[1]; #vertex
    		W:=[seq(primitive_vector(cone[2][j]), j=1..d)]; #rays of the cone.
	    	for a from 1 to nops(C) do
        		K:=C[a]; KK:=ComplementList(K,nops(W));
        		
        		#put the loop here: loop over every ell/xi
        		if useRealDilations then
        			output:=output + dilatedS_Ispace_Cone_real(new_n,s,W,KK,xi) ;
        		else
        			output:=output + dilatedS_Ispace_Cone(new_n,s,W,KK,xi) ;
        		fi;

        	od;
        od;
        partialF[j+1]:=eval(subs({TODD=Todd,EXP=exp},output));
        
        partialSeries[j+1]:=coeff(series(partialF[j+1],t=0,M+d+2),t,M);
        partialSeries[j+1]:=coeff(series(partialSeries[j+1],epsilon=0,d+2),epsilon,0);
    
        totalSeries:=0;
        for l from 0 to j do
    		newOrder:=j;
    		totalSeries:=totalSeries + (-1)^(newOrder-l)*binomial(d-l-1,d-newOrder-1)*partialSeries[l+1];
    	od; #for l
    	totalSeries:=coeff(totalSeries,new_n,M+d-j);
    	
    	#we need to mult. by M! because we are computing w/the weight 1/M!* ell^M
    	term:= subs({new_n=n, N=n},totalSeries)*factorial(M);
    	ehrhartPoly:=ehrhartPoly+term*nn^(M+d-j);    	
    od;

    return ehrhartPoly;
end:


#### LATTE INTERFACE FUNCTION:
##Used by latte to print ALL of the coefficients of the ehrhart polynomial coefficients incrementally
##But this function could compute just the top k ehrhart polynomial incrementally.
##If the user wants the top k, latte calls the function findEhrhartPolynomial
##  instead because it computes the polynomial per linear form, and hence has a small 
##  memory footprint because printIncrementalEhrhartPolynomial finds the 
##  current coefficient of all the linear forms at once.

#input
#	n: symbolic variable. the coefficients are functions of n. example: 3mod(n,2)^3
#	nn: symbolic variable. The coefficients are graded by N. example (3mod(n,2)^3 + 2)*nn^3
#	simpleCones: the polytope.
#	linearForms: list of powers of linear forms
#	d: dimension of the polytope
#	useRealDilations: ture=the polynomial can be evaluaded at rational dilations.
#	topK: find the top topK coefficients (note this function does not compute the top topk+1)
#	filename: if -1, the polynomial is not saved to a file. Else, the polynomial is saved to fileName.
printIncrementalEhrhartPolynomial:=proc(n,nn,simpleCones,linearForms, d, useRealDilations, topK, fileName) 
 local xi, numLinearForms, M, iLF, ell, minDegree, currentDegree, coef_current, M_current, currentDifference;
 local partialF, partialSeries, totalCoeffSum; 
 local term;
 local l, j, i, a, s, W, C, K, KK, reg, output;
 local cone, rays;
 local ehrhartPoly;
 local fPtr; #file pointer.
 local new_n;
 
 	new_n := `tools/gensym`('myn');
 	
 	#notes to self. iLF is "Index of Linear Form"
    #coef:=linearForms[iLF][1];
	#M   :=linearForms[iLF][2][1];
	#ell :=linearForms[iLF][2][2];
	
	#set up the integrand/residue direction.    
    numLinearForms:=nops(linearForms);
    M:=linearForms[1][2][1];
    xi:=Array(1..numLinearForms);
    for iLF from 1 to numLinearForms do
    	#construct a different residue direction for each linear form.
    	
    	reg:=random_vector(5000,d);
    	ell :=linearForms[iLF][2][2];
    	xi[iLF]:=[seq(t*(ell[i]+epsilon*reg[i]),i=1..d)];    	
    	
    	#find the largest degree of all the linear forms.
    	M:=max(M, linearForms[iLF][2][1]);
    od;

    #set up the output file if needed
    if fileName <> -1 then
    	fPtr:=fopen(fileName, WRITE, TEXT);
    	fprintf(fPtr, "epoly:= ");
    fi;
    
    ehrhartPoly:=0;

    #topK should be a natural number or -1.
    ASSERT(topK = -1 or topK > 0);
    
 	if topK > 0 then 
    	minDegree:=M+d-(topK-1);
    else
    	minDegree:=0;
    fi;
    
    #these are just placeholders for the dilatedS_Ispace_Cone_real() output.
	partialF:=0; 
	partialSeries:=Array(1..numLinearForms, 1..(M+d+1), 0); 
	#holds the partial taylor series of the rational functions.
		#Note: The index of array's start at ONE in maple.
		#the rows corresponds to the linear forms. Columns corresponds to the  difference between the current power and the current degree
		#example: partialSeries[1][2+1] is the rational functions when we are summing over choose(d,2)

		
    
	
	#for every degree of the ehrhart polynomial we are going to compute
	for currentDegree from M+d to minDegree by -1 do
	
		#for every linear form.
		for iLF from 1 to numLinearForms do
		    coef_current:=linearForms[iLF][1];
		    M_current   :=linearForms[iLF][2][1];
		    
			#does this linear form add a coefficient of currentDegree to the polynomial?
			if currentDegree > M_current + d then
				next; #skip this linear form.
			fi;
		
			#figure out how far the currentDegree is from M_current + d
			currentDifference:=M_current + d - currentDegree; 
			
			C:=choose(d, currentDifference);
			
			output:=0; #output is the "short rational" functions for this order
			
			
			#compute the integrand/valuation over each cone
			for cone in simpleCones do
			
    			s:=cone[1]; #vertex
    			W:=[seq(primitive_vector(cone[2][j]), j=1..d)]; #rays of the cone.
    			for a from 1 to nops(C) do
        			K:=C[a]; 
        			KK:=ComplementList(K,nops(W));
        		
        			
        			if useRealDilations then
        				output:=output + coef_current*factorial(M_current)*dilatedS_Ispace_Cone_real(new_n,s,W,KK,xi[iLF]) ;
        			else
        				output:=output + coef_current*factorial(M_current)*dilatedS_Ispace_Cone(new_n,s,W,KK,xi[iLF]) ;
        			fi;
        		
        		od; 
			od; #for every simple cone
			
			#partialF is now the sum of rational functions for every cone and the current linear form.
			partialF:=eval(subs({TODD=Todd,EXP=exp},output));
			
			#find the series expansion of partialF.        
			partialSeries[iLF, currentDifference+1]:=coeff(series(partialF,t=0,M_current+d+2),t,M_current);			
			partialSeries[iLF, currentDifference+1]:=coeff(series(partialSeries[iLF, currentDifference+1],epsilon=0,d+2),epsilon,0);

			
        od; #for every linear form
        
        #we are ready to find the coefficient of currentDegree
        
        totalCoeffSum:=0;
		for iLF from 1 to numLinearForms do
		    M_current   :=linearForms[iLF][2][1];
		    
			#does this linear form add a coefficient of currentDegree to the polynomial?
			if currentDegree > M_current + d then
				next; #skip this linear form.
			fi;
			
			#figure out how far the currentDegree is from M_current + d
			currentDifference:=M_current + d - currentDegree; 
		
			for l from 0 to currentDifference do
				
	    		totalCoeffSum:=totalCoeffSum + (-1)^(currentDifference-l)*binomial(d-l-1,d-currentDifference-1)*coeff(partialSeries[iLF, l+1], new_n, currentDegree);
    		od; #for l
    	
		od; #for every linear form.
        term:= subs({new_n=n},totalCoeffSum);
    	term:= subs({N=n},totalCoeffSum);
    	ehrhartPoly:=ehrhartPoly+term*nn^(currentDegree);
    	printf("+ %a\n", term*nn^(currentDegree));

    	#also print to output file if needed
    	if fileName <> -1 then
    		fprintf(fPtr, "\\ \n+ %a", term*nn^(currentDegree));
    	fi;
	od; #for every degree of the polynomial.
	
	
   	#close the file if needed.
  	if fileName <> -1 then
  		fprintf(fPtr, ";\n");
    	fclose(fPtr);
   	fi;

    
	#finaly, we are done!    
    return ehrhartPoly;
end:




##################################################################
### CODE FOR REAL DILATIONS
##################################################################

######################################################################""""
# Approximate  functions  S^L  for a  dilated  cone ns+Cone; HERE n is a real.
# Input: n a variable,  s a numeric vector in Q^d,
# W a cone, Ispace a subset of [1,2,...d];
# xi a list of variables.
# Output a function f(n,N,xi);

#

# Output a function f(n,N,xi);  where we emphasize the dependance in n;
# N is the same than n, but here the function of N are perodic;
# I did that,  as we will need to pick  up a polynomial term in n, while N are then considered as constants;
# EXAMPLE IS GIVEN  AFTER;
#
#
#
#
#

dilatedS_Ispace_Cone_real:=proc(n,s,W,ISpace,xi) local i,ss,uni_cones,function_on_Cspace,function_on_ISpace,W_projected,WW,WWW,signuni,signL,ts,j,Cspace,out1,out2,s_in_cone_coord,s_Cspace_in_cone_coord,s_small_move,M,newxi,dimL,g,testrank,newP,dilateds,
    s_Cspace_in_lattice_coord,news,
    ProjLattice;
    Cspace:=ComplementList(ISpace,nops(W));
    ProjLattice := projectedlattice(W,Cspace);
    dilateds:=[seq(n*s[i],i=1..nops(W))];
    s_Cspace_in_lattice_coord:=projectedvertexinbasislattice(W,Cspace,ProjLattice,s);# I keep n outside;
    function_on_ISpace:=functionIb(dilateds,W,ISpace,xi);
#from here express in terms of the basis lattice for projected cone.
    W_projected:=projectedconeinbasislattice(W,Cspace,ProjLattice):
    if W_projected=[] then
        out1:=function_on_ISpace[1]/function_on_ISpace[2];
    else
        newxi:=changeofcoordinates(W,Cspace,ProjLattice,xi);
        uni_cones:=cone_dec(W_projected);
        out1:=0;
        for j from 1 to nops(uni_cones) do
            WWW:=uni_cones[j][3];
            signuni:=uni_cones[j][1];
            ASSERT(abs(uni_cones[j][2])=1, "decomposition not unimodular");
            newP:=MatrixInverse(Transpose(Matrix(WWW)));
            news:=convert(Multiply(newP,Vector(s_Cspace_in_lattice_coord)),list); #print("news",news);

            s_small_move:=[seq(nfractionalpartreal(N,-numer(news[f]),denom(news[f])),f=1..nops(news))];  #print("smallmove",s_small_move);
            function_on_Cspace:=functionS(s_small_move,WWW,newxi);
            out1:=out1+signuni*function_on_Cspace[1]/function_on_Cspace[2]*function_on_ISpace[1]/function_on_ISpace[2];
        od:
    fi;
    EXP(add(n*s[i]*xi[i],i=1..nops(W)))*out1;
end:
if check_examples() then
    ASSERT(dilatedS_Ispace_Cone_real(n,[0,0],[[1,0],[1,2]],[1],xi)
           =
           -EXP(0)*TODD(0,1/2*xi[1]+xi[2])/(1/2*xi[1]+xi[2])/xi[1], ## result has not been checked
           "dilatedS_Ispace_Cone_real test #1");
fi;

######################################################################""""



##EXAMPLES:
#VERIFICATION FOR APPROXIMATION;
random_rational_vector:=proc(N,d) local R;
    R:=rand(N);
    [seq(R()/(R()+1),i=1..d)]:
end:

randomaffinecone:=proc(N,d) local S,i,c;
    c:=[];
    for i from 1 to d do
        c:=[op(c),primitive_vector(random_vector(N,d))];
    od;
    [random_rational_vector(10,d),c];
end:
#randomaffinecone(10,4);

random_rational_simplex:=proc(N,d) local S,i,c;
    c:=[];
    for i from 1 to d+1 do
        c:=[op(c),random_rational_vector(N,d)];
    od;end:

####################################################################
## EXAMPLES OF FORMULA a) versus FORMULA b)


#S_Ispace_Coneformulaa([s[1],s[2]],[[1,1],[1,-1]],[1],xi);

#S_Ispace_Coneformulab([s[1],s[2]],[[1,1],[1,-1]],[1],xi);



#S_Ispace_Coneformulaa([s[1],s[2],s[3]],[[1,1,1],[1,-1,0],[1,1,0]],[1],xi);

#S_Ispace_Coneformulab([s[1],s[2],s[3]],[[1,1,1],[1,-1,0],[1,1,0]],[1],xi);

###############################################################
##########################

#approx_Cone_formulab([s[1],s[2],s[3]],[[1,1,1],[1,-1,0],[1,1,0]],0,xi);
#approx_Cone_formulab([s[1],s[2],s[3]],[[1,1,1],[1,-1,0],[1,1,0]],1,xi);
#approx_Cone_formulab([s[1],s[2],s[3]],[[1,1,1],[1,-1,0],[1,1,0]],2,xi);


##########################
# VERIFICATION OF THE PROPERTY OF APPROXIMATION; BY EVALUATING IN A RANDOM VECTOR;
#SEEMS CORRECT;

# TODO: add automatic tests. --mkoeppe

checkapprox:=proc(s,Cone,k) local FFa,FFb,Fd,xx,xi;
    xi:=random_vector(100,nops(Cone));print(xi);
    xx:=[seq(t*xi[i],i=1..nops(Cone))];print(xx);
    FFa:=eval(subs({TODD=Todd,EXP=exp},approx_Cone_formulaa(s,Cone,k,xx)));#print(FFa);
    FFb:=eval(subs({TODD=Todd,EXP=exp},approx_Cone_formulab(s,Cone,k,xx)));#print(FFb);
    Fd:=eval(subs({TODD=Todd,EXP=exp},approx_Cone_formulab(s,Cone,nops(Cone),xx)));# print(Fd);
    [simplify(series(FFa-Fd,t=0,nops(Cone)+2)),simplify(series(FFb-Fd,t=0,nops(Cone)+2))];
end:

#checkapprox([1/2,1/2,1/3,1/4],[[1,1,-1,1],[1,2,0,0],[1,2,3,4],[1,4,5,7]],3);
#Coneindex2:=[[1,1,1],[1,-1,0],[1,1,0]];

#A0:=eval(subs({TODD=Todd,EXP=exp},approx_Cone_formulab([1/2,1/3,1/4],Coneindex2,0,xi)));
#A1:=eval(subs({TODD=Todd,EXP=exp},approx_Cone_formulab([1/2,1/3,1/4],Coneindex2,1,xi)));
#A2:=eval(subs({TODD=Todd,EXP=exp},approx_Cone_formulab([1/2,1/3,1/4],Coneindex2,2,xi)));
#A3:=eval(subs({TODD=Todd,EXP=exp},approx_Cone_formulab([1/2,1/3,1/4],Coneindex2,3,xi)));



#####################################################################
### LATTE INTERFACE HELPER FUNCTIONS:
####################################################################

#### LATTE INTERFACE FUNCTION:
# input:
#   lau: a triple [l, a, u], where a is symbolic and l, u are rational
#        such that a is assumed to lie in the interval between l and u.
#        (l does not have to be smaller than u.)
#   n:   base (modulus), a positive rational number
#   return an expression in a for the number between [0,n) that is
#        equal to a mod n, without using "mod".  If no such expression
#        can be given because [l, u] is too large, signal an error.
latteIntervalMod:=proc(lau, n)
	local r, l, a, u, lnfloor, unfloor;
	ASSERT(n > 0, "modulus must be positive");
    l, a, u := op(lau); 
    if (u < l) then  # this easily happens by multiplying lau by a negative scalar.
        u, l := l, u;
    end if;
    lnfloor := floor(l/n);
    unfloor := floor(u/n);
    if lnfloor <> unfloor then
        error "The range [%1, %2] is too large to simplify MOD(%3, %4).  Provide more precision.", l, u, a, n;
    end if;
    return a - lnfloor*n;
end:

#### LATTE INTERFACE FUNCTION:
# Return a rational interval l, u representing the range of numbers
# represented by a Maple software floating point number.
floatToInterval := proc(f)
    local mantissa, exponent;
    ASSERT(type(f, sfloat), cat("not a software float: ", f));
    mantissa, exponent := op(f);
    return (mantissa - 1/2 ) * 10^exponent, (mantissa + 1/2) * 10^exponent;
end:
if check_examples() then
    ASSERT([floatToInterval(2.0)] = [39/20, 41/20], "floatToInterval test #1");
    ASSERT([floatToInterval(2.00)] = [399/200, 401/200], "floatToInterval test #2");
fi;

#### LATTE INTERFACE FUNCTION:
evaluateEhrhart:=proc(epoly, a)
    local l, u, r;
    if type(a, float) then
        l, u := floatToInterval(a);
        try
             r := expand(eval(subs({N=n, n=[l, n, u], MOD=latteIntervalMod}, epoly)));
             # The result is a polynomial in n, which has to be constant if
             # epoly was the full Ehrhart quasi-polynomial.
             if type(r, rational) then
                 printf("# Exact answer (assuming that all provided digits of the floating-point dilation factor were correct):\n");
                 return r;
             else
                 printf("# Given Ehrhart quasi-polynomial was not complete; evaluating using floating point.  Increase Digits if you want more precision.\n");
                 #Digits := ilog10(1 + abs(op(1, a))) + 3;
                 return evalf(subs({n=a}, r));
             end if;
        catch "The range [%1, %2] is too large to simplify MOD(%3, %4).  Provide more precision.":
            printf("# The precision of the given floating point dilation factor is not large enough to allow exact computation.  Resorting to floating point evaluation.  Increase Digits if you want more precision in this evaluation.\n");
            #Digits := ilog10(1 + abs(op(1, a))) + 3;
            return eval(subs({N=a, n=a, MOD=latteMod}, epoly));
        end;
    else
        # `expand' helps in case that a is a symbolic expression
        return expand(eval(subs({N=a, n=a, MOD=latteMod}, epoly)));
    end if;
end:
if check_examples() then
    # example from cube_3 for real dilations; see LattE manual.
    epoly := 8*N^3+(12-24*MOD(n,1))*N^2+(6-12*MOD(n,1)+359/125245152*(-5521248*MOD(n,1)+2760624)*MOD(n,1)-503/125245152*(-6102432*MOD(n,1)+3051216)*MOD(n,1)-289/62622576*(145728*MOD(n,1)-72864)*MOD(n,1)-359/125245152*(-2623400*MOD(n,1)+1311700)*MOD(n,1)+289/62622576*(-187200*MOD(n,1)+93600)*MOD(n,1)+5/434879*(393984*MOD(n,1)-196992)*MOD(n,1)+1081/125245152*(1120600*MOD(n,1)-560300)*MOD(n,1)-1081/125245152*(-2606688*MOD(n,1)+1303344)*MOD(n,1)+937/125245152*(-2358432*MOD(n,1)+1179216)*MOD(n,1)+503/125245152*(-2042216*MOD(n,1)+1021108)*MOD(n,1)-5/434879*(-435456*MOD(n,1)+217728)*MOD(n,1)-937/125245152*(872344*MOD(n,1)-436172)*MOD(n,1))*N+1-1700/431*MOD(n,1)+1652/431*MOD(n,1)^2-18/434879*(-2018*MOD(n,1)+1009)*MOD(n,1)^2-1/62064*(1728+10368*MOD(n,1)^2-10368*MOD(n,1))*MOD(n,1)+18/434879*(2018*MOD(n,1)-1009)*MOD(n,1)^2+1/62064*(185761/3+371522*MOD(n,1)^2-371522*MOD(n,1))*MOD(n,1)-431/144*(-1081/1009*MOD(n,1)+1081/2018)*MOD(n,1)^2+1/144*(-72/1009*(2018*MOD(n,1)-1009)*MOD(n,1)-1241209/6054-1023265/1009*MOD(n,1)^2+1095913/1009*MOD(n,1))*MOD(n,1)-431/144*(-937/1009*MOD(n,1)+937/2018)*MOD(n,1)^2+1/144*(72/1009*(-2018*MOD(n,1)+1009)*MOD(n,1)+805321/6054+1023265/1009*MOD(n,1)^2-950617/1009*MOD(n,1))*MOD(n,1)+431/2018*(-1081/72*MOD(n,1)+1081/144)*MOD(n,1)^2-1/2018*(-1009/72*(-144*MOD(n,1)+72)*MOD(n,1)+1241209/432+1023265/72*MOD(n,1)^2-1095913/72*MOD(n,1))*MOD(n,1)-431/2018*(937/72*MOD(n,1)-937/144)*MOD(n,1)^2+1/2018*(1009/72*(-144*MOD(n,1)+72)*MOD(n,1)+805321/432+1023265/72*MOD(n,1)^2-950617/72*MOD(n,1))*MOD(n,1)-36/1009*(-1440/431*MOD(n,1)+720/431)*MOD(n,1)^2+1/2018*(-1009/431*(862*MOD(n,1)-431)*MOD(n,1)-2508479/2586-1203842/431*MOD(n,1)^2+1638721/431*MOD(n,1))*MOD(n,1)+1/72*(-72/431*(862*MOD(n,1)-431)*MOD(n,1)-284041/2586-190945/431*MOD(n,1)^2+221977/431*MOD(n,1))*MOD(n,1)+36/1009*(578/431*MOD(n,1)-289/431)*MOD(n,1)^2-1/2018*(1009/431*(862*MOD(n,1)-431)*MOD(n,1)+100795/2586-1203842/431*MOD(n,1)^2+768963/431*MOD(n,1))*MOD(n,1):
    ASSERT(evaluateEhrhart(epoly, 0) 
           = 1,
           "evaluateEhrhart test #1");
    ASSERT(evaluateEhrhart(epoly, 1)
           = 27,
           "evaluateEhrhart test #2");
    ASSERT(evaluateEhrhart(epoly, 3/2)
           = 27,
           "evaluateEhrhart test #3");
    ASSERT(evaluateEhrhart(epoly, 2-1/1000000)
           = 27, 
           "evaluateEhrhart test #4");
    ASSERT(evaluateEhrhart(epoly, 2)
           = 125);
    ## TODO: Add more examples; suppress output in the real case... --mkoeppe
fi:

#### LATTE INTERFACE FUNCTION:
#input:
#	simpleCones is a list of d+1 vertex-ray cones in the form
#		[[[vertex],[[ray1], ..., [rayn]]], ...]
#		We assume the cones are simple, that is, we have the tangent-cones of a simplex.
#	d: integer, dimension.
#
# return a list of just the d+1 vertices.
tangentConesToSimplex:=proc(simpleCones, d)
local Simplex, cone;

	Simplex:=[];
	for cone in simpleCones do
		Simplex:=[op(Simplex), cone[1]];
	end;

	ASSERT(nops(Simplex) = d+1); #make sure we have a simplex.
		
	return Simplex;
end:


#### LATTE INTERFACE FUNCTION:
#input
#	Simplex: list of d+1 verticies
#return maple-list of the tangent-cones.
SimplexToTangentCones:=proc(Simplex)
	local simpleCones;
	local cone, rays, d, i;
	
	simpleCones:=[];
	cone:=[];
	d:=nops(Simplex)-1;
	
	for i from 1 to d+1 do
		rays:= [seq(primitive_vector(Simplex[j]-Simplex[i]), j=1..(i-1)), seq(primitive_vector(Simplex[j]-Simplex[i]), j=(i+1)..d+1)] ;
		cone:=[Simplex[i], rays ];
		simpleCones:=[op(simpleCones), cone];
	end;
	
	return simpleCones;
end:


#####################################################################
### Functions I want to delete
#####################################################################

# Input: n a variable,  s a numeric vector in Q^d,
# W a cone, order is an integer;
# xi a list of variables.
# Output a function f(n,xi);
# This is the sulm of the approximate  function S^{Ispace}(ns+W)(xi), where we emphasize the dependance in n;

# EXAMPLE IS GIVEN  AFTER;

#
dilated_approxi_cone:=proc(n,s,W,order,xi) local output,d,j,C,a,K,KK,cc,P;
    #printf("##### dilated_approxi_cone: order = %d\n", order);
    output:=0;
    d:=nops(W);
    if order=d then
        # Fast path; general code below handles this case just fine.
        output:=dilatedS_Ispace_Cone(n,s,W,[],xi);
    else
        for j from 0 to order do
            C:=choose(d,j);
            #print("choose j, order, c", j, order, C);
            cc[j]:=(-1)^(order-j)*binomial(d-j-1,d-order-1);
            for a from 1 to nops(C) do
                K:=C[a]; KK:=ComplementList(K,nops(W));
                output:=output+cc[j]*dilatedS_Ispace_Cone(n,s,W,KK,xi) ;
            od;
        od:
    fi;
    output;
end:
if check_examples() then
    ASSERT(dilated_approxi_cone(n,[1/2,1/2],[[1,0],[1,2]], 1,xi)
           = -2*EXP((1/2)*n*xi[1]+(1/2)*n*xi[2])/(xi[1]*(xi[1]+2*xi[2]))+2*EXP((1/2)*n*xi[1]+(1/2)*n*xi[2])*TODD((1/2)*MOD(N, 2), (1/2)*xi[1])/(xi[1]*(-xi[1]-2*xi[2]))-EXP((1/2)*n*xi[1]+(1/2)*n*xi[2])*TODD((1/2)*MOD(N, 2), (1/2)*xi[1]+xi[2])/(((1/2)*xi[1]+xi[2])*xi[1]),
           "dilated_approxi_cone test #1");
fi;

# Input: n a variable,  Simplex  a numeric rational simplex ; given by a list of  rational  vectors in Q^d
# order is an integer;
# xi a list of variables. xi can be numeric but then there can be an error message;
# Output a function f(n,N,xi);
# This is the sum of the approximate  functions  over the tangent cones  for the dilated simplex nS; where we emphasize the dependance in n;
# N is the same as n, but here the functions of N are periodic;
# I did that, as we will need to pick up a polynomial term in n, while N are then considered as constants;
# EXAMPLE IS GIVEN  AFTER;

#
ApproxEhrhartSimplexgeneric:=proc(n,Simplex,order,xi) local F,W,i,st,d,S,y,P;
    F:=0;  S:=Simplex;
    d:=nops(S)-1;
    for i from 1 to nops(S) do
        W:=[seq(primitive_vector(S[j]-S[i]),j=1..i-1),seq(primitive_vector(S[j]-S[i]),j=i+1..nops(S))];
        F:=F+dilated_approxi_cone(n,S[i],W,order,xi) ;
    od:
    F:=eval(subs({TODD=Todd,EXP=exp},F));
    
    return F;
end:
if check_examples() then
    ASSERT(ApproxEhrhartSimplexgeneric(n,[[0,0],[1/2,0],[0,1/2]], 1,xi)
           = -1/xi[2]/xi[1]-1/(1-exp(xi[1]))/xi[2]-1/(1-exp(xi[2]))/xi[1]+exp(1/2*n*xi[1])/xi[1]/(-xi[1]+xi[2])+exp(1/2*n*xi[1])*exp(-1/2*MOD(N,2)*xi[1])/(1-exp(-xi[1]))/(xi[1]-xi[2])+exp(1/2*n*xi[1])/(1-exp(-xi[1]+xi[2]))/xi[1]+exp(1/2*n*xi[2])/xi[2]/(xi[1]-xi[2])+exp(1/2*n*xi[2])*exp(-1/2*MOD(N,2)*xi[2])/(1-exp(-xi[2]))/(-xi[1]+xi[2])+exp(1/2*n*xi[2])/(1-exp(xi[1]-xi[2]))/xi[2], ## result has not been checked
           "ApproxEhrhartSimplexgeneric test #1");
fi;

#WARNING; THIS WORKS ONLY IF ell is generic;
# Input; n is a variable, Simplex is a rational simplex, ell is a linear form fiven as a numeric list of d+1 rational numbers; M is in integer, m is an integer.
# The ouput is  a periodic function of n;
# Math: the output is the m coefficient Ehrhart polynomial E(n S, ell^M)
# Here I did not employ a deformation vector, so the procedure might return: error; diviasion by zero.
#
TopEhrhartweightedluckyell:=proc(n,Simplex,ell,M,m) local d,order,xx,AA,CC;
    d:=nops(Simplex)-1;
    order:=M+d-m;
    xx:=[seq(t*ell[i],i=1..d)];
    AA:=ApproxEhrhartSimplexgeneric(n,Simplex,order,xx);
    CC:=coeff(coeff(series(AA,t=0,M+d+2),t,M),n,m);
    subs({N=n},CC);
end:
# Input; n is a variable, Simplex is a rational simplex, ell is a linear form fiven as a numeric list of d+1 rational numbers; M is in integer, m is an integer.
# The ouput is  a periodic function of n;
# Math: the output is the m coefficient Ehrhart polynomial E(n S, ell^M)
# Here we employ a random deformation vector, so if the procedure might return: error: division by zero. RERUN:
#
#
TopEhrhartweighted:=proc(n,Simplex,ell,M,m) local d,order,xx,AA,CCt,CCeps,CCn,reg;
    d:=nops(Simplex)-1;
    order:=M+d-m;
    reg:=random_vector(5000,d);
    xx:=[seq(t*(ell[i]+epsilon*reg[i]),i=1..d)];
    AA:=ApproxEhrhartSimplexgeneric(n,Simplex,order,xx);
    CCt:=coeff(series(AA,t=0,M+d+2),t,M); #print(CCt);
    CCeps:=coeff(series(CCt,epsilon=0,d+2),epsilon,0);
    CCn:=coeff(CCeps,n,m);
    return subs({N=n},CCn);
end:
# Input; n is a variable, Simplex is a rational simplex, ell is a linear form fiven as a numeric list of d+1 rational numbers; M is in integer, m is an integer.
# The ouput is  a polynomial with coefficients  periodic function of n;
# Math: the output is Ehrhart polynomial E(n S, ell^M)
# Here we employ a random deformation vector, so if the procedure might return: error; diviasion by zero. RERUN:
#
#
#
CompleteEhrhartweighted:=proc(n,Simplex,ell,M) local d;
    d:=nops(Simplex)-1;
    add(TopEhrhartweighted(n,Simplex,ell,M,m)*n^m,m=0..M+d);
end:

#####################################################################
### real Functions I want to delete
#####################################################################


dilated_approxi_cone_real:=proc(n,s,W,order,xi) local output,d,j,C,a,K,KK,cc,P;
    output:=0;
    d:=nops(W);
    if order=d then
        output:=dilatedS_Ispace_Cone_real(n,s,W,[],xi);
    else
        for j from 0 to order do
            C:=choose(d,j);
            cc[j]:=(-1)^(order-j)*binomial(d-j-1,d-order-1);
            for a from 1 to nops(C) do
                K:=C[a]; KK:=ComplementList(K,nops(W));
                output:=output+cc[j]*dilatedS_Ispace_Cone_real(n,s,W,KK,xi) ;
            od;
        od:
    fi;
    output;
end:

ApproxEhrhartSimplexgeneric_real:=proc(n,Simplex,order,xi) local F,W,i,st,d,S,y,P;
    F:=0;  S:=Simplex;
    d:=nops(S)-1;
    for i from 1 to nops(S) do
        W:=[seq(primitive_vector(S[j]-S[i]),j=1..i-1),seq(primitive_vector(S[j]-S[i]),j=i+1..nops(S))];
        F:=F+dilated_approxi_cone_real(n,S[i],W,order,xi) ;
    od:
    F:=eval(subs({TODD=Todd,EXP=exp},F));
end:


#WARNING; THIS WORKS ONLY IF ell is generic;
TopEhrhartweightedluckyell_real:=proc(n,Simplex,ell,M,m) local d,order,xx,AA,CC;
    d:=nops(Simplex)-1;
    order:=M+d-m;
    xx:=[seq(t*ell[i],i=1..d)];
    AA:=ApproxEhrhartSimplexgeneric_real(n,Simplex,order,xx);
    CC:=coeff(coeff(series(AA,t=0,M+d+2),t,M),n,m);
    subs({N=n},CC);
end:

TopEhrhartweighted_real:=proc(n,Simplex,ell,M,m) local d,order,xx,AA,CCt,CCeps,CCn,reg;
    d:=nops(Simplex)-1;
    order:=M+d-m;
    reg:=random_vector(5000,d);
    xx:=[seq(t*(ell[i]+epsilon*reg[i]),i=1..d)];
    AA:=ApproxEhrhartSimplexgeneric_real(n,Simplex,order,xx);
    CCt:=coeff(series(AA,t=0,M+d+2),t,M); #print(CCt);
    CCeps:=coeff(series(CCt,epsilon=0,d+2),epsilon,0);
    CCn:=coeff(CCeps,n,m);
    subs({N=n},CCn);
end:
CompleteEhrhartweighted_real:=proc(n,nn,Simplex,ell,M) local d;
    d:=nops(Simplex)-1;
    add(TopEhrhartweighted_real(n,Simplex,ell,M,mTopEhrhartweightedPoly_real)*nn^m,m=0..M+d);
end:


