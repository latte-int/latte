with(linalg):with(LinearAlgebra):with(combinat):
kernelopts(assertlevel=1):       ### Enable checking ASSERTions

# PEDAGOCICAL PROGRAM FOR COMPUTING EXAMPLES FOR ARTICLE:
# HIGHEST EHRHART  COEFFICIENTS;
# Version of November 10 -2010 : I have added the Ehrhart polynomial over the reals.
#
# I PUT SOME EXAMPLES AFTER THE PROCEDURE.
#
#
#
# HERE IS THE LIST OF WHAT THE PROGRAM DOES.
#
#
# I
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# PROGRAM FOR  FORMULAE A and b for S^{L^I}(s+c) .
# This function  is expressed in terms of the "black box functions" TODD and EXP.
# TODD(s,x) representing the function e^(sx) x/(1-e^x); IF we want to evaluate S_C_k(x) at a regular element reg; the command eval and subs (TODD=Todd, EXP=exp) should be used.
#
# The commands are  EITHER;
# S_Ispace_Coneformulaa:=proc(vertex,cone, ISpace,xi);
# S_Ispace_Coneformulab:=proc(vertex,cone, ISpace,xi);
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
#
# Programs on lists: addition on lists, complement of a list, sublist,etc...
#
# Input: K a subset of integers, L a list.The output takes the elements of the list L in the position of the list K
#
insert:=proc(K,L) local out;
    out:=[seq(L[K[i]],i=1..nops(K))];
end:
#The output is the Complement  List, within the list [1,..,d]
ComplementList:=proc(K,d);
    RETURN([seq (`if` (member(i,K)=false, i, op({})),i=1..d)]);
end:
#The output is the Complement  List, within the list [a[1],..,a[d]]
GeneralComplementList:=proc(K,L)local d;d:=nops(L);
    RETURN([seq (`if` (member(L[i],K)=false, L[i], op({})),i=1..d)]);
end:
#GeneralComplementList([2,3],[1,2,3,7]);
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
# Example: #primitive_vector([0,-1/2])->[0,-1];

#
primitive_vector:=proc(A) local d,n,g;
    d:=nops(A);
    n:=ilcm(seq(denom(A[i]),i=1..d));
    g:=igcd(seq(n*A[i],i=1..d));if g<>0 then
                                    [seq(n*A[i]/g,i=1..d)];else [seq(n*A[i],i=1..d)];fi;
end:
ortho_basis:=proc(d) local i,v;
    for i from 1 to d do
        v[i]:=[seq(0,j=1..i-1),1,seq(0,j=i+1..d)]
    od;[seq(v[j],j=1..d)];
end:
fracpart:=proc(x); x-floor(x);
end:

fracpart(1);
#  Signed decomposition into unimodular cones
# A "simplicial cone" is a list of  d linearly independent  vectors in Z^d, sometimes assumed primitive.
#
# short_vector(A)
#
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
#
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
# #
# # Input   G  is a  "simplicial cone"
# # Output consists of 2 elements:
# #              V is a vector in Z^d.
# #               L=[ Lplus,Lminus,Lzero] is a partition of [1..d] into three sublists,
# #               according to the signs of the entries of the vector V. in the basis G.
good_vector:=proc(G) local n,A,Ainverse,B,sho,V,L;
    n:=nops(G);
    A:=Transpose(Matrix(G));
    Ainverse:=MatrixInverse(A);
    B:=[seq(convert(Ainverse[1..n,i],list),i=1..n)];
    sho:=short_vector(B);
    V :=[seq(add(G[j][i]*sho[j],j=1..n),i=1..n)];
    L:= sign_entries_vector(sho);
    [V,L];
end:
# # signed_decomp(eps,G,v,L)
#
# # Input :  eps = 1 or -1
# #             G  is a  "simplicial cone"
# #              V is a vector of dim d
# #              L= [ Lplus,Lminus,Lzero] is a partition of [1..d] into three sublists,
# #
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

            detC := Determinant(Matrix(C));
            Csigned:=[eps*(-1)^(i+kplus),detC,C];

            if abs(detC)>1 then
                Nonuni:=[op(Nonuni),Csigned] else Uni:=[op(Uni),Csigned];
            fi;
        od;
    fi;

    if kminus>0 then
        for i from 1 to kminus do
            C:=[seq(G[Lplus[j]],j=1..kplus),-v,seq(-G[Lminus[j]],j=1..i-1),seq(G[Lminus[j]],j=i+1..kminus),seq(G[Lzero[j]],j=1..kzero)];

            detC := Determinant(Matrix(C));
            Csigned:=[eps*(-1)^(i+1),detC,C];

            if abs(detC)>1 then
                Nonuni:=[op(Nonuni),Csigned] else Uni:=[op(Uni), Csigned];
            fi;
        od;
    end if;
    [Nonuni,Uni];
end:
#
#
# # good_cone_dec(eps,G)
# #  Input: eps = 1 or -1
# #             G  is a  simplicial cone
# #
# #  Output:  two lists [Nonuni,Uni] as in procedure signed_decomp:
# #
good_cone_dec:=proc(eps,G) local n,A,R,Output;
    n:=nops(G);  A:=Matrix([seq(G[i],i=1..n)]);
    if abs(Determinant(A))=1 then  Output:=[[],[[eps,Determinant(A),G]]];
    else R:=good_vector(G);
        Output:=signed_decomp(eps,G,R[1],R[2]);
    fi;
end:
# # more_decomposition_in_cones(cones)
#
# # Input:  cones =[cones[1],cones[2]] as in procedure signed_decomp
# # Output: [Newnonuni,Newuni] as in procedure signed_decomp
#
#
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
# #
# # Input:  G is a "simplicial cone"
# # Output: A list of  terms [eps,detG,G] where
# "               eps =1 or -1,
# #               detG is an integer ( hopefully 1 or -1),
# #               G  is a  "simplicial cone", (hopefully unimodular)
# #
cone_dec:=proc(G) local seed, i,ok;
    if G=[] then RETURN([[1,1,[]]]);fi:
    seed:=good_cone_dec(1,G);
    ok:=0;
    i:=1; while ok=0  do
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
# Example: projectedvector([[1,0,0],[0,1,2],[0,1,0]],[3],[0,0,1])->[0,-1/2,0];
projectedvector:=proc(W,Cspace,b) local M,S,j,v,V,m;
    M:=transpose(matrix([seq(W[i],i=1..nops(W))]));
    S:=linsolve(M,b);
    m:=det(M);
    for j from 1 to nops(W) do
        v[j]:=add(S[Cspace[i]]*W[Cspace[i]][j],i=1..nops(Cspace));
    od:
    V:=[seq(v[j],j=1..nops(W))];
end:
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
# EXAMPLE:
#projectedlattice([[1,3,0],[0,1,0],[0,0,2]],[1,3])-># [[1, 3, 0], [0, 0, 1]];
#
#
#
#
#
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
#projectedlattice([[1,3,0],[0,1,0],[0,0,2]],[1,3]);
# Projected cone and projected vertex (expressed in the lattice basis)
# Input: W is a Cone in Z^d and Cspace is a subset of [1,..,d] of cardinal k;
#  Output: A "Cone" in Z^k;

# Be careful: our input must have integral coordinates.
# The ouput then will have integral coordinates.
#
#
# Here W is the cone and we are projecting W over lin( Cspace) and expressing it in term of the basis H_1,H_2,...,H_k of projectedlattice(W,Cspace).
#  Example: projectedconeinbasislattice([[1,1,0],[0,1,0],[0,0,2]],[1,3])→[[1,0],[0,1]]
projectedconeinbasislattice:=proc(W,Cspace) local P,M,output,i,F;
    P:=projectedlattice(W,Cspace);
    M:=transpose(matrix([seq(P[i],i=1..nops(P))]));
    output:=[];
    for i from 1 to nops(Cspace) do
        F:=convert(linsolve(M,Vector(W[Cspace[i]])),list);
        output:=[op(output),primitive_vector(F)];
    od;
    output;
end:

#projectedconeinbasislattice([[1,1,0],[0,1,0],[0,0,2]],[1,3]);
#
# #Input; W a Cone in Z^d;
# Cspace a subset of [1,2,..d] of cardinal k;
# s a vector in R^d with rational coordinates (or symbolic coordinates);
# #Ouput: a vector in R^k with rational coordinates (or symbolic coordinates);

# Math: Here W is the cone and we are projecting V over lin( Cspace)  using  V:=lin(Cspace) oplus
#  lin(ISpace). We express the projection of s
# with respect to the basis of the projected lattice. If the ouput is [a1,a2], this means that our
# projected vertex is s_Cspace=a1*H1+a2*H2 where H1,H2 is the basis of the projected lattice computed before.
#
#
# Example: projectedvertexinbasislattice([[1,0,0],[0,2,1],[0,1,1]],[1,3],[s1,s2,s3]) ->[s1, 2*s3-s2];;
#
projectedvertexinbasislattice:=proc(W,Cspace,s) local m,P,M,output,i,F;
    P:=projectedlattice(W,Cspace);
    if Cspace=[] then RETURN([]);fi;
    M:=Transpose(Matrix([seq(P[i],i=1..nops(P))]));
    F:=convert(LinearSolve(M,Vector(projectedvector(W,Cspace,s))),list);
    output:=F;
end:
# Input: s a vector in R^d with rational coordinates (or symbolic).
# W a cone in Z^d;
# Ispace a subset of [1,2,...,d];
# Output:  a vector in R^d
#
# Math: We decompose V in lin(CSpace) oplus lin (ISpace), (with CSpace spanned by the v|i] in the complementary indices of Ispace );
#  and here we write s=s_Cspace+s_(ISpace): Here the output is
# s_(ISpace);
# Example: s_ISpace([s1,s2],[[1,0],[0,1]],[1])->[s1,0];
s_ISpace:=proc(s,W,ISpace) local M,s_in_cone_coord,s_ISpace;
    M:=Matrix([seq(Vector([W[i]]),i=1..nops(W))]);
    s_in_cone_coord:=convert(LinearSolve(M,Vector(s)),list);
    s_ISpace:=[seq(s_in_cone_coord[ISpace[k]],k=1..nops(ISpace))];
    special_lincomb_v(s_ISpace,[seq(W[ISpace[k]],k=1..nops(ISpace))],nops(W));
end:
#s_ISpace([s1,s2],[[1,0],[0,1]],[1]);
# Basic functions
#
#
# Todd(z,t):  the function (e^(zt)*t/(1-exp(t)));
Todd:=proc(z,t);
    exp(z*t)*t/(1-exp(t));
end:
#Todd(z,t);
# Input: a symbolic variable or a number; Output:This gives the formal ceil function or the ceil of a number.
# EXAMPLE : ourceil(t)->ceil(t); ourceil(1/2)->1;
ourceil:=proc(t) local u:
    if type(t,rational) then
        ceil(t);
    else ceil(t);
    fi;
end:

#ourceil(t);

# Input: a symbolic variable or a number; Output:This gives the formal  fractionalpart of a function or the fractional part (betwwen 0 and 1) of a number.
# EXAMPLE : fractionalpart(t)->{t};fractionalpart(3/2);->1/2
#
fractionalpart:=proc(s) local our,T;
    if  type(s,rational) then our:=s-floor(s);
    else our:={s};
    fi;
    RETURN(our);
end:
fractionalpart(3/2);
fmod:=proc(p,q,t) local u:
    u:=modp(p,q);##print(u);
    MOD(u*t,q);
end:
ourmod:=proc(p,q,t) local our,T;
    if q=1 or modp(p,q)=0 then our:=0;
    elif type(t,integer) then our:=modp(t*p,q);
    else our:=MOD(modp(p,q)*t,q);
    fi;
    RETURN(our);
end:
nfractionalpart:=proc(n,p,q) local our;
    if  type(n,rational) then our:=fractionalpart(p*n/q)
    else our:=1/q*ourmod(p,q,n);
    fi;
    our;
end:

ourmodreal:=proc(p,q,t) local our,T;
    if type(t,integer) then our:=modp(t*p,q);fi;
    if t=0 or p=0 then our:=0;
    else our:=MOD(p*t,q);
    fi;
    our;
end:
nfractionalpartreal:=proc(n,p,q) local our;
    if  type(n,rational) then our:=fractionalpart(p*n/q)
    else our:=1/q*ourmodreal(p,q,n);
    fi;
    our;
end:
#nfractionalpartreal(1,0,1);

# Relative volume
#
#
# Input: W is a Cone in R^d and Cspace is a subset of [1,..,d] of cardinal k;
# Ouput: a number;
#
# Math; the volume of the Box(v[i], i not in Cspace), with respect to the intersected lattice.
# Example: volume_ISpace([[1,1],[0,1]],[1])->1;
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
#volume_ISpace([[1,0],[0,1]],[1]);
# Necessary  functions to compute S_L
# Input: s a vector in R^d;  W a "Cone" in R^d; Ispace a subset of [1, 2,...,d];
# xi a variable (for a list of  d symbolic variables):

# Output: a list of two functions of xi;
# Math: #We compute integral over the affine cone s+c; sliced by subspaces parallel to ISpace of
# exp^(xi,x) ; the answer is given as [vol*exp (<q,xi>, product of linear forms]
# Representing separately the numerator and the denominator.
# Furthermore, we enter exp as a "black box" EXP(.); later on we might want to replace it.
# #Example functionIa([0,s], [[1,0],[1,2]],[1,2],xi)-> [2*EXP(s*xi[2]), xi[1]*(xi[1]+2*xi[2])];
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
# For formula b we do not enter theexponential inside.

# Input: s a vector in R^d;  W a "Cone" in R^d; Ispace a subset of [1, 2,...,d];
# xi a variable (for a list of  d symbolic variables):

# Output: a list of two functions of xi;
# Math: #We compute integral over the cone with vertex at zero  sliced by subspaces parallel to ISpace of
# exp^(xi,x) ; the answer is given as [vol, product of linear forms]
# Representing separately the numerator and the denominator.
# Furthermore, we enter exp as a "black box" EXP(); later on we might want to replace it.
# #Example functionIb([0,s], [[1,0],[1,2]],[1,2],xi)-> [2, xi[1]*(xi[1]+2*xi[2])];
#
#
#
#
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
functionIb([0,s], [[1,0],[1,2]],[1,2],xi);
# Input: z =[z1,...,zd], x=[x1,x2,..,xd];  two lists of symbolic expressions (or just z,x), W a cone in R^d.
# Output: a symbolic expression.
# Math: Our cone has generator w1,w2,...,wd.
# We replace x by <x,w_i> and we compute  the product of Todd(z_i,<x,w_i>);
# #Example: prod_Todd(z,[[1,0,0],[1,2,1]],x)-> TODD(z[1], x[1])*TODD(z[2], x[1]+2*x[2]+x[3]);
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
prod_Todd(z,[[1,0,0],[1,2,1]],xi);
#
#
# Input: z =[z1,...,zd], xi=[xi1,xi2,..,xid];  two lists of symbolic expression, or letters (z,xi); W a cone in R^d.
# Output: a list of two symbolic expressions [P1,Q1].
# Math: P1 is the   product of Todd(z_i,<xi,w_i>), while Q1 is  the product of the (<xi,wi>)
# Example: functionS(z,[[1,0,0],[1,1,2],[0,5,1]],xi) ->[TODD(z[1],xi[1]) TODD(z[2],xi[1]+xi[2]+2 xi[3]) TODD(z[3],5 xi[2]+xi[3]),xi[1] (xi[1]+xi[2]+2 xi[3]) (5 xi[2]+xi[3])];
# ;
#
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
functionS(z,[[1,0,0],[1,1,2],[0,5,1]],xi);
#
#
# Input: a Cone W;  Cspace a subset of [1..d] of cardinal k; xi a letter:
# Ouput: a list of  k linear forms in variables xi[1],...xi[d].
#
#
#  Math:
# We write R^d=lin(Cspace)+lin(ISpace). We computed a basis H1,H2...H_k of the projection of the lattice Z^d in lin(Cspace).
# Thus the output is the list i <xi,H_i> where H_i are the basis vectors of the projected lattice.
#
# Example: changeofcoordinates([[1,0,0],[1,1,2],[0,5,1]],[1,2],xi)->[xi[1], (1/9)*xi[2]+(2/9)*xi[3]];;

changeofcoordinates:=proc(W,Cspace,xi) local H,newxi,i,d;
    H:=projectedlattice(W,Cspace);
    d:=nops(W[1]);
    newxi:=[];
    for i from 1 to nops(H) do
        newxi:=[op(newxi),add(xi[j]*H[i][j],j=1..d)];
    od;
    newxi;
end:
changeofcoordinates([[1,0,0],[1,1,2],[0,5,1]],[1,2],xi);
# THE FUNCTION (S^Ispace) for a cone. Here we sum  the integrals of e^{xi,x}
#  on slices of the cone
# parallel to  L generated by w_i with i in Ispace.
# WE GIVE THE TWO FORMULAE A) and B)
# THESE ARE  THE MAIN  TECHNICAL PROCEDURES.
#
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
# Here we take out a function of s, ceil
# EXAMPLE: S_Ispace_Coneformulaa([s1,s2],[[1,0],[1,2]],[1],xi)->  -TODD(ceil(s2), (1/2)*xi[1]+xi[2])*EXP((s1-(1/2)*s2)*xi[1])/(((1/2)*xi[1]+xi[2])*xi[1]);
S_Ispace_Coneformulaa:=proc(s,W,ISpace,xi) local i,ss,uni_cones,function_on_Cspace,function_on_ISpace,W_projected,WW,WWW,signuni,signL,j,Cspace,out1,out2,s_in_cone_coord,s_Cspace_in_cone_coord,s_prime_Cspace,M,newxi,dimL,g,testrank,newP,
    s_Cspace_in_lattice_coord,news;
    Cspace:=ComplementList(ISpace,nops(W));
    s_Cspace_in_lattice_coord:=projectedvertexinbasislattice(W,Cspace,s);
    function_on_ISpace:=functionIa(s,W,ISpace,xi);
#from here express in terms of the basis lattice for projected cone.
    W_projected:=projectedconeinbasislattice(W,Cspace):
    if W_projected=[] then
        out1:=function_on_ISpace[1]/function_on_ISpace[2];
    else
        newxi:=changeofcoordinates(W,Cspace,xi);
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


S_Ispace_Coneformulaa([s1,s2],[[1,0],[1,2]],[1],xi);
# Input: s a vector in Q^d,  or a symbolic variable ; BUT THEN IT HAS TO BE ENTERED AS A LIST OF  d SYMBOLIC VARIABLES
# s:=[s1,s2,...,sd]; W a cone in Z^d, Ispace a subset of [1,...,d]
# xi a list of lenght d  of variable (or xi);

# The output is a function of xi[i].
#
# The subspace $L$ where we integrate is the following face of W: L is the linear span of
# <w[j]>, with j running of Ispace. (thus Ispace should be "big")
#
# Here we take out a function of s, fractionalpart
# ##EXAMPLE: S_Ispace_Coneformulab([s1,s2],[[1,0],[1,2]],[1],xi)-> -EXP(s1*xi[1]+s2*xi[2])*TODD({-s2}, (1/2)*xi[1]+xi[2])/(((1/2)*xi[1]+xi[2])*xi[1]);
#
S_Ispace_Coneformulab:=proc(s,W,ISpace,xi) local i,ss,uni_cones,function_on_Cspace,function_on_ISpace,W_projected,WW,WWW,signuni,signL,j,Cspace,out1,out2,s_in_cone_coord,s_Cspace_in_cone_coord,s_small_move,M,newxi,dimL,g,testrank,newP,
    s_Cspace_in_lattice_coord,news;
    Cspace:=ComplementList(ISpace,nops(W));
    s_Cspace_in_lattice_coord:=projectedvertexinbasislattice(W,Cspace,s);
    function_on_ISpace:=functionIb(s,W,ISpace,xi);
#from here express in terms of the basis lattice for projected cone.
    W_projected:=projectedconeinbasislattice(W,Cspace):
    if W_projected=[] then
        out1:=function_on_ISpace[1]/function_on_ISpace[2];
    else
        newxi:=changeofcoordinates(W,Cspace,xi);
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


S_Ispace_Coneformulab([s1,s2],[[1,0],[1,2]],[1],xi);
#
# WE WILL NOT USE THE FOLLOWING  PROCEDURE, AS OUR CHOICE OF REGULAR VECTOR WILL BE DONE WITH A RANDOM PROCEDURE.
# BUT WE SHOULD WHEN DETERMINING DETERMINISTICALLY A REGULAR VECTOR.
#
#
# Input: W a cone in Z^d, Cspace  a subset of [1,2,..,d]  ,x a variable.
# Output: a list of linear forms.
#  Math: this is  the forms in denominator of the   function S_Ispace_Cone(W,Cspace,x). In practice we will not use this procedure.
# This is useful to determine a "deterministic regular vector", but we will plug a random regular vector;
# EXAMPLE linindenom([[1,0],[1,2]],[1,2])->{x[1], x[2], x[1]+2*x[2]};
#

linindenom:=proc(W,Cspace) local YY,i,ISpace,g,WW,newx,d,a,z,cc,
    WW_projected,uni_cones,t,cleanYY,r;
    d:=nops(W);
    YY:={};
    WW_projected:=projectedconeinbasislattice(W,Cspace):
###print(WW_projected);
    newx:=changeofcoordinates(W,Cspace,x);
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


#linindenom([[1,0],[1,2]],[1,2]);
#  Approximation for a cone;
#
#
#
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
#approx_Cone_formulaa([s,1/2], [[1,0],[1,2]],2,xi);
#approx_Cone_formulaa([s1,s2], [[1,0],[1,2]],1,xi);
#
#
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
#approx_Cone_formulab([1/2,1/2], [[1,0],[1,2]],1,xi);
#approx_Cone_formulab([s1,s2], [[1,0],[1,2]],1,xi);
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
cone_by_cone_approxi_simplex_formulab([[0,0],[1,0],[0,1]], 1,xi);

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
    s_Cspace_in_lattice_coord,news;
    Cspace:=ComplementList(ISpace,nops(W));
    dilateds:=[seq(n*s[i],i=1..nops(W))];
    s_Cspace_in_lattice_coord:=projectedvertexinbasislattice(W,Cspace,s);# I keep n outside;
    function_on_ISpace:=functionIb(dilateds,W,ISpace,xi);
#from here express in terms of the basis lattice for projected cone.
    W_projected:=projectedconeinbasislattice(W,Cspace):
    if W_projected=[] then
        out1:=function_on_ISpace[1]/function_on_ISpace[2];
    else
        newxi:=changeofcoordinates(W,Cspace,xi);
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


#dilatedS_Ispace_Cone(n,[1/2,1/2],[[1,0],[1,2]],[1],xi); #Ouput;-EXP((1/2)*n*xi[1]+(1/2)*n*xi[2])*TODD((1/2)*MOD(N, 2), (1/2)*xi[1]+xi[2])/(((1/2)*xi[1]+xi[2])*xi[1])

# Input: n a variable,  s a numeric vector in Q^d,
# W a cone, order is an integer;
# xi a list of variables.
# Output a function f(n,xi);
# This is the sulm of the approximate  function S^{Isplace}(ns+W)(xi), where we emphasize the dependance in n;

# EXAMPLE IS GIVEN  AFTER;

#
dilated_approxi_cone:=proc(n,s,W,order,xi) local output,d,j,C,a,K,KK,cc,P;
    output:=0;
    d:=nops(W);
    if order=d then
        output:=dilatedS_Ispace_Cone(n,s,W,[],xi);
    else
        for j from 0 to order do
            C:=choose(d,j);
            cc[j]:=(-1)^(order-j)*binomial(d-j-1,d-order-1);
            for a from 1 to nops(C) do
                K:=C[a]; KK:=ComplementList(K,nops(W));
                output:=output+cc[j]*dilatedS_Ispace_Cone(n,s,W,KK,xi) ;
            od;
        od:
    fi;
    output;
end:

#dilated_approxi_cone(n,[1/2,1/2],[[1,0],[1,2]], 1,xi); # Ouput:-2*EXP((1/2)*n*xi[1]+(1/2)*n*xi[2])/(xi[1]*(xi[1]+2*xi[2]))+2*EXP((1/2)*n*xi[1]+(1/2)*n*xi[2])*TODD((1/2)*MOD(N, 2), (1/2)*xi[1])/(xi[1]*(-xi[1]-2*xi[2]))-EXP((1/2)*n*xi[1]+(1/2)*n*xi[2])*TODD((1/2)*MOD(N, 2), (1/2)*xi[1]+xi[2])/(((1/2)*xi[1]+xi[2])*xi[1])
# Input: n a variable,  Simplex  a numeric rational simplex ; given by a list of  rational  vectors in Q^d
# order is an integer;
# xi a list of variables. xi can be numeric but then there can be an error message;
# Output a function f(n,N,xi);
# This is the sum of the approximate  functions  over the tangent cones  for the dilated simplex nS; where we emphasize the dependance in n;
# N is the same than n, but here the function of N are perodic;
# I did that,  as we will need to pick  up a polynomial term in n, while N are then considered as constants;
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
end:



#ApproxEhrhartSimplexgeneric(n,[[0,0],[1/2,0],[0,1/2]], 1,xi);
random_vector:=proc(N,d) local R;
    R:=rand(N);
    [seq(R()+1,i=1..d)]:
end:

#WARNING; THIS WORKS ONLY IF ell is generic;
# Input; n is a variable, Simplex is a rational simplex, ell is a linear form fiven as a numeric list of d+1 rational numbers; M is in integer, m is an integer.
# The ouput is  a periodic function of n;
# Math: the output is the m coefficient Ehrhart polynomial E(n S, ell^M)
# Here I did not employ a deformation vector, so the procedure might return: error; diviasion by zero.
#
TopEhrhartweightedluckyell:=proc(n,Simplex,ell,M,m) local d,order,xx,AA,CC;
    d:=nops(Simplex)-1;
    order:=M+nops(Simplex)-m;
    xx:=[seq(t*ell[i],i=1..d)];
    AA:=ApproxEhrhartSimplexgeneric(n,Simplex,order,xx);
    CC:=coeff(coeff(series(AA,t=0,M+d+2),t,M),n,m);
    subs({N=n},CC);
end:
# Input; n is a variable, Simplex is a rational simplex, ell is a linear form fiven as a numeric list of d+1 rational numbers; M is in integer, m is an integer.
# The ouput is  a periodic function of n;
# Math: the output is the m coefficient Ehrhart polynomial E(n S, ell^M)
# Here we employ a random deformation vector, so if the procedure might return: error; diviasion by zero. RERUN:
#
#
TopEhrhartweighted:=proc(n,Simplex,ell,M,m) local d,order,xx,AA,CCt,CCeps,CCn,reg;
    d:=nops(Simplex)-1;
    order:=M+nops(Simplex)-m;
    reg:=random_vector(5000,d);
    xx:=[seq(t*(ell[i]+epsilon*reg[i]),i=1..d)];
    AA:=ApproxEhrhartSimplexgeneric(n,Simplex,order,xx);
    CCt:=coeff(series(AA,t=0,M+d+2),t,M); #print(CCt);
    CCeps:=coeff(series(CCt,epsilon=0,d+2),epsilon,0);
    CCn:=coeff(CCeps,n,m);
    subs({N=n},CCn);
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
# New functions: --Matthias
# Compute the highest (k+1) terms.
TopEhrhartweightedPoly:=proc(n,Simplex,ell,M,given_k) local k,d,order,xx,AA,CCt,CCeps,CCn,reg;
    d:=nops(Simplex)-1;
    k := min(M+d, given_k);
    order:=M+nops(Simplex)-k;
    reg:=random_vector(5000,d);
    xx:=[seq(t*(ell[i]+epsilon*reg[i]),i=1..d)];
    AA:=ApproxEhrhartSimplexgeneric(n,Simplex,order,xx);
    CCt:=coeff(series(AA,t=0,M+d+2),t,M); #print(CCt);
    CCeps:=coeff(series(CCt,epsilon=0,d+2),epsilon,0);
    CCn:=add(coeff(CCeps,n,m) * N^m, m=M+d-k..M+d);
    subs({N=n},CCn);
end:
#### LATTE INTERFACE FUNCTION:
printTopEhrhartweightedPoly:=proc(n,Simplex,ell,M,k)
    printf("%a\n", TopEhrhartweightedPoly(n,Simplex,ell,M,k));
end:
# Incrementally compute and print all terms 0f the weighted Ehrhart
# polynomial.  User can interrupt when computation takes too long.
#### LATTE INTERFACE FUNCTION:
printIncrementalEhrhartweightedPoly:=proc(n,Simplex,ell,M) local d;
    local m, term, poly;
    d:=nops(Simplex)-1;
    m := M+d;
    term := TopEhrhartweighted(n,Simplex,ell,M,m)*n^m;
    printf("%a\n", term);
    poly := term;
    for m from M+d-1 to 0 by -1 do
        term := TopEhrhartweighted(n,Simplex,ell,M,m)*n^m;
        printf("+ %a\n", term);
        poly := poly + term;
    od;
    printf("## Evaluation at n=1: %a\n",
           eval(subs(n=1,MOD=modp, poly)));
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
    s_Cspace_in_lattice_coord,news;
    Cspace:=ComplementList(ISpace,nops(W));
    dilateds:=[seq(n*s[i],i=1..nops(W))];
    s_Cspace_in_lattice_coord:=projectedvertexinbasislattice(W,Cspace,s);# I keep n outside;
    function_on_ISpace:=functionIb(dilateds,W,ISpace,xi);
#from here express in terms of the basis lattice for projected cone.
    W_projected:=projectedconeinbasislattice(W,Cspace):
    if W_projected=[] then
        out1:=function_on_ISpace[1]/function_on_ISpace[2];
    else
        newxi:=changeofcoordinates(W,Cspace,xi);
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


#dilatedS_Ispace_Cone_real(n,[0,0],[[1,0],[1,2]],[1],xi);
dilated_approxi_cone_real:=proc(n,s,W,order,xi) local output,d,j,C,a,K,KK,cc,P;
    output:=0;
    d:=nops(W);
    if order=d then
        output:=dilatedS_Ispace_Cone(n,s,W,[],xi);
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



random_vector:=proc(N,d) local R;
    R:=rand(N);
    [seq(R()+1,i=1..d)]:
end:

#WARNING; THIS WORKS ONLY IF ell is generic;
TopEhrhartweightedluckyell_real:=proc(n,Simplex,ell,M,m) local d,order,xx,AA,CC;
    d:=nops(Simplex)-1;
    order:=M+nops(Simplex)-m;
    xx:=[seq(t*ell[i],i=1..d)];
    AA:=ApproxEhrhartSimplexgeneric_real(n,Simplex,order,xx);
    CC:=coeff(coeff(series(AA,t=0,M+d+2),t,M),n,m);
    subs({N=n},CC);
end:
TopEhrhartweighted_real:=proc(n,Simplex,ell,M,m) local d,order,xx,AA,CCt,CCeps,CCn,reg;
    d:=nops(Simplex)-1;
    order:=M+nops(Simplex)-m;
    reg:=random_vector(5000,d);
    xx:=[seq(t*(ell[i]+epsilon*reg[i]),i=1..d)];
    AA:=ApproxEhrhartSimplexgeneric_real(n,Simplex,order,xx);
    CCt:=coeff(series(AA,t=0,M+d+2),t,M); #print(CCt);
    CCeps:=coeff(series(CCt,epsilon=0,d+2),epsilon,0);
    CCn:=coeff(CCeps,n,m);
    subs({N=n},CCn);
end:
CompleteEhrhartweighted_real:=proc(n,Simplex,ell,M) local d;
    d:=nops(Simplex)-1;
    add(TopEhrhartweighted_real(n,Simplex,ell,M,m)*n^m,m=0..M+d);
end:
# New functions: --Matthias
# Compute the highest (k+1) terms.
TopEhrhartweightedPoly_real:=proc(n,Simplex,ell,M,given_k) local k,d,order,xx,AA,CCt,CCeps,CCn,reg;
    d:=nops(Simplex)-1;
    k := min(M+d, given_k);
    order:=M+nops(Simplex)-k;
    reg:=random_vector(5000,d);
    xx:=[seq(t*(ell[i]+epsilon*reg[i]),i=1..d)];
    AA:=ApproxEhrhartSimplexgeneric_real(n,Simplex,order,xx);
    CCt:=coeff(series(AA,t=0,M+d+2),t,M); #print(CCt);
    CCeps:=coeff(series(CCt,epsilon=0,d+2),epsilon,0);
    CCn:=add(coeff(CCeps,n,m) * N^m, m=M+d-k..M+d);
    subs({N=n},CCn);
end:
#### LATTE INTERFACE FUNCTION:
printTopEhrhartweightedPoly_real:=proc(n,Simplex,ell,M,k)
    printf("%a\n", TopEhrhartweightedPoly_real(n,Simplex,ell,M,k));
end:
# Incrementally compute and print all terms 0f the weighted Ehrhart
# polynomial.  User can interrupt when computation takes too long.
#### LATTE INTERFACE FUNCTION:
printIncrementalEhrhartweightedPoly_real:=proc(n,Simplex,ell,M) local d;
    local m, term, poly;
    d:=nops(Simplex)-1;
    m := M+d;
    term := TopEhrhartweighted_real(n,Simplex,ell,M,m)*n^m;
    printf("%a\n", term);
    poly := term;
    for m from M+d-1 to 0 by -1 do
        term := TopEhrhartweighted_real(n,Simplex,ell,M,m)*n^m;
        printf("+ %a\n", term);
        poly := poly + term;
    od;
    printf("## Evaluation at n=1: %a\n",
           eval(subs(n=1,MOD=modp, poly)));
end:

######################################################################""""



EXAMPLES:
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
randomaffinecone(10,4);
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




