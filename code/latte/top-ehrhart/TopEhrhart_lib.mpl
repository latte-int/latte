
##########################################################################
# Program for computing  the first k+1 COEFFICIENTs of  t^(d-k) of the Ehrhart
# polynomial of a full dimensional SIMPLEX WITH INTEGRAL VERTICES in R^d
#
# October 2010
# Authors: V. Baldoni, N. Berline, J. De Loera, M. Koeppe, and M. Vergne.
# with help from B. Dutra and G. Pinto.

# SUMMARY AND GLOBAL DESCRIPTION:

# Let S be a full dimensional simplex in R^d with integral
# vertices. We consider its Ehrhart polynomial
# Eh(S)(t)=Cardinal(t S\cap Z^d).
# We know that
# Eh(S)(t):=vol(S)t^d+1/2 vol( boundary(S))t^(d-1)+
# e_[d-2]t^(d-2)+ e_[d-3]t^(d-3)+... +1

# The aim  of this program is to compute the coefficient
# in t^(d-j) of Eh(S)(t), from j=0..k.
# We follow the algorithm (and the notation) of the companion
# paper (Baldoni-Berline-Koeppe-DeLoera-Vergne)
#
#
# so the algorithm presented here is runs in polynomial time with
# respect to the data, provided k is fixed.
# IMPORTANT: the dimension d, and the vertices of S  are variables.

# The two MAIN COMMANDS are

#coeff_dminusk_Eh(S,k);

# here the input is the simplex S and the codimension k  This returns
# a single coefficient.
#  k is an integer between 0 and d and  the output is the coefficient of degree t^(d-k) of Eh(S)(t).

# EXAMPLE
# We input the standard triangle S:=[[0,0],[1,0],[0,1]];
# It is well known that the Ehrhart polunomial is
# (t+1)(t+2)/2=t^2/2+3/2 t+1;

#coeff_dminusk_Eh([[0,0],[1,0],[0,1]],0);
#coeff_dminusk_Eh([[0,0],[1,0],[0,1]],1);
#coeff_dminusk_Eh([[0,0,0],[1,0,0],[0,1,0],[0,0,1]],3);

# Another faster alternative is to compute all the coefficients UP to
# k directly.


#################################################
# General guide to notation:

# d always denotes the dimension.
# We work in R^d:

# A vector in R^d: a list of  d rational numbers:

# A SIMPLEX in R^d:
# entered as a list of d+1
# vectors.

# A Cone in R^d : a list of d vectors of lenght d: (we will only consider simplicial cones)
# when we say Cone in Z^d we mean that the vectors have integral coordinates;

# A Signed cone:  [epsilon, Cone] where epsilon is -1 or 1.

###############################################################################



# This program is extracted from the more general program:
# total_approx_weighted_Eh;
# The program  "total_approx_weighted_Eh" is (from our point of view) easier to read and to understand.
# Here however, we tried to  just do the necessary calculation for the
# term we want.

################################ BEGIN OF CODE ################################

with(linalg):
with(LinearAlgebra):
with(combinat):
kernelopts(assertlevel=1):      ### Enable checking ASSERTions


####################################################
# The following are minor utility subroutines to manipulate
# lists: addition on lists, complement of a list, sublist,etc...

# Input: K a subset of integers, L a list.The output takes the elements of the list L in the position of the list K
insert:=proc(K,L) local out;
    out:=[seq(L[K[i]],i=1..nops(K))];
end:

#The output is the Complement  List, within the list [1,..,d]
ComplementList:=proc(K,d);
    RETURN([seq (`if` (member(i,K)=false, i, op({})),i=1..d)]);
end:

#The output is the Complement  List, within the list [a[1],..,a[d]]
#i.e. GeneralComplementList([2,3],[1,2,3,7]);
GeneralComplementList:=proc(K,L)local d;d:=nops(L);
    RETURN([seq (`if` (member(L[i],K)=false, L[i], op({})),i=1..d)]);
end:

special_lincomb_v:=proc(a,v,n) local out;
    ASSERT(nops(a)=nops(v)," the number of coefficients and vectors do not match");
    if v=[]
    then out:=[seq(0,i=1..n)];
    else
        out:=[seq(add(a[i]*v[i][j],i=1..nops(v)),j=1..nops(v[1]))];
    fi;
    out;
end:

###########################
# Computing PRIMITIVE VECTOR along a ray
# Input: A :a vector with rational coordinates.
# Output: A vector with integral coordinates:
# Math: the primitive vector on the half line R^+A;
# Example: #primitive_vector([0,-1/2])->[0,-1];

primitive_vector:=proc(A) local d,n,g;
    d:=nops(A);
    n:=ilcm(seq(denom(A[i]),i=1..d));
    g:=igcd(seq(n*A[i],i=1..d));
    if g<>0 then
        [seq(n*A[i]/g,i=1..d)];else [seq(n*A[i],i=1..d)];
    fi;
end:
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

###############################
# Projections:
# Input: W is a list of vectors  of V , [v[1],..v[d]], of lenght d.
# iota =[i[1]..,i[s]], is a list of integers, b is a vector of lenght d.
# Output: a vector of lenght d,.
# Math:
# We decompose the space V in lin(II)+lin(iotac) where lin(iota) of the vectors v[i], i in iota, and lin(iota_c) of the vectors in the complement indices. We project a vector b on lin(iota)
# Thus we write b=b_iota+b_iota_c;
# Our output is b_iota;
# Example: projectedvector([[1,0,0],[0,1,2],[0,1,0]],[3],[0,0,1])->[0,-1/2,0];
projectedvector:=proc(M_inverse, W,iota,b) local S,j,v,V;

    #printf("forGregtests.268: size(M)= %d, %d\n", rowdim(M), coldim(M));
    S:=multiply(M_inverse,b);
    #m:=det(M);
    for j from 1 to nops(W) do
        v[j]:=add(S[iota[i]]*W[iota[i]][j],i=1..nops(iota));
    od:
    V:=[seq(v[j],j=1..nops(W))];
end:


# Projected lattice
# # Input:  W=[v1,v2,.., vd];  a "Cone"  in  R^d;
# BE CAREFUl: The vectors in W must have integral coordinates.
#  iota a subset of [1,2..d] of cardinal k;
# # Output a list [H1,H2,...,Hk] of vectors in R^d with k terms.
# projectedlattice:
# Math: we
# decompose V in lin(iota)+lin(iota_c);
#  we project the standard lattice (that is Ze[1]+..+Ze[d], that is  Z[1,0,0..0]+... Z.[0,0,0..,1]])
# on lin[iota] which is a  subspace of dimension k  of a space of dim d.
# output: (using ihermite) a basis of k elements (of lenght d) of the projected lattice  on lin(iota).
# We will use over and over again this list H1,H2,..., Hk, so that we will work in Z^k  (embedded in R^d via H1,H2,..Hk).
# EXAMPLE:
#projectedlattice([[1,3,0],[0,1,0],[0,0,2]],[1,3])-># [[0, 1/2, 0]];
projectedlattice:=proc(W,iota) local m,B, d,k,i,r,S,IS,List,M_inverse, temp_projectedVectors;
    d:=nops(W);
    B:=ortho_basis(d);
    k:=nops(iota);
    m:=abs(det(transpose(matrix([seq(W[i],i=1..nops(W))]))));

    M_inverse:=inverse(transpose(matrix([seq(W[i],i=1..nops(W))])));
    for i from 1 to d do

        temp_projectedVectors:=m*projectedvector(M_inverse, W,iota,B[i]);
        r[i]:=[seq(temp_projectedVectors[j],j=1..nops(W))];
    od;
    S:=matrix([seq(r[i],i=1..d)]);;
    IS:=ihermite(S);
    List:=[seq(1/m*convert(row(IS,j),list),j=1..k)];
    List;
end:

# Projected cone.
# Input: W is a Cone in Z^d and iota is a subset of [1,..,d] of cardinal k;
#  Output: A "Cone" in Z^k;
# Be careful: our input must have integral coordinates.
# The ouput then will have integral coordinates.
# Here W is the cone and we are projecting W over lin( iota) and expressing it in term of the standard projectedlattice(W,iota).
#  Example: projectedconeinbasislattice([[1,1,0],[0,1,0],[0,0,2]],[1,3])→[[1,0],[0,1]]
projectedconeinbasislattice:=proc(W,iota) local P,M,output,i,F;
    P:=projectedlattice(W,iota);
    M:=transpose(matrix([seq(P[i],i=1..nops(P))]));
    output:=[];
    for i from 1 to nops(iota) do
        F:=convert(linsolve(M,Vector(W[iota[i]])),list);
        output:=[op(output),primitive_vector(F)];
    od;
    output;
end:

# Todd  function
# Toddzero(x):  the function x/(1-exp(x)));
Toddzero:=proc(x);
    x/(1-exp(x));
end:

# Relative volume
# Input: W is a Cone in R^d and iota is a subset of [1,..,d] ( of cardinal k);
# Ouput: a number;
# Math; the volume of the Box(v[i], i not in iota), with respect to the intersected lattice.
# Example: relativevolumeoffaceiotac([[1,1],[0,1]],[1])->1;
relativevolumeoffaceiotac:=proc(W,iota) local DD,iotac,P,M,H,MM,output;
    DD:=[seq(i,i=1..nops(W))];
    iotac:=GeneralComplementList(iota,DD);
    if iotac=[] then
        output:=1;
    else
        P:=matrix([seq(W[iotac[i]],i=1..nops(iotac))]);
        M:=transpose(matrix(P));
        H:=ihermite(M);
        MM:=matrix([seq(row(H,i),i=1..nops(iotac))]);
        output:=det(MM);
    fi;
    output;
end:

#relativevolumeoffaceiotac([[1,0],[0,1]],[1]);
# The 2 functions to compute S_L
# Input:   W a "Cone" in R^d; iota a subset of [1, 2,...,d];
# x a variable:
# Output: a  function of x;
# Math: #We compute integral over the cone iotac of
# exp^(csi,x) ; the answer is,  volume/ product of linear forms;
functionIzero:=proc(W,iota,x)
local DD,iotac,d,T,i,y,r,out;
    d:=nops(W);
    DD:=[seq(i,i=1..d)];
    if nops(iota)=nops(W) then
        out:=1;
    else
        iotac:=GeneralComplementList(iota,DD);
        r:=relativevolumeoffaceiotac(W,iota);
        T:=1;
        for i from 1 to nops(iotac) do
            y:=add(W[iotac[i]][j]*x[j],j=1..d);
            T:=T*y;
        od;
        T:=(-1)^(nops(iotac))*T;
        out:=r/T;
    fi;
    out;
end:

# Input:  x=[x1,x2,..,xd];  a list of symbolic expressions, W a cone in R^d.
# Output: a symbolic expression.
# Math: Our cone has generator w1,w2,...,wd.
# We replace x by <x,w_i> and we compute  the product of Todd0(<x,w_i>);
prod_Toddzero:=proc(x,W) local d,E,i,T,y;
    d:=nops(W);
    T:=1;
    for i from 1 to d do
        ASSERT(nops(W[i])=nops(x),"W[i], x need to be of the same length");
        y:=add(W[i][j]*x[j],j=1..nops(W[i]));
        T:=T*Toddzero(y);
    od;
    T;
end:

# Input:  x=[x1,x2,..,xd]  a list of symbolic expression, W a cone in R^d.
# Output: a  function of x.
# Math: P1 is the   product of Todd0(<x,w_i>), while Q1 is  the product of the (<x,wi>)
functionSzero:=proc(x,W) local P,Q,y,i;
    P:=prod_Toddzero(x,W);
    Q:=1;
    for i from 1 to nops(W) do
        ASSERT(nops(W[i])=nops(x),"W[i], x need to be of the same length");
        y:=add(W[i][j]*x[j],j=1..nops(W[i]));
        Q:=Q*y;
    od;
    P/Q;
end:

# Input: a Cone W;  iota a subset of [1..d] of cardinal k; x a list [x1,x2,...,xd]:
# Ouput: a list of  k linear forms
#  Math:
# We write R^d=V(iota)+V(iota_c). We computed a basis H1,H2n...H_k of the projection of the lattice Z^d in V(iota).
# Thus the output is the list is <x,h_i> where H_i are the basis of the projected lattice
# Example: changeofcoordinates([[1,0,0],[0,1,0],[1,2,3]],[1,2],[x1,x2,x3])->
changeofcoordinates:=proc(W,iota,x) local H,newx,i;
    H:=projectedlattice(W,iota);
    newx:=[];
    for i from 1 to nops(H) do
        newx:=[op(newx),innerprod(x,H[i])];
    od;
    newx;
end:


#####################################################
# Computation of the the function  S_iota for a cone.
# THE I of the paper is called here iota: reason I in maple means squareroot of minus  -1
# THIS IS THE MAIN TECHNICAL  PROCEDURE:

# Input: W a cone in R^d, iota a subset of [1,...,d]/
# reg a list of length d , k an integer
# Output: COMPUTE THE TERM IN  delta^(d-k) of function_S_iota(delta xi).

# Math: we will have to add cone by cone
# e^(ts,xi)S_C(xi); and compute the coefficient homogeneous of debree 0 in xi.
# and t^(d-k) in t; Thus we need
#  t^(d-k)<s,delta*xi>^(d-k) S_C(delta xi);
# Thus we need to compute the term in delta^(-d+k) of S_C(delta)(xi).
# The output is a number. We will have to multiply this number by <s,xi>^(d-k) LATER.
# The subspace $L$ where we integrate is the following face of W: L is the linear span of
# <w[j]>, with j running in THE COMPLEMENT of iota.
coeff_minusdplusk_iota_function_S:=proc(W,iota,reg,k) local d,x,i,uni_cones,function_on_iota,function_on_iotac,W_projected,WW,WWW,signuni,signL,j,iotac,
    functionconeuni, conff,out1,out2,M,newx,newP, seriesff, tt; #added tt and seriesff to local
    d:=nops(W);
    x:=[seq(delta*reg[i],i=1..d)];
    function_on_iotac:=functionIzero(W,iota,x);
    #from here express in terms of the basis lattice for projected cone.
    W_projected:=projectedconeinbasislattice(W,iota):
    if W_projected=[] then
        out1:=1
    else
        newx:=changeofcoordinates(W,iota,x);tt:=time();
        ##print('timebeforeconedec',tt);
        uni_cones:=cone_dec(W_projected): #print(nops(uni_cones),time()-tt);
        out1:=0;
        for j from 1 to nops(uni_cones) do
            WWW:=uni_cones[j][3];
            signuni:=uni_cones[j][1];
            ASSERT(abs(uni_cones[j][2])=1, "decomposition not unimodular");
            newP:=inverse(transpose(matrix(WWW)));
            function_on_iota:=functionSzero(newx,WWW);
            functionconeuni:=function_on_iotac*function_on_iota;
            ##print('functionconeuni',functionconeuni);
            seriesff:=convert(series(functionconeuni,delta=0,k+2),polynom);
            conff:=coeff(convert(series(functionconeuni,delta=0,k+2),polynom),delta,-d+k);

            out1:=out1+signuni*conff;
        od:
    fi;
    out1;
end:

#############################################################
#  Approximation "a  la Barvinok" for a simplicial cone;
# Input: W a cone in R^d,  reg a vector in R^d,k an integer;
# Output: an integer;
# Math;  We compute the coefficient in delta^(-d+k) of S(W)(delta reg), via the "special" Barvinok linear combination of
# intermediate valuations (at level k), that is summing over faces of dimension less than k and integrating over complement face;
# sum_iota cc_iota Sum(W_iota) Integral (W_iotac), where iota varies over all subsets of cardinality less or equal than k.
# We should include iota the empty set, but we dont, as we know that this corresponds (after summing over vertices) to the Lawrence calculation for the volume. So this we will give directly later.
coeff_minusdplusk_S:=proc(W,k,reg) local output,d,j,C,a,K,cc;
    d:=nops(W);
    if k=d then
        RETURN(coeff_minusdplusk_iota_function_S(W,[seq(i,i=1..d)],reg,d));
    fi;
    output:=0;
    for j from 1 to k do
        C:=choose(d,j);
        cc[j]:=(-1)^(k-j)*binomial(d-j-1,d-k-1);
        for a from 1 to nops(C) do
            output:=output+cc[j]*coeff_minusdplusk_iota_function_S(W,C[a],reg,k);
        od;
    od:
    output;
end:

#######################################################
#  Top Ehrhart coefficients for S a vertex with integral vertices.
# Input:  S a simplex, k an integer, reg  a list of lenght d.
# Output: a number (or an error message if reg is not regular).
# This is the coefficient of  t^(d-k) of the Ehrhart polynomial.
# Math: we collect all terms along the vertices of our 'aproximate' S(W_s), where
# W_s is the tangent cone at s to the simplex S.
coeff_dminusk_Eh_with_reg:=proc(S,k,reg) local M,out,F,W,i,st,d,y;
    F:=0;
    d:=nops(S)-1;
    if k=0 then
        M:=matrix([seq(S[j]-S[1],j=2..d+1)]);
        out:=abs(det(M)/d!);
    else
        for i from 1 to nops(S) do
            W:=[seq(primitive_vector(S[j]-S[i]),j=1..i-1),seq(primitive_vector(S[j]-S[i]),j=i+1..nops(S))];
            st:=[seq(S[i][u],u=1..d)];
            y:=add(st[j]*reg[j],j=1..d); ##print('i,W,k,reg',coeff_minusdplusk_S(W,k,reg));
            F:=F+1/(d-k)!*y^(d-k)*coeff_minusdplusk_S(W,k,reg);
        od:
        out:=F;
    fi;
    out;
end:


###################################################
# We require a random vector; in principle we could compute
# (in polynomial time) a regular vector (i.e, not orthogonal to any edge)
# deterministically, but the deterministic calculation of a regular
# vectors would take more time and in practice a random choice works.
#
# Note that you can get an error message
# if the random vector chosen in the procedure below is not regular.
# This happens with low probability, but if it does RERUN.
#
# INPUT: N a large integer; d  an integer.
# OUTPUT : a random list of lenght d with coefficients between 1 and N.
random_vector:=proc(N,d) local R;
    R:=rand(N);
    [seq(R()+1,i=1..d)]:
end:


##########################################################
# THIS IS ONE OF THE First of TWO MAIN PROCEDURES:
#
# INPUT: A SIMPLEX full dimensional in R^d; k an integer;
# OUTPUT:  the coefficient of degree t^(d-k) of the Ehrhart polynomial (unweighted) of S;
# or an error message if the random vector chosen in the procedure is not regular.
# In this case RERUN.
coeff_dminusk_Eh:=proc(S,k) local d,reg;
    d:=nops(S)-1;
    reg:=random_vector(5000,d);
    coeff_dminusk_Eh_with_reg(S,k,reg);
end:


##ATTENTION
#-------------------------BRANDON/GREG MAKE SURE THESE TWO SIMPLE
# EXAMPLES ARE TESTED ALWAYS FIRST, it make sense because we know the
# coefficients of these examples

# SOME EASY EXAMPLES.
#
Sstandard:=proc(n) local ze, S,j,zej; ze:=[seq(0,i=1..n)];
    S:=[ze];
    for j from 1 to n do zej:=subsop(j=1,ze);;
        S:=[op(S),zej];
    od;
end:

#Sstandard(3);
#S123:=[[0, 0, 0], [1, 0, 0], [1, 2, 0], [1, 2, 3]];
Srandom:=proc(d,N) local S,i;
    S:=[];
    for i from 1 to d+1 do
        S:=[op(S),random_vector(N,d)];
    od;
    S;
end:



Checkrandom:=proc(d) local S,reg,CC,k,tk;
    S:=Srandom(d,100); print(S);
    CC:=[];
    for k from 0 to d do
        tk:=time(); print(coeff_dminusk_Eh(S,k),time()-tk);
        CC:=[op(CC),[k,coeff_dminusk_Eh(S,k)] ];
    od;
    CC;
end:

#Checkrandom(4);

Sou:=proc(n) local ze, S,j,zej;
    ze:=[seq(0,i=1..n)];
    S:=[ze];
    for j from 1 to n do zej:=subsop(j=j,ze);
        S:=[op(S),zej];
    od;
end:

CheckSou:=proc(d) local S,k,tk,CC;
    CC:=[];S:=Sou(d);
    for k from 0 to d do
        tk:=time();
        print(coeff_dminusk_Eh(S,k),time()-tk);
        CC:=[op(CC),[k,coeff_dminusk_Eh(S,k)] ];
    od;
    CC;
end:

# EXAMPLE: when one runs
#CheckSou(5) it should give the result:
#1+(16/3)*t+(23/4)*t^4+t^5+(73/6)*t^3+(47/4)*t^2



##ATTENTION This is new code created with the help of Vergne et al.
################NEW CODE starts here it is faster!!


#########################################################
#
#
#
#################################################################

Topk_function_S_Ispace:=proc(W,iota,reg,k) local d,x,i,uni_cones,function_on_iota,function_on_iotac,W_projected,WW,WWW,signuni,signL,j,iotac,
    functionconeuni, conff,out1,out2,M,newx,newP, seriesff, tt; #added tt and seriesff to local
    d:=nops(W);
    x:=[seq(delta*reg[i],i=1..d)];
    function_on_iotac:=functionIzero(W,iota,x);
    #from here express in terms of the basis lattice for projected cone.
    W_projected:=projectedconeinbasislattice(W,iota):
    if W_projected=[] then
        out1:=1
    else
        newx:=changeofcoordinates(W,iota,x);tt:=time();
        uni_cones:=cone_dec(W_projected): #print(nops(uni_cones),time()-tt);
        out1:=0;
        for j from 1 to nops(uni_cones) do
            WWW:=uni_cones[j][3];
            signuni:=uni_cones[j][1];
            ASSERT(abs(uni_cones[j][2])=1, "decomposition not unimodular");
            newP:=inverse(transpose(matrix(WWW)));
            function_on_iota:=functionSzero(newx,WWW);
            functionconeuni:=function_on_iotac*function_on_iota;
            seriesff:=convert(series(functionconeuni,delta=0,k+2),polynom);
            #conff:=coeff(convert(series(functionconeuni,delta=0,k+2),polynom),delta,-d+k);

            out1:=out1+signuni*seriesff;
        od:
    fi; out1;
end:

#######################
# Input is a cone and value k
# It expresses the valuation at that cone as
# a linear combination (poset Moebius) of simpler
# valuations using the projections to low dimension spaces.
#
##########################
Topk_S_Cone:=proc(W,k,reg) local output,d,j,C,a,K,cc;
    d:=nops(W);
    if k=d then
        RETURN(Topk_function_S_Ispace(W,[seq(i,i=1..d)],reg,d));
    fi;
    output:=0;
    for j from 1 to k do
        C:=choose(d,j);
        cc[j]:=(-1)^(k-j)*binomial(d-j-1,d-k-1);
        for a from 1 to nops(C) do
            output:=output+cc[j]*Topk_function_S_Ispace(W,C[a],reg,k);
        od;
    od:
    output;
end:


#####################################################
# This is the second MAIN PROCEDURE.

# The code Topk_Eh outputs a POLYNOMIAL:
# If Eh(t)=\sum_{i=0}^d e_i t^(d-i);
# the ouput is sum_(i=0}^k e_i t^(d-i);
# Top Ehrhart coefficients for S a vertex with integral vertices.
# are precisely the few coefficients of the polynomial.
#
# Input;  S a simplex, k an integer, reg  a list of lenght d.
# Output: a polynomial in t.
#
# This function relies on Brion's theorem that allows to
# decompose the expression as a sum over each tangent cone at vertices
# of the simplex. Then each cone will be treated separately (call to Topk_S_Cone)
#####################################
Topk_Eh:=proc(S,k,t) local M,out,F,W,i,st,d,y,reg,AA;
    d:=nops(S)-1;
    reg:=random_vector(5000,d);
    M:=matrix([seq(S[j]-S[1],j=2..d+1)]);
    F:=t^d/d!*abs(det(M));
    for i from 1 to nops(S) do
        W:=[seq(primitive_vector(S[j]-S[i]),j=1..i-1),seq(primitive_vector(S[j]-S[i]),j=i+1..nops(S))];
        st:=[seq(S[i][u],u=1..d)];
        y:=add(st[j]*reg[j],j=1..d);
        F:=F+coeff((add(1/(d-i)!*y^(d-i)*delta^(d-i)*t^(d-i),i=1..k)*Topk_S_Cone(W,k,reg)),delta,0):
    od:
    F;
end:

##ATTENTION EXAMPLES OF HOW TO RUN
# UNCOMMENT
##############

# Either done coefficient by coefficient

#told6:=time(): coeff_dminusk_Eh(Sou(6),0); coeff_dminusk_Eh(Sou(6),1); coeff_dminusk_Eh(Sou(6),2); coeff_dminusk_Eh(Sou(6),3); time()-told6;

# Or done 3 coefficients in one shot!

#tnew6:=time(): Topk_Eh(Sou(6),3,t);time()-tnew6;



