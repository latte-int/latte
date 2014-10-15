with(linalg):with(LinearAlgebra):with(combinat):
kernelopts(assertlevel=1):       ### Enable checking ASSERTions

#############################################################
bonfrac:=proc(t); frac(t-floor(t));end:  
fractionalpart:=proc(s) local our,T;
    if t=0 then our:=0; else our:=ourfrac(t);
    our;fi;end:

 
  
ourmod:=proc(p,q,t) local our,T;
    if q=1 or modp(p,q)=0 then our:=0;
    elif type(t,integer) then our:=modp(t*p,q);
    else our:=q*ourfrac(modp(p,q)*t/q);
    fi;
    RETURN(our); end:

nfractionalpart:=proc(n,p,q) local our;
    if  type(n,rational) then our:=fractionalpart(p*n/q)
    else our:=1/q*ourmod(p,q,n);
    fi; our; end:


ourmodreal:=proc(p,q,t) local our,T;
    if type(t,integer) then our:=modp(t*p,q);fi;
    if t=0 or p=0 then our:=0;
    else our:=q*ourfrac(p*t/q);
    fi;
    our;
end:
ourmodreal(1,2,-t);
nfractionalpartreal:=proc(n,p,q) local our;
    if  type(n,rational) then our:=fractionalpart(p*n/q)
    else our:=1/q*ourmodreal(p,q,n);
    fi;
    our;
end:


##############################################################

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
special_lincomb_v:=proc(a,v,n) local out;
    ASSERT(nops(a)=nops(v)," the number of coefficients and vectors do not match");
    if v=[]   then out:=[seq(0,i=1..n)];else
        out:=[seq(add(a[i]*v[i][j],i=1..nops(v)),j=1..nops(v[1]))];
    fi;
    out;
end:
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
ortho_basis:=proc(d) local i,v;
    for i from 1 to d do
        v[i]:=[seq(0,j=1..i-1),1,seq(0,j=i+1..d)]
    od;
    [seq(v[j],j=1..d)];
end:
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

good_vector:=proc(G) local n,A,Ainverse,B,sho,V,L;
    n:=nops(G);
    Ainverse:=inverse(matrix(G));
    B:=convert(Ainverse, listlist);
    sho:=short_vector(B);
    V :=[seq(add(G[j][i]*sho[j],j=1..n),i=1..n)];
    L:= sign_entries_vector(sho);
    [V,L];
end:

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

good_cone_dec:=proc(eps,G) local n,A,R,Output, det_A;
    n:=nops(G);  A:=matrix([seq(G[i],i=1..n)]);
    det_A:=det(A);
    if abs(det_A)=1 then
        Output:=[[],[[eps,det_A,G]]];
    else R:=good_vector(G);
        Output:=signed_decomp(eps,G,R[1],R[2]);
    fi;
end:

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
projectedvector:=proc(W,Cspace,b) local M,S,j,v,V,m;
    M:=transpose(matrix([seq(W[i],i=1..nops(W))]));
    S:=linsolve(M,b);
    m:=det(M);
    for j from 1 to nops(W) do
        v[j]:=add(S[Cspace[i]]*W[Cspace[i]][j],i=1..nops(Cspace));
    od:
    V:=[seq(v[j],j=1..nops(W))];
end:
projectedvector_with_inverse:=proc(M_inverse, W,Cspace,b) local S,j,v,V;
    S:=multiply(M_inverse,b);
    for j from 1 to nops(W) do
        v[j]:=add(S[Cspace[i]]*W[Cspace[i]][j],i=1..nops(Cspace));
    od:
    V:=[seq(v[j],j=1..nops(W))];
end:
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

projectedvertexinbasislattice:=proc(W,Cspace,ProjLattice,s) local m,P,M,output,i,F;
    P:=ProjLattice;
    if Cspace=[] then RETURN([]);fi;
    M:=Transpose(Matrix([seq(P[i],i=1..nops(P))]));
    F:=convert(LinearSolve(M,Vector(projectedvector(W,Cspace,s))),list);
    output:=F;
end:
s_ISpace:=proc(s,W,ISpace) local M,s_in_cone_coord,s_ISpace;
    M:=Matrix([seq(Vector([W[i]]),i=1..nops(W))]);
    s_in_cone_coord:=convert(LinearSolve(M,Vector(s)),list);
    s_ISpace:=[seq(s_in_cone_coord[ISpace[k]],k=1..nops(ISpace))];
    special_lincomb_v(s_ISpace,[seq(W[ISpace[k]],k=1..nops(ISpace))],nops(W));
end:
Todd:=proc(z,t);
    exp(z*t)*t/(1-exp(t));
end:
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
prod_Todd:=proc(z,W,xi) local d,E,i,T,y;
    d:=nops(W);
    T:=1;
    for i from 1 to d do
        y:=add(W[i][j]*xi[j],j=1..nops(W[i]));
        T:=T*TODD(z[i],y);
    od;
    T;
end:
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

changeofcoordinates:=proc(W,Cspace,ProjLattice,xi) local H,newxi,i,d;
    H:=ProjLattice;
    d:=nops(W[1]);
    newxi:=[];
    for i from 1 to nops(H) do
        newxi:=[op(newxi),add(xi[j]*H[i][j],j=1..d)];
    od;
    newxi;
end:
#changeofcoordinates([[1,0,0],[1,1,2],[0,5,1]],[1,2],...,xi);

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
# EXAMPLE: S_Ispace_Coneformulaa([s1,s2],[[1,0],[1,2]],[1],xi)->  -TODD(ceil(s2), (1/2)*xi[1]+xi[2])*EXP((s1-(1/2)*s2)*#xi[1])/(((1/2)*xi[1]+xi[2])*xi[1]);

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
            s_prime_Cspace:=[seq(ceil(news[f]),f=1..nops(news))];
            function_on_Cspace:=functionS(s_prime_Cspace,WWW,newxi);
            out1:=out1+signuni*function_on_Cspace[1]/function_on_Cspace[2]*function_on_ISpace[1]/function_on_ISpace[2];
        od:
    fi;
    out1;
end:


# Input: s a vector in Q^d,  or a symbolic variable ; BUT THEN IT HAS TO BE ENTERED AS A LIST OF  d SYMBOLIC VARIABLES
# s:=[s1,s2,...,sd]; W a cone in Z^d, Ispace a subset of [1,...,d]
# xi a list of lenght d  of variable (or xi);

# The output is a function of xi[i].
#
# The subspace $L$ where we integrate is the following face of W: L is the linear span of
# <w[j]>, with j running of Ispace. (thus Ispace should be "big")
#
# Here we take out a function of s, fractionalpart
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






#ApproxEhrhartSimplexgeneric(n,[[0,0],[1/2,0],[0,1/2]], 1,xi);
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
    	fprintf(fPtr, "epoly:= %a;\\n", ehrhartPoly);
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
			od; 			
partialF:=eval(subs({TODD=Todd,EXP=exp},output)):    partialSeries[iLF, currentDifference+1]:=coeff(series(partialF,t=0,M_current+d+2),t,M_current):			
partialSeries[iLF, currentDifference+1]:=coeff(series(partialSeries[iLF, currentDifference+1],epsilon=0,d+2),epsilon,0):			
od; 
        
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


#####################################################################
### LATTE INTERFACE HELPER FUNCTIONS:
####################################################################

#### LATTE INTERFACE FUNCTION:
# input:
#	a: any rational number
#	n: base
#	return the number between [0,n) that is equal to  a mod n
latteMod:=proc(a, n)
	local r, x;
	ASSERT(n > 0);
	x:=a;
	while ( x >= n or x < 0) do
		x:= x - floor(x/n)*n;
	end;
	return x;
end:


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

    F:=eval(subs({TODD=Todd,EXP=exp},F)):    
    return F;
end:

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
thefirstk0Ehrhartweighted:=proc(n,Simplex,ell,M,k0) local d;
    d:=nops(Simplex)-1;
    add(TopEhrhartweighted(n,Simplex,ell,M,m)*n^m,m=M+d-k0..M+d);
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
    F:=eval(subs({TODD=Todd,EXP=exp},F)):
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
    reg:=random_vector(1000,d);
    xx:=[seq(t*(ell[i]+epsilon*reg[i]),i=1..d)];
    AA:=ApproxEhrhartSimplexgeneric_real(n,Simplex,order,xx);
    CCt:=coeff(series(AA,t=0,M+d+2),t,M); #print(CCt);
    CCeps:=coeff(series(CCt,epsilon=0,d+2),epsilon,0);
    CCn:=coeff(CCeps,n,m);
    subs({N=n},CCn);
end:
CompleteEhrhartweighted_real:=proc(n,nn,Simplex,ell,M) local d;
    d:=nops(Simplex)-1;
    add(TopEhrhartweighted_real(n,Simplex,ell,M,m)*nn^m,m=0..M+d);
end:





cone_by_cone:=proc(Simplex,ell,M,order) local reg,d,xx,AA,CCt,CCeps,CCn;
d:=nops(Simplex)-1;CCn:=0;
#order:=M+nops(Simplex)-codim; 
reg:=random_vector(5000,d); 
xx:=[seq(t*(ell[i]+epsilon*reg[i]),i=1..d)];
AA:=ApproxEhrhartSimplexgeneric(n,Simplex,order,xx);
CCt:=coeff(series(AA,t=0,M+d+2),t,M); 
CCeps:=coeff(series(CCt,epsilon=0,d+2),epsilon,0);
CCn:=add(coeff(CCeps,n,m)*t^(m),m=0..M+d);
subs({N=n},CCn);
end:

cone_by_cone_real:=proc(Simplex,ell,M,order) local reg,d,xx,AA,CCt,CCeps,CCn,newCCn;
d:=nops(Simplex)-1;CCn:=0;
reg:=random_vector(5000,d); 
xx:=[seq(t*(ell[i]+epsilon*reg[i]),i=1..d)];
AA:=ApproxEhrhartSimplexgeneric_real(n,Simplex,order,xx);
CCt:=coeff(series(AA,t=0,M+d+2),t,M); 
CCeps:=coeff(series(CCt,epsilon=0,d+2),epsilon,0);
CCn:=add(coeff(CCeps,n,m)*t^(m),m=0..M+d);
newCCn:=subs({N=n},CCn); newCCn:=subs(n=t,newCCn);
subs(ourfrac=bonfrac,newCCn);
end:
testfrac := [[1, 1, 1], [4, 2, 1], [1, 1, 2], [2, 2, 2]]; testdim2:=[[1, 1], [1, 0], [0, 1]];testdim4:=[[6,9,5,1],[3,5,4,0],[7,4,9,1],[1,3,7,2],[8,9,1,7]];
conetestfrac0:=cone_by_cone_real(testfrac,[1,5,7],0,0);
conetestfrac1:=cone_by_cone_real(testfrac,[1,5,7],0,1);
conetestfrac2:=cone_by_cone_real(testfrac,[1,5,7],0,2);
conetestfrac3:=cone_by_cone_real(testfrac,[1,5,7],0,3);
#conetestdim2:=cone_by_cone_real(testdim2,[1,5],0,2);
#conetestdim4:=cone_by_cone_real(testdim4,[1,1,1,1],0,4);

#plot(conetestdim4,t=0..1);
plot([conetestfrac0,conetestfrac1,conetestfrac2,conetestfrac3],t=1..2,color=[green,blue,red,black], axis[1]=[gridlines=[20, thickness=1, subticks=false]],axis[2]=[gridlines=[7, thickness=1, subticks=false]],thickness=2);
seq(eval(subs(t=i/24, conetestfrac3)),i=0..50);
#conetestfrac3,  conetestfrac2,conetestfrac1,conetestfrac0
#testdim2:=[[2, 7], [6, 8], [9, 0]];
#SS:=subs(t=2/5+eps,conetestfrac3);

#seq(eval(subs(t=1+i/10,conetestdim2)),i=-10..10);


