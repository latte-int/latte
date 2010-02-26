with(linalg):with(LinearAlgebra):
with(numapprox,laurent):
read("integration/createLinear.mpl"):
# Integral of a  power of a linear form over a simplex.
# Our Notations;

# 
# 

# The  integer d is the dimension;
# A vector in Q^d is a list of d rational numbers.
# A vertex is a vector in Q^d;
# The simplex  S is the convex hull of its vertices s_i. 
# Thus S is encoded as a list of vectors
# in Q^d.
# If the simplex is of full dimension, we have (d+1) vertices.  

#  A linear forms is called alpha, it is  represented by  a vector in Q^d.
# A monomial m is a list of d integers
#  A polynomial represented  in a sparse way;
#  Exemple x*y^2+2*x^2 with be given as a list of lists [[1,[1,2]],[2,[2,0]]. Each list represents a monomial with his coefficients. 
# Thus a sparse polynomial is represented as a list of lists.
# 
# .
# 
# Simplex and multiplicities.
# 
# INPUT: d an integer, S  list of d+1 lists of lenght d, alpha: list of lenght d . 
# OUTPUT: set of lists {[a_S], [m_S]} 
# MATH: S a simplex of dimension d+1, alpha a linear form, m_S is the list of the number of vertices S where <\alpha,S> = a_S.
# 
#    
multiplicity_alpha_simplex:=proc(S,d,alpha) local i,n,VS,m,Mult,j,c;
n:=nops(S);
VS:={seq(add(alpha[s]*S[i][s],s=1..d),i=1..nops(S))};
Mult:={};
for i from 1 to nops(VS) do 
m:=0; 
for j from 1 to nops(S) do 
c[j]:=add(alpha[s]*S[j][s],s=1..d);
if c[j]=VS[i] then m:=m+1;
else m:=m;
fi;
od;
Mult:={op(Mult),[VS[i],m]};
od;
Mult:
end:
# Computation of a coefficient
# 
# 
# 
# 
# M an integer, d an integer,  starting is a list of two elements  [v1,order1], where v1 is a rational number and order1 is an integer ;
#      L is a sequence of elements [ai,mi] where mi are integers and ai are numbers.  The output is 
# coeff((epsilon+v1)^(M+d)*1
# /(epsilon+a1)^m1*1/(epsilon+a2)^m2*...*1/(epsilon +aK)^mK; epsilon,order1-1)
prototype_residue:=proc(M,d,starting,L) local f,m,LL;
f:=(epsilon+starting[1])^(M+d); #print(f);
for m from 1 to nops(L) do
f:=f*1/(L[m][1])^(L[m][2])*1/(1+epsilon/L[m][1])^(L[m][2]);
od; 
f;
LL:=convert(laurent(f,epsilon,starting[2]),polynom);
coeff(LL,epsilon,starting[2]-1);
end:
# INPUT: S a simplex, alpha a linear form, d an integer, M an integer 
# OUTPUT: a number int_S alpha^M
# MATH:;
# See the manual for the formula.  
# 
integral_power_linear_form:=proc(S,d,M,alpha) local int,Mult,output,v,i,L,B,R,m,starting;
v:=abs(Determinant(Matrix([seq(S[j]-S[1],j=2..d+1)])));
if v=0 then int:=0;
else 
 Mult:=multiplicity_alpha_simplex(S,d,alpha); 
int:=0;
for  i from 1 to nops(Mult) do 
starting:=Mult[i];
 L:=
[seq([Mult[i][1]-Mult[j][1],Mult[j][2]],j=1..(i-1)),
seq([Mult[i][1]-Mult[j][1],Mult[j][2]],j=(i+1)..nops(Mult))];
R:=prototype_residue(M,d,starting,L);
int:=int+R;   
od;
fi;
M!/(M+d)!*int*v:
end:

# Decomposing  a monomial in powers of linear forms.
# INPUT:L list of integers, m a list of integers,M an integer
# OUTPUT: a nonnegative number
# MATH: The list $m$ represents the monomial x^m=x_1^{m[1]}x_2^{m[2]}\cdots x_d^{m[d], 
# the list L represents the linear form L=L[1]x[1]+..L[d]x_d. 
# The output is the coefficient of L^M,  M=M:=add(m[i],i=1..nops(m)), in the expansion of x^m in linear form using the formula in userguide. 
# CAUTION: The algorithm works only if we take "primitive" linear form as in the formula. 
coeff_linear_form_expansion:=proc(L,m) local p,j,c,s,M,out;
p:=[];
M:=add(m[i],i=1..nops(m));
    if igcd(seq(L[i],i=1..nops(L)))<>1 then print(trouble); 
    else
        for j from 1 to nops(L) do #print("j",j);
           if L[j]<>0 then
           p:=[op(p),iquo(m[j],L[j])];
           else p:=p;
           fi;
         od;

    c:=min(seq(p[i],i=1..nops(p)));
    s:=add(L[i],i=1..nops(L));#print("s,m,c",s,m,c);
    out:=1/M!*add((-1)^(M-k*s)*product(binomial(m[i],k*L[i]),i=1..nops(L))*k^M,k=1..c);
    fi;
out;
end:
# The input is a monomial.
# The output is the list of  elements [p1,....,pn] where$pi<=mi$ and mutually prime.
list_for_simplex_integral:=proc(m) local j,F,L,i,f,newL;newL:=[];
i:=1; 
if m[1]=0 then L:=[[0]]; else
L:=[seq([j],j=0..m[1])];fi;

     for i from 2 to nops(m) do #print(i,nops(L));
     L:=[seq(seq([op(L[s]),j],j=0..m[i]),s=1..nops(L))];
     od;
#print(L);
     for j from 1 to nops(L) 
     do F:=L[j];

       if
       igcd(seq(F[i],i=1..nops(F)))<>1
       then newL:=newL; 
       else newL:=[op(newL),F];
       fi;
     newL:
     od;
newL;
end:
# INPUT: m a list of integers, coe a number
# OUTPUT: a list of lists of lenght nops(m)
# MATH:  The list $m$ represents the monomial x^m=x_1^{m[1]}x_2^{m[2]}\cdots x_d^{m[d].
# The output is a list of lists. Each element in the list represents a linear form ([1,2]=x+2y). The output exausts all the linear form with exponents M=m[1]+..+m[d] which appear when  expressing   x^m as  linear combinations of linear form with exponent M. The first element is the coefficient multiplied by coe;
list_and_coeff_for_monome:=proc(m,coe) local M,L,out:
M:=add(m[i],i=1..nops(m));
L:=list_for_simplex_integral(m);#print(L);
out:=[ seq([coe*coeff_linear_form_expansion(L[j],m),[M,L[j]]],j=1..nops(L))];
end:
#  The input: S is  a simplex, d a number, m a monomial. The output is  a number, the integral of x^m over S.
integral_monome_via_waring:=proc(S,d,m) local out,M,L,i;
out:=0;
M:=add(m[i],i=1..nops(m));
L:=list_and_coeff_for_monome(m,1);
for i from 1 to nops(L) do
   out:=out
+L[i][1]*integral_power_linear_form(S,d,M,L[i][2][2]);
   od;
out;
end:
##integral_monome_via_waring([[0,0],[0,1],[1,0]],2,[9,2]);
# Integral of a polynomial via Waring.
# We give a polynomial in a sparse way; Exemple x*y^2+2*x^2 with be given as a list of lists [[1,[1,2]],[2,[2,0]]. Each list represents a monomial with his coefficients. 
#  We start by cleaning the sets for example we replace [[1,[1,2]],[1,[1,2]] by [2,[1,2]];
# Input  L; a list of lists [[a,alpha],[b,beta],...]. here a is a number, alpha is a list. The ouptput is
# a set of lists.
# If alpha=beta, we replace by [a+b,alpha]; if a=0 we skip;
# The input is a list  L of lists [a,\alpha] where a is a number. The output is of the same kind.
cleaned_set:=proc(L) local newL,subL,X,i;
newL:=[]; 
for i from 1 to nops(L) do 
if L[i][1]<>0 then
newL:=[op(newL),L[i]];
fi;
od;
subL:={seq(newL[s][2],s=1..nops(newL))};
X:=add(newL[s][1]*x[newL[s][2]],s=1..nops(newL));
{seq([coeff(X,x[subL[i]],1),subL[i]],i=1..nops(subL))};
end:
# The input is  a sparse polynomial represented as a lists [c,m] of monomials with coefficients. The output is a list of lists [coe,[M,L]], coe a number, M an integer, L a linear form. It represents coe*L^M. 
list_integral_via_waring:=proc(sparse_poly) 
local Y, new_sparse_poly,n,i,Z;
new_sparse_poly:=cleaned_set(sparse_poly);
n:=nops(new_sparse_poly);
Y:=list_and_coeff_for_monome(new_sparse_poly[1][2],new_sparse_poly[1][1]); 
 for i from 2 to n do 
Z:=list_and_coeff_for_monome(new_sparse_poly[i][2],new_sparse_poly[i][1]); 
Y:=cleaned_set([op(Y),op(Z)]);
od; 
Y;
end:
# The input is a simplex S, d the dimension, sparse_poly a sparse polynomial. 
# The ouput is a number; the integral over S of the polynomial.
integral_via_waring:=proc(S,d,sparse_poly) local output,i, L ;output:=0;
L:=list_integral_via_waring(sparse_poly);
for i from 1 to nops(L) do 
output:=output+L[i][1]*
integral_power_linear_form(S,d,L[i][2][1],L[i][2][2]);
   od;
end:

#

with(combinat):

lattice_random_simplex:=proc(d,N) local R,U;
  R := rand(N):
  U:=proc()[seq(R(),i=1..d)] end proc:
  [ seq(U(), i=1..d+1) ];
end:

## Computing a uniformly random monomial
## of prescribed degree.

## The number of compositions of M into n non-negative integers can be
## computed by the Maple function numcomp(M+n, n).
## There is the recursion numbcomp(M+n, n) = sum(numbcomp(M+n-1-i, n-1), i
## = 0...M), obtained by running the first number (m_1) from 0 to M.
## 
## So you can choose the first number m_1 between 0 and M with
## probabilities given by how many random choices there are. Then recurse
## to choose the remaining numbers.

### This is very very slow. --Matthias

## partial_sums := proc(degree, dimension)
##   #option remember;
##   counts := [ seq(numbcomp(degree+dimension-1-i, dimension-1),
##                     i = 0..degree) ];
##     #print(["counts", counts]);
##     total_count := add(c, c in counts);
##     #print(["total_count", total_count]);
##     pos := rand(total_count)();
##     #print(["pos", pos]);
##     partial := 0;
##     partials := [];
##     choice := -1;
##     for i from 0 to degree do
##       #print(i);
##       partial := partial + counts[i+1];
##       partials := [ op(partials), partial ];
##       #print(partial);
##     od;
##   partials;
## end:

## random_monomial_of_given_degree := proc(dimension, degree)
##   local x;
##   #print("random_monomial_of_given_degree", dimension, degree);
##   if dimension = 1 then
##     [degree]
##   elif dimension = 2 then
##     x := rand(degree + 1)();
##     [x, degree - x]
##   else
##     partials := partial_sums(degree, dimension);
##     total_count := partials[degree+1];
##     #print(["total_count", total_count]);
##     pos := rand(total_count)();
##     #print(["pos", pos]);
##     choice := -1;
##     for i from 0 to degree do
##       #print(i);
##       partial := partials[i+1];
##       #print(partial);
##       if pos < partial then
##         choice := i;
##         break;
##       end if;
##     od;
##     assert(choice >= 0);
##     [ choice, 
##       op(random_monomial_of_given_degree(dimension - 1, degree - choice)) ]
##   end if
## end:


## Converting from Maple polynomials to our sparse format

polynomial_to_sparsepoly := proc(p, dimension)
  local variables,coefficients, monomials, exponents; 
  variables := [ seq(x[i], i=1..dimension) ];     
  #print (variables);     
  coefficients := coeffs(p, variables, 'monomials');
  #print (coefficients);
  #print (monomials);
  exponents := map( monomial -> map( variable -> degree(monomial, variable),
                            variables),
                     [monomials]);
  #print(exponents);
  zip((c,m) -> [c,m],
      [coefficients], 
      exponents);
end:

#polynomial_to_sparsepoly(x[1]*x[2]^2+2*x[2]-3*x[1], 2);

# We rather use a built-in procedure of Maple,
# and convert the result to our format.

##
## Creating random polynomials in sparse format
##

# N - use coefficients up to this number
# d - dimension
# M - degree
# r - number of monomials
## random_sparse_homogeneous_polynomial_with_degree:=proc(N,d,M,r) 
##   local MM, out,R,U,A,m,coe,k,mono;out:=[];
##   out:=[];
##   R := rand(N):
##   U:=proc()[seq(R(),i=1..r)] end proc:
##   A:= U();#print(U);
##   for k from 1 to r do
##   mono := random_monomial_of_given_degree(d,M);
##   out:=[op(out),[A[k],mono]];
##   od;
##   out;
## end:

random_sparse_homogeneous_polynomial_with_degree:=proc(N,d,M,r) 
  local p, R;
  ## Give up if too large polynomials requested
  if (r > 500000) then
    error "Too large a polynomial requested"
  fi;
  R := rand(N);
  p := randpoly([ seq(x[i], i=1..d) ], 
                homogeneous, degree = M, terms = r, coeffs = proc() R() + 1; end);
  polynomial_to_sparsepoly(p, d);
end:

##random_sparse_homogeneous_polynomial_with_degree(100000, 50, 1000, 1);


# One monomial
random_sparse_homogeneous_polynomial_with_degree_1:=proc(d,M)
  random_sparse_homogeneous_polynomial_with_degree(100, d, M, 1);
end:

# Very sparse
random_sparse_homogeneous_polynomial_with_degree_2:=proc(d,M)
  random_sparse_homogeneous_polynomial_with_degree(100, d, M, d+M);
end:


random_sparse_homogeneous_polynomial_with_degree_3:=proc(d,M)
  random_sparse_homogeneous_polynomial_with_degree(100, d, M, 
    floor(evalf(root(binomial(M+d-1, d-1), 4))));
end:

random_sparse_homogeneous_polynomial_with_degree_4:=proc(d,M)
  random_sparse_homogeneous_polynomial_with_degree(100, d, M, 
    floor(evalf(root(binomial(M+d-1, d-1), 2))));
end:

random_sparse_homogeneous_polynomial_with_degree_5:=proc(d,M)
  random_sparse_homogeneous_polynomial_with_degree(100, d, M, 
    floor(evalf(root(binomial(M+d-1, d-1), 1))));
end:



random_sparse_homogeneous_polynomial_with_degree_and_eff_num_vars:=
proc(N,d,eff_num_vars,M,r) 
  local p, R;
  ## Give up if too large polynomials requested
  if (r > 500000) then
    error "Too large a polynomial requested"
  fi;
  R := rand(N);
  p := randpoly([ seq(x[i], i=1..eff_num_vars) ], 
                homogeneous, degree = M, terms = r, coeffs = proc() R() + 1; end);
  polynomial_to_sparsepoly(p, d);
end:

# One monomial, constant (2/4/6) number of effective variables!
random_sparse_homogeneous_polynomial_with_degree_6:=proc(d,M)
  random_sparse_homogeneous_polynomial_with_degree_and_eff_num_vars(100, d, 
                                                        min(2, d),
                                                   M, 1);
end:

random_sparse_homogeneous_polynomial_with_degree_7:=proc(d,M)
  random_sparse_homogeneous_polynomial_with_degree_and_eff_num_vars(100, d, 
                                                        min(4, d),
                                                   M, 1);
end:

random_sparse_homogeneous_polynomial_with_degree_8:=proc(d,M)
  random_sparse_homogeneous_polynomial_with_degree_and_eff_num_vars(100, d, 
                                                        min(6, d),
                                                   M, 1);
end:

test_integration:=proc(polyCount, bigConstant, numTerms, dimension, myDegree, decomposing)
  local errors, wrong:
  local myMonomials, mySimplices, myLinForms, mapleLinForms, myResults, mapleResults:
  local curForms, curTerm, curSet:
  local myIndex, formIndex, i, j:
  local myTime, temp, intTime, L:
  local inputFile, outputFile, errorFile:
  
  #get polynomials
  myTime:=0:
  intTime:=0:
  for myIndex from 1 to polyCount do
    mySimplices[myIndex]:=lattice_random_simplex(dimension, bigConstant);
    if decomposing = 1 then
      myMonomials[myIndex]:=random_sparse_homogeneous_polynomial_with_degree(bigConstant, dimension, myDegree, numTerms);
      temp:=time():
      mapleLinForms[myIndex]:=list_integral_via_waring(myMonomials[myIndex]):
      myTime:=myTime + time() - temp:
      temp:=time():
    else
      #print(myDegree, dimension, bigConstant, bigConstant, numTerms);
      mapleLinForms[myIndex]:=random_linearform_given_degree_dimension_maxcoef_componentmax_maxterm(myDegree, dimension, bigConstant, bigConstant, numTerms):
      print(random_linearform_given_degree_dimension_maxcoef_componentmax_maxterm(myDegree, dimension, bigConstant, bigConstant, numTerms));
    end if:
    mapleResults[myIndex]:=0:
    for formIndex from 1 to nops(mapleLinForms[myIndex]) do
      mapleResults[myIndex]:=mapleResults[myIndex]+mapleLinForms[myIndex][formIndex][1]*integral_power_linear_form(mySimplices[myIndex],dimension,mapleLinForms[myIndex][formIndex][2][1],mapleLinForms[myIndex][formIndex][2][2]):
    od:
    intTime:=intTime + time() - temp:
  od:
  myTime:=myTime / polyCount:
  intTime:=intTime / polyCount:
  
  #write to file
  if decomposing = 1 then
    inputFile:=fopen("integration/randomPolys.txt",WRITE,TEXT):
    for i from 1 to polyCount do
      writeline(inputFile, convert(myMonomials[i], string));
      writeline(inputFile, StringTools[DeleteSpace](convert(mySimplices[i], string)));
    od:
    close(inputFile):
    
    #run the integrate program
    system("./integrate_test integration/randomPolys.txt integration/forms.txt 1"):
    
    print(StringTools[Join]([convert(myTime,string),"s. avg. spent on Maple decomposition."], " "));
    
    #read forms in maple notation
    outputFile:=fopen("integration/forms.txt",READ,TEXT):
    myLinForms[1]:=readline(outputFile):
    myResults[1]:=readline(outputFile):
    i:=1:
    while (myLinForms[i] <> 0) do
      i:=i+1:
      myLinForms[i]:=readline(outputFile): #decomposition into linear forms
      myResults[i]:=readline(outputFile): #integral result
    od:
    close(outputFile):
  else
    inputFile:=fopen("integration/randomForms.txt",WRITE,TEXT):
    for i from 1 to polyCount do
      writeline(inputFile, convert(myLinForms[i], string));
      writeline(inputFile, StringTools[DeleteSpace](convert(mySimplices[i], string)));
    od:
    close(inputFile):
    
    #run the integrate program
    system("./integrate_test integration/randomForms.txt integration/results.txt 0"):
    
    #read results
    outputFile:=fopen("integration/results.txt",READ,TEXT):
    myResults[1]:=readline(outputFile):
    i:=1:
    while (myResults[i] <> 0) do
      i:=i+1:
      myResults[i]:=readline(outputFile): #integral result
    od:
    close(outputFile):
  end if:
  print(StringTools[Join]([convert(intTime,string),"s. avg. spent on Maple integration."], " "));
  
  #compare the forms
  errors:=0:
  errorFile:=fopen("integration/errors.log",APPEND,TEXT):
  for i from 1 to polyCount do
    #print(parse(myLinForms[i]));
    curForms:=Array(parse(myLinForms[i])):
    #print(mapleLinForms[i]);
    myResults[i]:=parse(myResults[i]):
    wrong:=0: #prevents double counting errors, hopefully
    if decomposing = 1 then #check that decomposition is correct
      if nops(parse(myLinForms[i])) <> nops(mapleLinForms[i]) then
        print("Different number of powers of linear forms.");
        print(nops(parse(myLinForms[i])), nops(mapleLinForms[i]));
        errors:= errors + 1;
        wrong:=1:
      else
        mapleLinForms[i]:=convert(mapleLinForms[i], 'set');
        curTerm:={};
        for j from 1 to nops(parse(myLinForms[i])) do
          #print(curForms[j][1], curForms[j][2][1]);
          curTerm:=curTerm union {[curForms[j][1] / curForms[j][2][1]!, curForms[j][2]]};
          #print({[curForms[j][1] / curForms[j][2][1]!, curForms[j][2]]});
        od:
        if curTerm <> mapleLinForms[i] then
          print("Powers of linear forms don't match.");
          print(curTerm minus mapleLinForms[i], mapleLinForms[i] minus curTerm);
          errors:=errors + 1;
          wrong:=1:
        end if:
      end if:
    end if:
    if wrong = 0 then
      if mapleResults[i] <> simplify(myResults[i][1] / myResults[i][2]) then
        writeline(errorFile, "Integral calculation mismatch.");
        writeline(errorFile, "Forms:");
        writeline(errorFile, convert(mapleLinForms[i], string));
        writeline(errorFile, "Simplex:");
        writeline(errorFile, convert(mySimplices[i], string));
        writeline(errorFile, "Maple result:");
        writeline(errorFile, convert(mapleResults[i], string));
        writeline(errorFile, "C++ result:");
        writeline(errorFile, convert(simplify(myResults[i][1] / myResults[i][2]), string));
        errors:=errors + 1;
      end if:
    end if;  
  od:
  close(errorFile):
  
  print(StringTools[Join]([convert(errors,string),"tests failed."], " "));
end:

local benchmarks:
local myDim, myDegree:
benchmarks:=fopen("integration/benchmarks.txt",WRITE,TEXT):
errorFile:=fopen("integration/errors.log",WRITE,TEXT):
close(errorFile):
writeline(benchmarks, "Integration of monomial sums over random simplices (average time over 50 random monomials)"):
writeline(benchmarks, ""):
writeline(benchmarks, "                                   Degree                                             "):
writeline(benchmarks, "     _________________________________________________________________________________"):
writeline(benchmarks, "  n  |       1         2         5         10        20        30        40"):
writeline(benchmarks, "______________________________________________________________________________________"):
fprintf(benchmarks, "%3.3s  |", "2"):
close(benchmarks):
for myDim from 2 to 40 do
  for myDegree from 1 to 40 do
    #samplesize, bigConstant, numTerms, dimension, myDegree, decomposing
    #test_integration(50, 1000, 1, myDim, myDegree, 1):
    test_integration(50, 100, 1, myDim, myDegree, 0):
    if myDegree = 2 then
      myDegree:= 4:
    elif myDegree = 5 then
      myDegree:=9:
    elif myDegree = 10 then
      myDegree:=19:
    elif myDegree = 20 then
      myDegree:=29:
    elif myDegree = 30 then
      myDegree:=39:
    end if:
  od:
  if myDim = 8 then
      myDim:= 9:
  elif myDim = 10 then
    myDim:=14:
  elif myDim = 15 then
    myDim:=19:
  elif myDim = 20 then
    myDim:=29:
  elif myDim = 30 then
    myDim:=39:
  end if:
  benchmarks:=fopen("integration/benchmarks.txt",APPEND,TEXT):
  writeline(benchmarks, ""):
  fprintf(benchmarks, "%3.3s  |", convert(myDim + 1, string)):
  close(benchmarks):
od:

benchmarks:=fopen("integration/benchmarks.txt",APPEND,TEXT):
writeline(benchmarks, "Integration of a power of a linear form over random simplices (average time over 50 random forms)"):
writeline(benchmarks, ""):
writeline(benchmarks, "                                    Degree                                             "):
writeline(benchmarks, "      _________________________________________________________________________________"):
writeline(benchmarks, "   n  |       2        10        20         50       100       300      1000"):
writeline(benchmarks, "_______________________________________________________________________________________"):
fprintf(benchmarks, "%4.4s  |", "10"):
close(benchmarks):
for myDim from 10 to 400 do
  for myDegree from 2 to 1000 do
    #samplesize, bigConstant, numTerms, dimension, myDegree, decomposing
    test_integration(50, 100, 1, myDim, myDegree, 0):
    if myDegree = 2 then
      myDegree:= 9:
    elif myDegree = 10 then
      myDegree:=19:
    elif myDegree = 20 then
      myDegree:=49:
    elif myDegree = 50 then
      myDegree:=99:
    elif myDegree = 100 then
      myDegree:=299:
    elif myDegree = 300 then
      myDegree:=999:
    end if:
  od:
  if myDim = 10 then
      myDim:= 19:
  elif myDim = 20 then
    myDim:=49:
  elif myDim = 50 then
    myDim:=99:
  elif myDim = 100 then
    myDim:=299:
  elif myDim = 300 then
    myDim:=399:
  end if:
  benchmarks:=fopen("integration/benchmarks.txt",APPEND,TEXT):
  writeline(benchmarks, ""):
  fprintf(benchmarks, "%3.3s  |", convert(myDim + 1, string)):
  close(benchmarks):
od:
