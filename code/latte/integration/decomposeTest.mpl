with(linalg):with(LinearAlgebra):

random_sparse_homogeneous_polynomial_with_degree:=proc(N,d,M,r) 
  local p, R;
  ## Give up if too large polynomials requested
  if (r > 500000) then
    error "Too large a polynomial requested"
  fi;
  R := rand(N);
  p := randpoly([ seq(x[i], i=1..d) ], 
                homogeneous, degree = M, terms = r, coeffs = proc() R() + 1; end);
end:

# Decomposing  a monomial in powers of linear forms.
# INPUT:L list of integers, m a list of integers
# OUTPUT: a nonnegative number
# MATH: The list $m$ represents the monomial x^m=x_1^{m[1]}x_2^{m[2]}\cdots x_d^{m[d], 
# the list L represents the linear form L=L[1]x[1]+..L[d]x_d. 
# The output is the coefficient of L^M,  M=M:=add(m[i],i=1..nops(m)), in the expansion of x^m in linear form using the formula in userguide. 
# CAUTION: The algorithm works only if we take "primitive" linear form as in the formula. 
coeff_linear_form_expansion:=proc(L,m) local p,i,j,k,c,s,M,out;
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
    s:=add(L[i],i=1..nops(L)); #print("s,m,c",s,m,c);
    #out:=1/M!*add((-1)^(M-k*s)*product(binomial(m[i], k*L[i]), i=1..nops(L))*(k^M), k=1..c);
    out:=add((-1)^(M-k*s)*product(binomial(m[i], k*L[i]), i=1..nops(L))*(k^M), k=1..c); #1/M! is implied in the C++ implementation
    end if;
#out;
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
# OUTPUT: a list of lists of length nops(m)
# MATH:  The list $m$ represents the monomial x^m=x_1^{m[1]}x_2^{m[2]}\cdots x_d^{m[d].
# The output is a list of lists. Each element in the list represents a linear form ([1,2]=x+2y).
# The output exhausts all the linear forms with exponents M=m[1]+..+m[d] which appear when  expressing   x^m as  linear combinations of linear forms with exponent M. The first element is the coefficient multiplied by coe
list_and_coeff_for_monome:=proc(m,coe) local M,L,out:
M:=add(m[i],i=1..nops(m));
L:=list_for_simplex_integral(m); #print("permutations", L);
out:=[ seq([coe*coeff_linear_form_expansion(L[j],m),[M,L[j]]],j=1..nops(L))];
end:

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

local polyCount:=10:
local bigConstant:=10000:
local numTerms:=10:
local dimension:=5:
local myDegree:=10:
local errors, curTerm, termCount, formList, termList, curPoly, myPolys, myList, myForms:
local myTime, temp:

#get polynomials
myTime:=0:
for myIndex from 1 to polyCount do
  curPoly:=random_sparse_homogeneous_polynomial_with_degree(bigConstant, dimension, myDegree, numTerms);
  if type(curPoly,`+`) then myList:=convert(curPoly,list) else myList:=[curPoly] end if;
  termCount:=nops(myList);
  formList[myIndex]:=[]:
  for i from 1 to termCount do
    constants:=subs(seq(k=1,k=[seq(x[i],i=1..dimension)]),myList[i]);
    termList:=[]:
    for j from 1 to dimension do
      termList:=ListTools[FlattenOnce]([termList, degree(myList[i], x[j])]);
    od;
    termList;
    formList[myIndex]:=ListTools[FlattenOnce]([formList[myIndex], [[constants, termList]]]);
  od;
  myPolys[myIndex]:=curPoly;
  temp:=time();
  myForms[myIndex]:={op(list_integral_via_waring(formList[myIndex]))};
  myTime:=myTime + time() - temp:
od:
myTime:=myTime / polyCount:

#write to file
polyFile:=fopen("integration/randomPolys.txt",WRITE,TEXT):
for i from 1 to polyCount do
  writeline(polyFile, StringTools[DeleteSpace](convert(formList[i], string)));
od:
close(polyFile):

#run the integrate program
system("./integrate_test integration/randomPolys.txt integration/forms.txt"):

print(StringTools[Join]([convert(myTime,string),"s. avg. spent on Maple decomposition."], " "));

#read forms in maple notation
formFile:=fopen("integration/forms.txt",READ,TEXT):
formList[1]:=readline(formFile):
i:=1:
while (formList[i] <> 0) do
  i:=i+1:
  formList[i]:=readline(formFile):
od:
close(formFile):

#compare the forms
errors:=0:
for i from 1 to polyCount do
#we can use the line below to make maple expand our linear forms and compare them to the original polynomial
#  if myPolys[i] = simplify(parse(formList[i])) then print("All Clear") else print("Issue here.") end if;
#  if {op(myForms[i])} = {op(parse(formList[i]))} then print("Equal") else print("Not Equal") end if;
   if {op(myForms[i])} <> {op(parse(formList[i]))} then errors:=errors + 1 end if;
od:

print(StringTools[Join]([convert(errors,string),"tests failed."], " "));
