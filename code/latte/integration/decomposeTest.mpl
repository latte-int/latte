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

local polyCount:=1:
local bigConstant:=10:
local numTerms:=10:
local maxDim:=50:
local maxDegree:=10:
local dimension:=5:
local myDegree:=20:
local curTerm, termCount, formList, termList, curPoly, myPolys, myList:

#get polynomials
for k from 1 to polyCount do
  curPoly:=random_sparse_homogeneous_polynomial_with_degree(bigConstant, dimension, myDegree, numTerms);
  if type(curPoly,`+`) then myList:=convert(curPoly,list) else myList:=[curPoly] end if;
  termCount:=nops(myList);
  formList[k]:=[]:
  for i from 1 to termCount do
    constants:=subs(seq(k=1,k=[seq(x[i],i=1..dimension)]),myList[i]);
    termList:=[]:
    for j from 1 to dimension do
      termList:=ListTools[FlattenOnce]([termList, degree(myList[i], x[j])]);
    od;
    termList;
    formList[k]:=ListTools[FlattenOnce]([formList[k], [[constants, termList]]]);
  od;
  myPolys[k]:=curPoly;
od:
#write to file
polyFile:=fopen("code/latte/integration/randomPolys.txt",WRITE,TEXT):
for i from 1 to polyCount do
  writeline(polyFile, StringTools[DeleteSpace](convert(formList[i], string)));
od:
close(polyFile):

#run the integrate program
system("code/latte/integrate code/latte/integration/randomPolys.txt code/latte/integration/forms.txt"):

#read forms
formFile:=fopen("code/latte/integration/forms.txt",READ,TEXT):
formList[1]:=readline(formFile):
i:=1:
while (formList[i] <> 0) do
  i:=i+1:
  formList[i]:=readline(formFile):
od:

#compare the forms
for i from 1 to polyCount do
  if myPolys[i] = simplify(parse(formList[i])) then print("All Clear") else print("Issue here.") end if;
#  myPolys[i];
#  simplify(parse(formList[i]));
od;
close(formFile):
