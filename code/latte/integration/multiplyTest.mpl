with(linalg):with(LinearAlgebra):
with(numapprox,laurent):

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

# We rather use a built-in procedure of Maple,
# and convert the result to our format.

##
## Creating random polynomials in sparse format
##

# N - use coefficients up to this number
# d - dimension
# M - degree
# r - number of monomials

random_sparse_homogeneous_polynomial_with_degree:=proc(N,d,M,r) 
  local p, R;
  ## Give up if too large polynomials requested
  if (r > 500000) then
    error "Too large a polynomial requested"
  fi;
  R := rand(N);
  p := randpoly([ seq(x[i], i=1..d) ], 
                homogeneous, degree = M, terms = r, coeffs = proc() R() + 1; end);
  #polynomial_to_sparsepoly(p, d);
end:

local myPolys, myResults, myIndex:
local benchmarks:
local myDim, myDegree:
errorFile:=fopen("integration/errors.log",WRITE,TEXT):
close(errorFile):
local polyCount:=10:
local temp, polyFile, formFile, i:

myIndex:=1:
polyFile:=fopen("integration/randomPolys.txt",WRITE,TEXT):
for myDim from 10 to 30 do
  for i from myIndex to myIndex + 2*polyCount - 1 do
  #polyCount, bigConstant, numTerms, dimension, myDegree
  print(StringTools[Join](["Multiplying monomials of dimension", convert(myDim,string)], " ")):
  myPolys[i]:=random_sparse_homogeneous_polynomial_with_degree(1000, myDim, myDim*3, 10);
  writeline(polyFile, convert(polynomial_to_sparsepoly(myPolys[i], myDim), string));
  myIndex:=myIndex+1:
  od:
od:
close(polyFile):

#run the integrate program
system("./multiply_test integration/randomPolys.txt integration/results.txt"):

#read results in maple notation
formFile:=fopen("integration/results.txt",READ,TEXT):
myResults[1]:=readline(formFile):
i:=1:
while (myResults[i] <> 0) do
  i:=i+1:
  myResults[i]:=readline(formFile):
od:
close(formFile):

#compare the forms
errors:=0:
myIndex:=1:
for myDim from 10 to 30 do
  for i from myIndex to polyCount - 1 do
  curPoly:=expand(myPolys[2*i - 1] * myPolys[2*i]);
  curPoly:=polynomial_to_sparsepoly(curPoly, myDim);
  if nops(curPoly) <> nops(parse(myResults[i])) then
    print("Different number of terms.");
    print(nops(curPoly), nops(myResults[i]));
    errors:= errors + 1;
  else
    myResults[i]:=convert(parse(myResults[i]), 'set');
    curPoly:=convert(curPoly, 'set');
    if curPoly <> myResults[i] then
      print("Products don't match.");
      print(curPoly, myResults[i]);
      errors:=errors + 1;
    end if;
  end if;
  myIndex:=myIndex+1:
  od:
od:

if errors <> 0 then
print("Multiply error!"):
end if;
