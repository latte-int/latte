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

local benchmarks:
local polyCount:=50:
local myDim, myDegree:
errorFile:=fopen("integration/errors.log",WRITE,TEXT):
close(errorFile):

for myDim from 1 to 10 do
  print(StringTools[Join](["Multiplying monomials of dimension", convert(myDim*5,string)], " ")):
  polyFile:=fopen("integration/randomPolys.txt",WRITE,TEXT):
  for i from 1 to polyCount do
    writeline(polyFile, convert(polynomial_to_sparsepoly(random_sparse_homogeneous_polynomial_with_degree(10000, myDim*5, 10, 50), myDim*5), string)); 
  od:
  close(polyFile):
  system("./compare -m integration/randomPolys.txt"):
od:
  
for myDim from 1 to 10 do
  print(StringTools[Join](["Decomposing monomials of dimension", convert(myDim,string)], " ")):
  polyFile:=fopen("integration/randomPolys.txt",WRITE,TEXT):
  for i from 1 to polyCount do
    writeline(polyFile, convert(polynomial_to_sparsepoly(random_sparse_homogeneous_polynomial_with_degree(10000, myDim, 10, 5), myDim), string)); 
  od:
  close(polyFile):
  system("./compare -d integration/randomPolys.txt"):  
od:
