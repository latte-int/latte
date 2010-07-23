read("integration/integrationTestsLib.mpl");

global totalErrors:
local benchmarks:
local myDim, myDegree:
errorFile:=fopen("integration/errors.log",WRITE,TEXT):
close(errorFile):
totalErrors:= 0:
for myDim from 2 to maxDim do
  for myDegree from 2 to maxDegree do
    
    printf("Integrating monomials of degree %d, dimension %d...\n", myDegree, myDim):
    #parameters:            polyCount, bigConstant, numTerms, dimension, myDegree, decomposing, randomGen)
    totalErrors:=test_integration(10, 1000, 10, myDim, myDegree, 1, random_sparse_homogeneous_polynomial_with_degree):
    
    if (totalErrors > 0) then
      quit:
    end if:
    printf("Finished integrating monomials of degree %d, dimension %d...\n", myDegree, myDim):
    printf("Integrating powers of linear forms of degree %d, dimension %d...\n", myDegree, myDim):
    totalErrors:=test_integration(10, 1000, 10, myDim, myDegree, 0, random_linearform_given_degree_dimension_maxcoef_componentmax_maxterm):
    
    if (totalErrors > 0) then
      quit:
    end if:
    printf("Finished integrating powers of linear forms of degree %d, dimension %d...\n", myDegree, myDim):
  od:
od:
#print(integral_power_linear_form([[0],[1]],1,1,[1]);
printf("%d total errors.\n", totalErrors):

 
