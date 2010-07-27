read("integration/integrationTestsLib.mpl");

#comment this line out if you do not want a random value.
seed := randomize(): #calls to rand() will now be different on different maple runs. 
printf("seed = %d\n\n", seed);


global totalErrors:
local benchmarks:
local myDim, myDegree:
errorFile:=fopen("integration/errors.log",WRITE,TEXT):
close(errorFile):
totalErrors:= 0:
for myDim from 2 to maxDim do
  for myDegree from 2 to maxDegree do
    
    printf("Integrating monomials of degree %d, dimension %d...\n", myDegree, myDim):
    #parameters:            polyCount, bigConstant, numTerms, dimension, myDegree, decomposing, randomGen, rationalCoeff)
    totalErrors:=test_integration(10, 1000, 10, myDim, myDegree, 1, random_sparse_homogeneous_polynomial_with_degree, 1):
    
    if (totalErrors > 0) then
      printf("seed = %d\n\n", seed);
      quit:
    end if:
    
    printf("\n");
    
    printf("Integrating powers of linear forms of degree %d, dimension %d...\n", myDegree, myDim):
    totalErrors:=test_integration(10, 1000, 10, myDim, myDegree, 0, random_linearform_given_degree_dimension_maxcoef_componentmax_maxterm, 1):
    
    if (totalErrors > 0) then
      printf("seed = %d\n\n", seed);
      quit:
    end if:
    
    printf("\n");
  od:
od:

printf("seed = %d\n\n", seed);
#print(integral_power_linear_form([[0],[1]],1,1,[1]);
printf("%d total errors.\n", totalErrors);

 
