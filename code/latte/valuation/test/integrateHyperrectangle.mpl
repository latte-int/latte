read("valuation/testCorrectnessAndBenchmark/valuationTestsLib.mpl"):


#	integrateHyperrectangle.mpl
#   Created on: July 22, 2010
#      Author: Brandon Dutra and Gregory Pinto
$
# This maple scripts is a driver for test_hyperrectangle_integtation
# many times, which is defined in valuationTestsLib.mpl.
#
# test_hyperrectangle_integtation() makes a n-dim hyper-rectangle,
# and a random polynomial. We compute the integral of the
# polynomial over the polytope in maple and in the PolytopeValuation
# class, and compare for any differences. 
#
# Note that we assume polynomial dimension = dimension of the polytope.
#

print("Testing integration of hyper-rectangles...");

local myDim, myDegree, fileName, totalErrors, returnStatus, seed:

 
seed := randomize():
printf("random seed = %d \n", seed);


fileName:="valuation/hyerrectangleIntegrationTest.latte":
for myDim from 2 to 5 do
	for myDegree from 2 to 10 do
	
		printf("Testing polynomials of dimension %d, and at most degree %d\n", myDim, myDegree);
  
#                   the parameters are:  polyMaxDegree, polytopeDimension, maxNumberOfTermsPerDegree, integrationLimit, fileName, rationalCoefficents)    
    	returnStatus:=test_hyperrectangle_integtation_polynomials(myDegree, myDim, 10, 30000, fileName, 1):
   
		if ( not(returnStatus = 0)) then 
			printf("\n\n An error was reported in integrating polynomials, please see the polytope and polynomial files.\n");
			`quit`(5); 
		end if:
		
		printf("Testing powers of linear forms of dimension %d, and at most degree %d\n", myDim, myDegree);
  
#                   the parameters are:  formMaxDegree, polytopeDimension, maxNumberOfTermsPerDegree, integrationLimit, fileName, rationalCoefficents)    
    	returnStatus:=test_hyperrectangle_integtation_linear_forms(myDegree, myDim, 10, 30000, fileName, 1):
   
		if ( not(returnStatus = 0)) then 
			printf("\n\n An error was reported in integrating powers of linear forms, please see the polytope and polynomial files.\n");
			`quit`(5); 
		end if:

		
	od:
od:

