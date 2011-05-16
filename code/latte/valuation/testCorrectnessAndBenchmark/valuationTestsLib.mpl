#	ValuationTestsLib.mpl
#   Created on: July 22, 2010
#      Author: Brandon Dutra and Gregory Pinto
#
# Goal: This maple file should contain all the functions need for testing the valuation class
#  through maple.
# Currently, we only test integration of hyper-rectangles.
#
# Assumes this script is called/exicuted in the top-level code director (two folders above it's location).
#


#dependent files: (used for making the polynomials.)
read("integration/integrationTestsLib.mpl"):


#Start of functions:


# Makes a latte file for the hyper-rectangle defined by the lower/upper integration bounds.
# @parm lowerIntegrationBound, upperIntegrationBound: vectors such that lowerIB[i] <= x[i] <= upperIB[i]
# @parm polytopeDimension: number of vars in the polynomial.
# @parm latteFileName: where we save the latte file for the hyper-rectangle.
# @return none
# @side affect: makes (or overwrites) a file on the system.
make_hyperrectangle_latte_file:=proc(lowerIntegrationBound, upperIntegrationBound, polytopeDimension, latteFileName)
	local latteFile, i, k;

	latteFile:=fopen(latteFileName,WRITE,TEXT):
	fprintf(latteFile, "%d %d\n", 2*polytopeDimension, polytopeDimension + 1);
 

	for i from 1 to polytopeDimension do
  		# write the lowerbound < x_i constraint. if a/b < x, write 0 < -a + bx.
		fprintf(latteFile, "%d ", -1 * numer(lowerIntegrationBound[i]));
  	
		for k from 1 to polytopeDimension do
			if k = i then
				fprintf(latteFile, "%d ",  denom(lowerIntegrationBound[i]));
			else
				fprintf(latteFile, "0 ");
			end if:
		od:
  	
		fprintf(latteFile, "\n");
  	
		#write the x_i < upperbound constraint. if x < a/b, write 0 < a - bx
		fprintf(latteFile, "%d ", numer(upperIntegrationBound[i]));
  	
		for k from 1 to polytopeDimension do
			if k = i then
				fprintf(latteFile, "%d ", -1*denom(upperIntegrationBound[i]));
			else
				fprintf(latteFile, "0 ");
			end if:
		od:
  	
		fprintf(latteFile, "\n");
	od: #write two lines to the latte file.
	close(latteFile);
end: 




# Integrates a maple polynomial over a rectangle defined in lower/upperIntegrationBound vector.
# @parm upperIntegrationBound: A vector, bounds such that x[i] <= upperIntegrationBound[i]
# @parm lowerIntegrationBound: A vector, bounds such that lowerIntegrationBound[i] <= x[i]
# @parm thePolynomial: A regular maple polynomial
# @parm polytopeDimension: the number of variables in the polynomial.
# @return a number (rational)
integrate_polynomial_over_rectangle:=proc(upperIntegrationBound, lowerIntegrationBound, thePolynomial, polytopeDimension)
	local integratedPoly, i;

	integratedPoly:= int(thePolynomial, x[1]=lowerIntegrationBound[1]..upperIntegrationBound[1]);
  
  	for i from 2 to polytopeDimension do
		integratedPoly:= int(integratedPoly, x[i]=lowerIntegrationBound[i]..upperIntegrationBound[i]);
	od:
  
	#printf("Inetegral of the rand poly:\n\n");
	#print(integratedPoly);	
	
	integratedPoly;
end:	


# Saves a maple polynomial to a file after converting it to the [ [c, [eee]]...] form.
# @parm maplePolynomial: a regular polynomial (not the special [ coef [exps]] list)
# @parm polynomialFileName the name of the file to save the polynomial after converting it to the list form.
# @parm dimension: number of variables in the polynomial.
# @return none
# @side affect: makes (or overwrites) a file on the system.
save_polynomial_to_file_mapleEncoded:=proc(maplePolynomial, polynomialFileName, dimension)
	local polyFile, stringPolynomial;
	
	polyFile:=fopen(polynomialFileName,WRITE,TEXT):
	
	stringPolynomial:=polynomial_to_sparsepoly(maplePolynomial, dimension);
	
	writeline(polyFile, convert(stringPolynomial, string));
	
	close(polyFile);
	
end:


#
# Makes a random rectangular polytope in polytopeDimension-dimension. The rational vertices
# are no more than integrationLimit. 
# Then we make a random polynomial of dimension=polytopeDimension and at most polyMaxDegree degree.
# For each r <= polyMaxDegree, the polynomial will have at most maxNumberOfTermsPerDegree-degree monomials whose rational coefficients are no more than 5000 to -5000.
#
# filename is the name of the latte file we save the polytope in, and the polynomial is saved in filename.polynomial
#
# Once that is done, we compute the real integral and call ./testPolytopeIntegration for comparison.
#
# if rationalCoeff is 1, then the random polynomial has rational coefficients, if it is 0, it has integer coefficients.
# @return int: the system status.
test_hyperrectangle_integtation:=proc(polyMaxDegree, polytopeDimension, maxNumberOfTermsPerDegree, integrationLimit, fileName, rationalCoeff) 
 
	local randomPoly, lowerIntegrationBound, upperIntegrationBound, randNumber, randNumber2, positive, i, integratedPoly, correctAnswer, polynomialFileName, systemCommand, correctAnswerString, status:
	  
	randNumber:= rand(integrationLimit + 1);
	randNumber2:= rand(100); #used for the denominator.
	positive := rand(2); #positive() is 0 or 1.
  
   
	#get a random polynomial.
	randomPoly:=random_sparse_nonhomogeneous_polynomial_with_degree_mapleEncoded(5000, polytopeDimension, polyMaxDegree, maxNumberOfTermsPerDegree, rationalCoeff):
  
	lowerIntegrationBound:=Vector(polytopeDimension);
	upperIntegrationBound:=Vector(polytopeDimension);
  

	#find the upper and lower integration bounds.
	for i from 1 to polytopeDimension do
		lowerIntegrationBound[i] := (randNumber() / (randNumber2() + 1)) * (-1)^positive();
		upperIntegrationBound[i] := lowerIntegrationBound[i] + randNumber() + 1;
	od:
  
  
	print(Transpose(lowerIntegrationBound));
	print(Transpose(upperIntegrationBound));
  
  
	#integrate the polynomial.
	correctAnswer := integrate_polynomial_over_rectangle(upperIntegrationBound, lowerIntegrationBound, randomPoly, polytopeDimension):
    #printf("Found the correct answer\n");
    
    printf("Making latte file\n");
	#now make the latte file.
	make_hyperrectangle_latte_file(lowerIntegrationBound, upperIntegrationBound, polytopeDimension, fileName);
	#printf("made the latte file\n");
	
	#make the polynomial file
	polynomialFileName:=fileName||".polynomial"; #concat. the strings.
	save_polynomial_to_file_mapleEncoded(randomPoly, polynomialFileName, polytopeDimension);
	#print(polynomialFileName);
	
	#Finally, now test our code.
	correctAnswerString :=convert(correctAnswer, string):
	systemCommand:= "./test-hyperrectangle-integration " || correctAnswerString || " " || polynomialFileName || " " || fileName :
	#print(systemCommand);
	status:=system(systemCommand):
	printf("status=%d\n", status);
	status; #return the status.
end:
