with(linalg):
with(LinearAlgebra):
with(numapprox,laurent):
read("integration/createLinear.mpl"): #load the function to make a random linear form
read("integration/integrationTestsLib.mpl"):#integration functions.
read("testingSL/SL_lib.mpl");	#load the SL functions.

# A B C D E F G H I J K L M N O P Q R S T U V W X Y Z


#Input: simplexDim: the amb. dim of the simplex.
#Output: returns a list of simplexDim+1 vectors in R^(simplexDim).
create_random_simplex:=proc(simplexDim)
	local i, j, M, checkRankMatrix:

	do
		M:=RandomMatrix(simplexDim+1, simplexDim) /~ ( abs(RandomMatrix(simplexDim+1, simplexDim)) +~ 1);
		
		checkRankMatrix:=RandomMatrix(simplexDim, simplexDim);
		for i from 1 to simplexDim do
			for j from 1 to simplexDim do
				checkRankMatrix[i, j] := M[i, j] - M[simplexDim+1, j];
			od;
		od;
		#print("got here");
		#print(checkRankMatrix, "=checkRandMatrix");
		
		if Rank(checkRankMatrix) = simplexDim then:
			break;
		fi;
	end do;
	#print(M, "=M");
	
	return(convert(M, listlist));
end:



#Input: A list of numSimplex many simplex
#output: A list of the volume of each as computed by the SL functions.
find_volume_using_SL:=proc(simplexList, numSimplex, simplexDim)
	local volumeList, oneVolume, i, dummyLinearForm, startTime:
	
	volumeList:=[];
	dummyLinearForm:=[];
	#set dummyLinearForm=[0,0....0].
	for i from 1 to simplexDim do:
		dummyLinearForm:=[ 0, op(dummyLinearForm)];	
	od;
	
	#Find the volume of each simplex with the SL functions.
	for i from 1 to numSimplex do:
		printf("SL Valuation: going to compute volume of simplex %d out of %d\n", i, numSimplex);
		#tfunction_SL_simplex_ell(t,S,L,ell,M)
		startTime:=time();
		oneVolume:=tfunction_SL_simplex_ell(1,simplexList[i],simplexList[i], dummyLinearForm,0);
		printf("SL Volume: %d/%d\n", numer(oneVolume), denom(oneVolume));
		printf("SL Time: %f sec\n", time() - startTime);
		volumeList:=[oneVolume, op(volumeList)];
	od;
	return(volumeList);
end:




#Input: n+1 points in R^n
#Output: the facet equations that define the simplex of the n+1 points.
simplex_to_hyperplanes:=proc(simplex)
local i, b_ax, equations, n_points;

	equations:= [];
	for i from 1 to nops(simplex) do:
		n_points:= subsop(i = NULL, simplex);
		b_ax:= facets_equation(n_points, simplex[i]); # 0 < b - a*x
		equations:=[ op(b_ax), op(equations)];
	od;
	return(equations);
end:
	




test_sl_countLattice:=proc(simplexDim, numTests, fileName)
	local i, filePtr, mySimplices, stringCommand;
	local countAnswers, SLAnswers;
	
	countAnswers:=[];
	SLAnswers:=[];
	latteFacetFile:=filName||".latte";
	countFile:=fileName||".count";
	
	filePtr:=fopen(countFile, WRITE, TEXT);
	fprintf(filePtr, " ");
	fclose(filePtr);
	
	#make the simplex
	for i from 1 to numTests do:
		mySimplices[i]:=lattice_random_simplex(simplexDim, 10000);
	od:
	
	#call count, save answer in new line of the 'count file'
	for i from 1 to numTests do:
		printf("Calling count simplex %d of %d\n", i, numTests);
		equations:=simplex_to_hyperplanes(mySimplices[i]);
		write_facets_to_file(equations, fileName, simplexDim);
		stringCommand:="./count --redundancy-check=none "||fileName||" 2>/dev/null | tail -n 1 >> "||countFile;
		print(stringCommand);
		system(stringCommand);
	od:
	
	#read latte's answers.
	filePtr:=fopen(countFile, READ, TEXT);
	for i from 1 to numTests do:
		countAnswers:=[parse(readline(filePtr)), op(countAnswers)];
	od:
	fclose(filePtr);

	#find the count using SL.
	for i from 1 to numTests do:
		SLAnswers:= [1, op(SLAnswers)]; #TO FIX
	od: 
	
	#compare the answers.
	if nops(countAnswers) <> nops(SLAnswers) then:
		print("Different number of answers reported by countAnswers and SLAnswers");
		print(countAnswers, "=countAnswers");
		print(SLAnswers, "=SLAnswers");
		print(numTests, "=numTests");
		quit;
	fi;
	for i from 1 to numTests do:
		if countAnswers[i] <> SLAnswers[i] then:
			print("Different answers: count vs sl:", countAnswers[i], "vs.", SLAnswers[i]);
			print(countAnswers, "=countAnswers");
			print(SLAnswers, "=SLAnswers");
			print(numTests, "=numTests");
		fi;
	od:
end:



#Input:
#@parm: simplexDim: the abm. dim of the simplix
#@parm: numTests: how many simpleices you want to test at once
#@parm: baseFileName: string. File names used for saving latte's facet equations and latte's integration answer. ex:"testingSL/testingSL_volume"
test_sl_integration:=proc(simplexDim, numTests, degreeL)
############################################################################################################
	local randomGen:
	local mySimplices, mapleLinForms,  mapleResults, integrationSLanswer:
	local myIndex, formIndex, i, j:
	local rationalCoeff, bigConstant, seed, numTerms:
	local mapleTime, SLTime:
	
	randomGen:=random_linearform_given_degree_dimension_maxcoef_componentmax_maxterm; 

	bigConstant:=10000; #sets the max size of the coeff. and terms in the random linear form.
	numTerms:=1; #we only want to test 1 linear form at a time.
	rationalCoeff:=1;#the coeff. of the linear form should be rational.
	mapleTime:= 0;
	SLTime:=0;
  
	seed:=randomize();
	print(seed, "seed");
  
  	#make the simplex and linear forms.
	for myIndex from 1 to numTests do:
		mySimplices[myIndex]:=lattice_random_simplex(simplexDim, bigConstant);		
      	mapleLinForms[myIndex]:=randomGen(degreeL, simplexDim, bigConstant, bigConstant, numTerms, rationalCoeff):
	od:
	
	#print(mySimplices, "=mySimplices");
	#print(mapleLinForms, "=mapleLinForms");
	
  
  	mapleTime:=time();
  	#find maple and SL integraton. 
	for myIndex from 1 to numTests do
		printf("Integrating %d-linear form to the %d power using origional maple code. Test %d out of %d\n", simplexDim, degreeL, myIndex, numTests);
    	mapleResults[myIndex]:=0:
    	
    	for formIndex from 1 to nops(mapleLinForms[myIndex]) do
      		mapleResults[myIndex]:=mapleResults[myIndex]+mapleLinForms[myIndex][formIndex][1]*integral_power_linear_form(mySimplices[myIndex],simplexDim,mapleLinForms[myIndex][formIndex][2][1],mapleLinForms[myIndex][formIndex][2][2]):
    	od:
    	#print("[1]",mapleLinForms[myIndex][1]);
 		#print("[1][1]",mapleLinForms[myIndex][1][1]);
 		#print("[1][2]",mapleLinForms[myIndex][1][2]);	    
    od:
    mapleTime:=time()-mapleTime;
    printf("Total time using origional maple code: %f\n", mapleTime);

	SLTime:=time();    
    for myIndex from 1 to numTests do
		#tfunction_SL_simplex_ell(t,S,L,ell,M)
		printf("Integrating %d-linear form to the %d power using SL code. Test %d out of %d\n", simplexDim, degreeL, myIndex, numTests);
    	integrationSLanswer[myIndex]:=mapleLinForms[myIndex][1][1]*tfunction_SL_simplex_ell(1,mySimplices[myIndex],mySimplices[myIndex],mapleLinForms[myIndex][1][2][2],mapleLinForms[myIndex][1][2][1]);
	od:
  	SLTime:=time()-SLTime;
  	printf("Total time using SL code: %f\n", SLTime);
  
  
  #compare the forms
	if nops(mapleResults) <> nops(integrationSLanswer) then:
		print("Different number of test cases for maple results and SL results");
		quit;
	fi; 
	for i from 1 to numTests do
		if  mapleResults[i] <> integrationSLanswer[i]*degreeL! then:
		  	print("Test: ", i);
  			print("mapleResults[i]=\n", mapleResults[i]);
  			print("integrationSLanswer[i] =\n", integrationSLanswer[i]);
  			print("integrationSLanswer[i]*m! =\n", integrationSLanswer[i]*degreeL!);
  			print(integrationSLanswer[i]/mapleResults[i]);
    		quit;
    	fi;
  	od;
  	
  	#if got here, no errrors.
  	printf("NO INTEGRATION ERRORS\n");
end:




#This is the main function for testing volume between the SL functions and latte.
#Input: dim. of the simplex, number of simplices you want to test, and base file name for savine latte's facets and volume answers.
#Output: creates latte-style facet files and a volume file.
test_sl_volume:=proc(simplexDim, numTests, baseFileName)
local simplexList, i, volumeAnswersLatte, volumeAnswersSL, fileNameSimplex, fileNameVolume, systemCommand, status, equations:
local fileNameSimplex_i, stringI, DD, stringNum, seed;

	seed:=randomize();

	fileNameSimplex:=baseFileName || ".simplex" :
	fileNameVolume:=baseFileName || ".volume":
	
	#print(fileNameSimplex, "=fileNameSimplex");

	#make the test simplex
	simplexList:= [];
	
	for i from 1 to numTests do:
		DD:=[[0,0,0,0],[0,0,1,0],[1,0,0,0],[0,1,0,0],[0,0,0,1]]; #used for debugging.
		simplexList:=[create_random_simplex(simplexDim), op(simplexList)];
		#simplexList:=[DD, op(simplexList)];
	od;
	
	#print(simplexList, "=simplexList");
	
	
	for i from 1 to numTests do:
		printf("%d: Writing simplex facet equations out to file\n", i);
		equations:=simplex_to_hyperplanes(simplexList[i]);
		stringI:=convert(i, string);
		fileNameSimplex_i:=cat(fileNameSimplex,stringI);
		write_facets_to_file(equations, fileNameSimplex_i, simplexDim);	
	od;
	
	
	#Ask C++ exe to compute the volumes of the simplex.
	stringNum := convert(nops(simplexList), string);
	systemCommand:= "./testVolumeForSL " || fileNameSimplex  || " " || fileNameVolume || " " || stringNum || " 2>/dev/null":
	#print(systemCommand);
	status:=system(systemCommand):

	#read latte's answers
	volumeAnswersLatte:=read_volume_from_file(fileNameVolume, numTests);
	
	#compute SL's volume.
	volumeAnswersSL:=find_volume_using_SL(simplexList, numTests, simplexDim);
	
	#test the answers
	if nops(volumeAnswersLatte) <> nops(volumeAnswersSL)  or nops(volumeAnswersLatte) <> numTests then:
		print("Error, latte and SL functions have different number of results.");
		print(volumeAnswersLatte);
		print(volumeAnswersSL);
		exit(1);
	fi;
	
	for i from 1 to numTests do:
		if  volumeAnswersLatte[i] <> volumeAnswersSL[i] then:
			print("The ", i, "st test does not agree");
			print(volumeAnswersLatte);
			print(volumeAnswersSL);
			exit(1);
		fi;
	od;
	
	#if made it this far, the test went error free!!!
	printf("Seed: %d\n", seed);
	printf("Tested %d dim-%d simplices\n", numTests, simplexDim);
	printf("NO ERRORS!\n");
	
end:





#input: Read numberAnswers many rows from the file fileName. So the volumes should all be on their own lines
#output: List of volumes as read in by the file.
#usage: This function reads the output of the testVolumeForSL C++ exe.
read_volume_from_file:=proc(fileName, numberAnswers)
	local filePtr, volumeString, volumeList, i;

	volumeList:=[];
	filePtr:=fopen(fileName, READ, TEXT);
	for i from 1 to numberAnswers do:
		volumeString:=readline(filePtr);
		volumeList:=[ parse(volumeString), op(volumeList)];
	od:

	close(filePtr);
	#print(volumeList, "=read_volume.volumeList");
	return(volumeList);
end:

#Input:facet equations in maple syntax.
#Output: Writes to file fileName the latte-style facet equations.
write_facets_to_file:=proc(equations, fileName, simplexDim)
	local i, j, filePtr, lcmDenom, M:
	
	filePtr:=fopen(fileName, WRITE, TEXT);
	
	#print(equations, "=equations");
	#print(fileName, "=fileName");
	
	
	
	#sorry, I had to do a bunch of converting because I couldn't pull the elements in the origional structure.
	M:=convert(equations, Matrix);
	#print(M, "=M1");
	M:=convert(M, list);
	#print(M, "=M3");
	#print(simplexDim, "=simplexDim");
	M:=matrix(simplexDim+1, simplexDim+1, M);
	#print(M, "=M2");
	#print(M[1], "=M2[1,1]");
	
	#write the latte file.
	fprintf(filePtr,"%d %d\n", simplexDim+1, simplexDim+1);
	for i from 1 to simplexDim+1 do:
		lcmDenom:=1;
		
		for j from 1 to simplexDim+1 do:
			lcmDenom:=lcm(lcmDenom, denom(M[i, j]));
		od:

		for j from 1 to simplexDim+1 do:
			#print(M[i, j], "M[i, j]");
			fprintf(filePtr, "%d ", M[i, j]*lcmDenom);
		od:
		fprintf(filePtr, "\n");
	od:
	close(filePtr);
end:

#input:filename and list of simplices
#output: writes each simplex to a new line of the file in maple syntax.
write_simplex_to_file:=proc(simplexList, fileName)
	local filePtr, i:
	filePtr:=fopen(fileName,WRITE,TEXT):
	
	for i from 1 to nops(simplexList) do:	
		writeline(filePtr, convert(simplexList[i], string));
	od;	
	close(filePtr);
end:



                                                          
#Input:List of n points of a simplex and the last point.
#Output: the Facet of the equation containing the n points with the sign such that the last point is in the halfspace. 
#Description: Feeding the n points of a n-simplex and the 1 point that
# is not in their facet it generates ONE equation of the simplex in LattE format
#for those points.
facets_equation:=proc(L,notinL)
	local x,aux1,aux2,M,i;
	M:=Matrix([op(1,L)]);
	for i from 2 to nops(L) do
		M:=stackmatrix(M,op(i,L));
	od:
	M:=stackmatrix(M,[x[j] $j=1..nops(L)]);
	M:=transpose(M);
	M:=stackmatrix(M,[1 $j=1..nops(L)+1]);
	aux1:=det(M):
	for i from 1 to nops(notinL) do
		aux1:=subs(x[i]=op(i,notinL),aux1):
	od;
	#print(sign(aux1),"hi I am here");
	aux2:=sign(aux1)*det(M);
	for i from 1 to nops(notinL) do
		aux2:=subs(x[i]=0,aux2):
	od;
	return(genmatrix({aux2*x[0]+sign(aux1)*det(M)},[x[j] $j=0..nops(L)]));
end:

