read("testingTopEhrhart_lib.mpl"); #load the test functions
read("BirkhoffGenerator.mpl"); #read the Birkhoff Simplex Generator Code.



#
#  Oct 12, 2010
#  Author: Brandon, Gregory
#  Description: Uses the main TopEhrhart lib to compute the top-three ehrhart coeff of the Birkhoff Simplex.
#  Writes a log of the Ehrhart polytopes computed along with the non-projected polytopes.
#  The non-projected polytopes are saved in a "database" file so that we do not test the same polytope twice.
#


#Parameters:
#	baseFileName: base name of the data files.
#	dim: the full-dimension 
#	numberToTest: number of test cases to run in the same dimension.
baseFileName:="debug":#"birkhoffNonProjectDim5";
numberToTest:=5000;
dim:=4;

#set the other parameters that depend on the last ones above.
fileNameDatabaseSimplex:=cat(baseFileName,".database");
fileNameDatabaseZeroDivider:=cat(baseFileName, ".zerodivider");
fileNameLog:=cat(baseFileName,".log");
logFile_ptr:=fopen(fileNameLog, APPEND, TEXT);
fprintf(logFile_ptr, "******************************************\n");
numFailedTest:=0;
numFinishedTest:=0;


#Important:
#	if this is the first time testing a new dimension (a database does not exist), then please set 1=1 so that 
#	this if-statement runs. This if-statement makes a new database file and inserts a dummy value into it
#	because I could not get maple to read an empty file.
if 1 = 1 then:
	#first, make a dummy database file.
	#the database cannot start empty, so lets add a "polytope" that cannot exist: a polytope who's matrix is all 1.
	fakePolytopeVertex:=[]:
	for i from 1 to dim*dim do:
		fakePolytopeVertex:=[op(fakePolytopeVertex), 1]:
	end: #fakePolytopeVertex==[1,1,1,1...1];
	nonProjectedDatabase:=[[{fakePolytopeVertex}]];
	save nonProjectedDatabase, fileNameDatabaseSimplex: #this also makes/clears the db file.
end: #used if want to make a new db file.
 

#the rest of this script does not need to be edited.

#This while loop does all the work, here is its logic:
#	for ever test case 
#		make a new simplex. save the full-dimension and non-projected simplex.
#		
while numFinishedTest < numberToTest do;
	#print to screen
	printf("Number of repeaded slimplex: %d\n", numFailedTest);
	printf("Starting test %d out of %d\n", numFinishedTest+1, numberToTest);

	#print to log
	fprintf(logFile_ptr,"Number of repeaded slimplex: %d\n", numFailedTest);
	fprintf(logFile_ptr,"Starting test %d out of %d\n", numFinishedTest+1, numberToTest);
	fflush(logFile_ptr);

	#find a simplex. save the full-dimension and non-projected simplex.
	twoSimplex:=gensimplexB_n_nonProjected(dim,(dim -1)^2 +1,[]):
	simplexForTopEhrhart:=twoSimplex[1]:
	simplexNonProjected:=twoSimplex[2]:

	simplexNonProjected:=convert(simplexNonProjected, set):
	#print("simplexNonProjected=", simplexNonProjected);


	#read the database.
	nonProjectedDatabase:=0:
	read fileNameDatabaseSimplex:

	#check that our simplex is not a repeat. 
	repeadedSimplex:= 0:
	for i from 1 to nops(nonProjectedDatabase) do:
		if nonProjectedDatabase[i][1] subset simplexNonProjected  and simplexNonProjected subset nonProjectedDatabase[i][1] then:
			repeadedSimplex:=1:
			break;
		end: #if repeaded simplex.
	end:

	if repeadedSimplex = 1 then:
		numFailedTest:=numFailedTest+1:
		next: #try again. This is like a c/c++ "continue"
	end:

	#simp is a new simplex!
	#save the new simplex.
	nonProjectedDatabase:=[op(nonProjectedDatabase), [simplexNonProjected]]:
	save nonProjectedDatabase, fileNameDatabaseSimplex:

	#run the top-ehrhart functions on it.
	printf("Calling top ehrhart functions\n");
	fprintf(logFile_ptr, "calling top ehrhart function\n");
	fflush(logFile_ptr);
	divideByZero:=1: #start off as true.
	while divideByZero = 1 do
		try  
			ans:=test_top_ehrhart_given_simplex_v2(simplexForTopEhrhart):
			divideByZero:=0: #it worked ok.
		catch :
			fprintf(logFile_ptr, "I think we divided by zero: %q\n", lastexception);
		end try;
	end; #end loop.

	#print to screen
	print("this simplex", simplexForTopEhrhart);
	print(" gave the ans", ans);

	#print to log.
	fprintf(logFile_ptr, "Test Finished: %d out of %d\n",numFinishedTest+1, numberToTest);
	fprintf(logFile_ptr, "%q\n nonprojected: %q\n", convert(ans, string), convert(simplexNonProjected, string));
	fflush(logFile_ptr);

	#now test to see if this simplex is special or not.
	polynomial:=ans[2]:
	polynomialCoeff:=[coeffs(polynomial, t)]:

	for i from 1 to nops(polynomialCoeff) do
		if polynomialCoeff[i] < 0 then
			printf("OH MY GOSH, IT IS NEGATIVE!!!!\n");
			quit;
		end;
	end:

	numFinishedTest:=numFinishedTest+1;
end; #end while loop 

fclose(logFile_ptr);


