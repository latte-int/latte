read("testingTopEhrhart_lib.mpl"); #load the test functions
with(CodeTools[Profiling]):
#
#  Sept. 9, 2010
#  Author: Brandon
#  Description: Uses the main SL testing lib to compute the volume of simplies and compare the answer's with latte's tirangulation method.
#

#Saves the average time it takes to find the top-3 ehrhart coeff in a file.
#For each random simplex of dim between startingDim to infinity, test numTests many ehrhart top 3 coefficients.
# The time of each simplex coeff. is saved in the log file, and the average of the numTests tests is saved in the table file.
# If the first test takes more than 1/2 hour, we stop the function.
#Input:
#@parm: startingDim: starting dim. of the simplex.
#@parm: numTests: how many simpleices you want to test at once
#@parm: baseFileName: string. File names used for saving the time, average, and latte's facet equations and latte's volume answer. ex:"testingSL/testingSL_volume"
table_time_sl_ehrhart:=proc(startingDim, numTests, fileBaseName)
	local currentDim, currentTestNumber, time_coeff, fileNameLog, fileNameTable, filePtrLog, filePtrTable;
	local failFlag, totalTime;
	
	#file names
	fileNameLog:=fileBaseName||".log"; #hold times of everything.
	fileNameTable:=fileBaseName||".average";	#hold times only for the average.
	
	#file ptr.
	filePtrLog:=fopen(fileNameLog, WRITE, TEXT);
	filePtrTable:=fopen(fileNameTable, WRITE, TEXT);
	fprintf(filePtrTable, "Number of tests | dim | average time\n");
	
	
	currentDim:=startingDim;
	#this loop is broken if a test takes more than 30mins.
	while 1 = 1 do
		printf("currently starting dim %d\n", currentDim);
		currentTestNumber:=1;
		totalTime:=0;
		
		#we first find the time for the currentDim test.
		#if this 1 test takes more than 30min, we stop the function.
		
		failFlag:= 1;
		while failFlag <> 0 do
			try
				time_coeff:=test_sl_ehrhart(currentDim, fileBaseName||".debug");
				failFlag:=0; #no erros. break out of while loop.				
			catch:
				printf("Something went wrong: %q\n",lastexception);
				failFlag:=1; #try again.
			end try;
		end; #while the first test case did not fail.

		#print to log.
		fprintf(filePtrLog, "Dim %d:\t test %d out of %d:\t time %f\n", currentDim, currentTestNumber, numTests, time_coeff[1]);
		fprintf(filePtrLog, "coeff: %q\n", convert(time_coeff[2], string));
		fflush(filePtrLog);
		#also print to screen.  
		printf("Dim %d: test %d out of %d: time %f\n", currentDim, currentTestNumber, numTests, time_coeff[1]);
		print(time_coeff);
		
		#finally, how long did 1 test take? If more than 30min*60sec/min, stop.
		if time_coeff[1] > 60*30 then:
			printf("Test took longer than 1/2 hour, stopping script\n");
			currentDim:=currentDim + 1;
			continue;
			#break; #break the while 1 = 1 loop.
		fi;
		
		#save the time.
		totalTime:=time_coeff[1];
		
		#get times for numTest -1 more examples!!!
		for currentTestNumber from 2 to numTests do:
			failFlag:= 1;
			while failFlag <> 0 do
				try
					time_coeff:=test_sl_ehrhart(currentDim, fileBaseName||".debug");
					failFlag:=0; #no erros. break out of while loop.
					
					#print to log.
					fprintf(filePtrLog, "Dim %d:\t test %d out of %d:\t time %f\n", currentDim, currentTestNumber, numTests, time_coeff[1]);
					fprintf(filePtrLog, "coeff: %q\n", convert(time_coeff[2], string));
					fflush(filePtrLog);
					#also print to screen.  
					printf("Dim %d: test %d out of %d: time %f\n", currentDim, currentTestNumber, numTests, time_coeff[1]);
					print(time_coeff);
		
					#save the time.
					totalTime:=totalTime + time_coeff[1];
				catch:
					printf("Something went wrong: %q\n",lastexception);
					failFlag:=1; #try again.
				end try;
			end; #while the current test case did not fail.				
		od; #end for. for every test.
		
		#now, we just did numTest many test cases. so find the average and print it.
		#print to log.
		fprintf(filePtrLog, "Average %f by %d tests:\t %f\n",totalTime, numTests, totalTime/numTests);
		fflush(filePtrLog);
		#print to table 
		fprintf(filePtrTable, "%d\t%d\t%f\n", numTests, currentDim, totalTime/numTests);
		fflush(filePtrTable);
		
		currentDim:=currentDim + 1;
	end; #while. Keep trying higher dim. simplex.\
	
	fclose(filePtrLog);
	fclose(filePtrTable);
end:


#Find the average time to compute the top 3 coeff. of many simplices and increasing dim.
				#starting dim, number of tests, file base name.
table_time_sl_ehrhart(2, 10, "tableTimeSLTopEhrhart_oneTest");
#13, 50,

#Find the top 3 coeff. of 1 simplex.
#test_sl_ehrhart(6, "debug");



quit;
Profile(primitive_vector);
Profile(short_vector);
Profile(sign_entries_vector);
Profile(good_vector);
Profile(signed_decomp);
Profile(good_cone_dec); 
Profile(more_decomposition_in_cones);
Profile(cone_dec);
Profile(projectedvector);
Profile(projectedlattice);
Profile(projectedconeinbasislattice);
Profile(Toddzero);
Profile(relativevolumeoffaceiotac);
Profile(functionIzero);
Profile(prod_Toddzero);
Profile(functionSzero);
Profile(changeofcoordinates);
Profile(coeff_minusdplusk_iota_function_S);
Profile(coeff_minusdplusk_S);
Profile(coeff_dminusk_Eh_with_reg);
Profile(random_vector);
Profile(coeff_dminusk_Eh);



PrintProfiles(primitive_vector);
PrintProfiles(short_vector);
PrintProfiles(sign_entries_vector);
PrintProfiles(good_vector);
PrintProfiles(signed_decomp);
PrintProfiles(good_cone_dec);
PrintProfiles(more_decomposition_in_cones);
PrintProfiles(cone_dec);
PrintProfiles(projectedvector);
PrintProfiles(projectedlattice);
PrintProfiles(projectedconeinbasislattice);
PrintProfiles(projectedconeinbasislattice);
PrintProfiles(Toddzero);
PrintProfiles(relativevolumeoffaceiotac);

PrintProfiles(changeofcoordinates);
PrintProfiles(coeff_minusdplusk_iota_function_S);
