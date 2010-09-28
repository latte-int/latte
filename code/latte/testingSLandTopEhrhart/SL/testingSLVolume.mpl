read("testingSL_lib.mpl");	#load the testing functions.
with(CodeTools[Profiling]):
#
#  Sept. 9, 2010
#  Author: Brandon
#  Description: Uses the main SL testing lib to compute the volume of simplies and compare the answer's with latte's tirangulation method.
#



#Input:
#@parm: simplexDim: the abm. dim of the simplix
#@parm: numTests: how many simpleices you want to test at once
#@parm: baseFileName: string. File names used for saving latte's facet equations, latte's volume answer, the log file, and the average file. 
#Finds the average time of finding the volume of a random simplex.
#If the time to find the volume is larger than 1/2 hour, then only 1 test is done instead of numTests many.
table_time_sl_volume:=proc(startDim, numTests, fileBaseName)
local fileNameLog, filenameTable, filePtrLog, filePtrTable, currentTest, i, totalTime, timeRest;
	
	fileNameLog:=fileBaseName||".log";
	fileNameTable:=fileBaseName||".average";
	filePtrLog:=fopen(fileNameLog, WRITE, TEXT);
	filePtrTable:=fopen(fileNameTable, WRITE, TEXT);
	fprintf(filePtrLog, "Testing SL volume on simplex of dim %d and up with %d samples\n", startDim, numTests); 

	fprintf(filePtrTable, "Number of tests | Simplex Dim | average time\n");
	for i from startDim to 10000 do:
		totalTime:=time_sl_volume(i, 1); #for now, test only 1 at a time.
		fprintf(filePtrLog, "Dim %d, test 1 out of %d, time %f\n", i, numTests, totalTime);
		fflush(filePtrLog);
		
		if totalTime < 60*60 and 1 < numTests then:
			for currentTest from 2 to numTests do:
				timeRest:=time_sl_volume(i, 1);
				totalTime:=totalTime + timeRest;
				fprintf(filePtrLog, "Dim %d, test %d out of %d, time %f\n", i, currentTest, numTests, totalTime);
				printf("Dim %d, test %d out of %d, time %f\n", i, currentTest, numTests, totalTime);
			end; 
			fprintf(filePtrTable, "%d\t %d\t %f\n", i, numTests, totalTime/numTests);
			fflush(filePtrTable);
		fi; #if under an hour, do the rest of the tests.		
	od;
	fclose(filePtrTable);
	fclose(filePtrLog);
end:

table_time_sl_volume(6, 3, "tableTimeSLVolume");


#test_sl_volume(6, 10, "testingSL_volumeCorrectness");

quit;
Profile(denomWL);
Profile(functionI);
Profile(linindenom);
Profile(projectedconeinbasislattice);
Profile(projectedlattice);
Profile(projectedvector); 
Profile(projectedvertexinbasislattice);
Profile(regularSL);
Profile(tfunction_SL);
Profile(ttruncatedSL);
Profile(tSLell);


PrintProfiles(denomWL);
PrintProfiles(functionI);
PrintProfiles(linindenom);
PrintProfiles(projectedconeinbasislattice);
PrintProfiles(projectedlattice);
PrintProfiles(projectedvector);
PrintProfiles(projectedvertexinbasislattice);
PrintProfiles(regularSL);
PrintProfiles(tfunction_SL);
PrintProfiles(ttruncatedSL);
PrintProfiles(tSLell);


