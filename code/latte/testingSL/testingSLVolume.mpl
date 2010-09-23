read("testingSL/testingSL_lib.mpl");	#load the testing functions.
with(CodeTools[Profiling]):
#
#  Sept. 9, 2010
#  Author: Brandon
#  Description: Uses the main SL testing lib to compute the volume of simplies and compare the answer's with latte's tirangulation method.
#

#Input:
#@parm: simplexDim: the abm. dim of the simplix
#@parm: numTests: how many simpleices you want to test at once
#@parm: baseFileName: string. File names used for saving latte's facet equations and latte's volume answer. ex:"testingSL/testingSL_volume"
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


table_time_sl_volume:=proc(startDim, numTests, fileBaseName)
local fileName, filePtr, i, totalTime;
	
	fileName:=fileBaseName||"2_"||startDim||"_"||numTests||".txt";
	filePtr:=fopen(fileName, WRITE, TEXT);
	fprintf(filePtr, "Testing SL volume on simplex of dim %d and up with %d samples\n", startDim, numTests); 

	fprintf(filePtr, "Simplex Dim | time of first case | average time\n");
	for i from startDim to 10000 do:
		fprintf(filePtr, "Dim %d ", i);
		print("Dim %d ", i);
		totalTime:=time_sl_volume(i, 1);#print("hot here", temp);
		fprintf(filePtr, ",%f", totalTime);
		fflush(filePtr);
		
		if totalTime < 60*60 and 1 < numTests then:
			timeRest:=time_sl_volume(i, numTests - 1);
			totalTime:=totalTime + timeRest; 
			fprintf(filePtr, ",%f", totalTime/numTests);
		fi; #if under an hour, do the rest of the tests.
		fprintf(filePtr, " Done.\n");
		fflush(filePtr);
	od;
	fclose(filePtr);
end:

#table_time_sl_volume(6, 3, "testingSL/tableTimeSLVolume");


test_sl_volume(6, 10, "testingSL/testingSL_volume");



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


