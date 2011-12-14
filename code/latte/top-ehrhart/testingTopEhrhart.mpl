read("testingTopEhrhart_lib.mpl"); #load the test functions
interface(quiet=true);
#with(CodeTools[Profiling]):
#
#  Sept. 28, 2010
#  Author: Brandon, Gregory, Jesus
#  Description: Uses the main TopEhrhart lib to compute the top-three ehrhart coeff.
#	This script can either 1)generate a table of time test for finding the top ehrhart coeff. or 2) test
#

#Saves the average time it takes to find the top-3 ehrhart coeff in a file.
#For each random simplex of dim between startingDim to infinity, test numTests many ehrhart top 3 coefficients.
# The time of each simplex coeff. is saved in the log file, and the average of the numTests tests is saved in the table file.
# If the first test takes more than 1/2 hour, we stop the function.
# So the basic idea is:
#
#  for i from startindDim to infinity
#	find the time to do 1 test
#	if the test took longer than 1/2 hour, stop script.
#	else find the time for numTests-1 more tests and average them.
#  end for.
#Input:
#@parm: startingDim: starting dim. of the simplex.
#@parm: numTests: how many simpleices you want to test at once
#@parm: baseFileName: string. File names used for saving the time, average, and latte's facet equations and latte's volume answer. ex:"testingSL/testingSL_volume"
table_time_top_ehrhart:=proc(startingDim, numTests, fileBaseName, test_procedure)
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
        time_coeff := [0.0, ERROR];
		while failFlag <> 0 do
			try
				#time_coeff:=test_top_ehrhart_v1(currentDim, fileBaseName||".debug");
				#time_coeff:=test_top_ehrhart_v2(currentDim, fileBaseName||".debug");
				#time_coeff:=test_top_ehrhart_v3(currentDim, fileBaseName||".debug");
                time_coeff:=test_procedure(currentDim, fileBaseName||".debug");
                print(time_coeff):

				failFlag:=0; #no erros. break out of while loop.
			catch:
				printf("Something went wrong: %q\n",lastexception);
				failFlag:=1; #try again.
			end try;
		od; #while the first test case did not fail.

		#print to log.
		fprintf(filePtrLog, "Dim %d:\t test %d out of %d:\t time %f\n", currentDim, currentTestNumber, numTests, time_coeff[1]);
		fprintf(filePtrLog, "coeff: %q\n", convert(time_coeff[2], string));
		fflush(filePtrLog);
		#also print to screen.
		printf("Dim %d: test %d out of %d: time %f\n", currentDim, currentTestNumber, numTests, time_coeff[1]);
		print(time_coeff);

		#finally, how long did 1 test take? If more than 30min*60sec/min, stop.
		if time_coeff[1] > 60*30 then:
			printf("Test took longer than 1/2 hour\n");
			#currentDim:=currentDim + 1;
			#next;
			break; #break the while 1 = 1 loop.
		fi;

		#save the time.
		totalTime:=time_coeff[1];

		#get times for numTest -1 more examples!!!
		for currentTestNumber from 2 to numTests do:
			failFlag:= 1;
			while failFlag <> 0 do
				try
					#time_coeff:=test_top_ehrhart_v1(currentDim, fileBaseName||".debug");
					#time_coeff:=test_top_ehrhart_v2(currentDim, fileBaseName||".debug");
					#time_coeff:=test_top_ehrhart_v3(currentDim, fileBaseName||".debug");
                    time_coeff:=test_procedure(currentDim, fileBaseName||".debug");
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
			od; #while the current test case did not fail.
		od; #end for. for every test.

		#now, we just did numTest many test cases. so find the average and print it.
		#print to log.
		fprintf(filePtrLog, "Average %f by %d tests:\t %f\n",totalTime, numTests, totalTime/numTests);
		fflush(filePtrLog);
		#print to table
		fprintf(filePtrTable, "%d\t%d\t%f\n", numTests, currentDim, totalTime/numTests);
		fflush(filePtrTable);

		currentDim:=currentDim + 1;
	od; #while. Keep trying higher dim. simplex.

	fclose(filePtrLog);
	fclose(filePtrTable);
end:


#Find the average time to compute the top 3 coeff. of many simplices and increasing dim.
				#starting dim, number of tests, file base name.
table_time_top_ehrhart(3, 50, "tableTimeSLTopEhrhart-3", test_top_ehrhart_v2);


#test_top_ehrhart_compare_v1_v2(4, "compareV1V2");


#Find the top 3 coeff. of 1 simplex.
#test_top_ehrhart(6, "oneTopEhrhartTest");



quit;
##resetprofile();
profile(primitive_vector);
profile(short_vector);
profile(sign_entries_vector);
profile(good_vector);
profile(signed_decomp);
profile(good_cone_dec);
profile(more_decomposition_in_cones);
profile(cone_dec);
profile(projectedvector);
profile(projectedlattice);
profile(projectedconeinbasislattice);
profile(Toddzero);
profile(relativevolumeoffaceiotac);
profile(functionIzero);
profile(prod_Toddzero);
profile(functionSzero);
profile(coeff_minusdplusk_S);
profile(coeff_dminusk_Eh_with_reg);
profile(random_vector);
profile(coeff_dminusk_Eh);

profile(changeofcoordinates);
profile(coeff_minusdplusk_iota_function_S);

profile(projectedvertexinbasislattice);





# All procedures in Conebyconeapproximations...:
profile(insert);
profile(ComplementList);
profile(GeneralComplementList);
profile(special_lincomb_v);
profile(primitive_vector);
profile(ortho_basis);
profile(fracpart);
profile(short_vector);
profile(sign_entries_vector);
profile(good_vector);
profile(signed_decomp);
profile(good_cone_dec);
profile(more_decomposition_in_cones);
profile(cone_dec);
profile(projectedvector);
profile(projectedvector_with_inverse);
profile(projectedlattice);
profile(projectedconeinbasislattice);
profile(projectedvertexinbasislattice);
profile(s_ISpace);
profile(Todd);
profile(ourceil);
profile(fractionalpart);
profile(fmod);
profile(ourmod);
profile(nfractionalpart);
profile(ourmodreal);
profile(nfractionalpartreal);
profile(volume_ISpace);
profile(functionIa);
profile(functionIb);
profile(prod_Todd);
profile(functionS);
profile(changeofcoordinates);
profile(S_Ispace_Coneformulaa);
profile(S_Ispace_Coneformulab);
profile(linindenom);
profile(approx_Cone_formulaa);
profile(approx_Cone_formulab);
profile(cone_by_cone_approxi_simplex_formulaa);
profile(cone_by_cone_approxi_simplex_formulab);
profile(dilatedS_Ispace_Cone);
profile(dilated_approxi_cone);
profile(ApproxEhrhartSimplexgeneric);
profile(random_vector);
profile(TopEhrhartweightedluckyell);
profile(TopEhrhartweighted);
profile(CompleteEhrhartweighted);
profile(TopEhrhartweightedPoly);
profile(printTopEhrhartweightedPoly);
profile(printIncrementalEhrhartweightedPoly);
profile(dilatedS_Ispace_Cone_real);
profile(dilated_approxi_cone_real);
profile(ApproxEhrhartSimplexgeneric_real);
profile(random_vector);
profile(TopEhrhartweightedluckyell_real);
profile(TopEhrhartweighted_real);
profile(CompleteEhrhartweighted_real);
profile(TopEhrhartweightedPoly_real);
profile(printTopEhrhartweightedPoly_real);
profile(printIncrementalEhrhartweightedPoly_real);
profile(random_rational_vector);
profile(randomaffinecone);
profile(random_rational_simplex);
profile(checkapprox);


#TopEhrhartweightedPoly(n, [[-9, 40, -73], [-37, -46, -7], [70, -85, 31], [59, 99, 68]], [seq(0,i=1..3)], 0, 2);
#TopEhrhartweightedPoly(n, [[-9, 40, -73], [-37, -46, -7], [70, -85, 31], [59, 99, 68]], [seq(0,i=1..3)], 0, 2);
#table_time_top_ehrhart(4, 3, "tableTimeSLTopEhrhart-v3-4", test_top_ehrhart_v3);

Simplex := [[-49, -89, -14, -1], [77, -46, 26, 6], [44, 97, -85, 93],
    [94, -7, -47, -72], [30, 32, 94, 19]];
#TopEhrhartweightedPoly(n, Simplex, [seq(0,i=1..4)], 0, 2);
#Topk_Eh(Simplex,2,t);
Simplex6 := [[56, -47, 27, 95, 79, -30], [88, -16, 18, -13, 77, -86], [-52, -89, 67, 21, 64, 22], [-31, 71, 78, 81, -90, 64], [-37, -31, -33, 82, 17, -22], [48, 29, -19, 72, -79, -16], [45, -70, -41, -51, -91, -66]];
TopEhrhartweightedPoly(n, Simplex6 , [seq(0,i=1..6)], 0, 2);
#Topk_Eh(Simplex6,2,t);

showprofile();
resetprofile();

