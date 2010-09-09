

with(linalg):with(LinearAlgebra):
read("integration/createLinear.mpl"): #load the function to make a random linear form
read("testingSL/SL_lib.mpl");	#load the SL functions.


#insert making  simplex H-rep here.


test_sl_volume:=proc(simplexDim, numTests, fileBaseName)
local simplexList, i, volumeAnswersLatte, fileNameSimplex, fileNameVolume, systemCommand, status:

	fileNameSimplex:=baseFileName || ".simplex" :
	fileNameVolume:=baseFileName || ".volume":

	simplexList:= [];
	for i from 1 to numTests do:
		simplexList:=[[op(make_random_simplex(simplexDim))], op(simplexList)];
	od;

	write_simplex_to_file(simplexList, fileNameSimplex);

	#Ask C++ exe to compute the volumes of the simplex.
	systemCommand:= "./testVolumeForSL " || fileNameSimplex  || " " || fileNameVolume :
	#print(systemCommand);
	#status:=system(systemCommand):
	#printf("status=%d\n", status);


	volumeAnswersLatte:=read_volume_from_file(fileNameVolume, numTests);
end:


read_volume_from_file:=proc(fileName, numberAnswers)
	local filePtr, volumeList, i;

	volumeList:=[];
	filePtr:=fopen(fileName, READ, TEXT);
	for i from 1 to numberAnswers do:
		volumeList:=[ readline(filePtr), volumeList];
	od:

	close(filePtr);
	print(volumeList, "=read_volume.volumeList");
end:

write_simplex_to_file:=proc(simplexList, fileName)
	local filePtr:
	filePtr:=fopen(fileName,WRITE,TEXT):
	
	for i from 1 to nops(simplexList) do:	
		writeline(filePtr, convert(simplexList[i], string));
	od;	
	close(filePtr);
end:




#tfunction_SL_simplex_ell(1,[[0, 0], [1/7*5/4, 0], [1/7*5/4, 1/7*5/2]],[[1,0]],[1,1],1); 
	

test_sl_volume(5, 3, "testingSL_volume");

