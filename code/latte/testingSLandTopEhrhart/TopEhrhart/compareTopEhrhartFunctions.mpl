read("testingTopEhrhart_lib.mpl"); #load the test functions
#
#  Oct 7, 2010.
#  Author: Brandon, Gregory
#  Description: 
#


for dim from 3 to 10 do:
	for numTests from 1 to 3 do:
		test_top_ehrhart_compare_v1_v2(dim, "compareV1V2");
	end;
end;
 

