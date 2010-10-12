read("testingTopEhrhart_lib.mpl"); #load the test functions
read("BirkhoffGenerator.mpl"); #read the Birkhoff Simplex Generator Code.



#
#  Oct 12, 2010
#  Author: Brandon, Gregory
#  Description: Uses the main TopEhrhart lib to compute the top-three ehrhart coeff of the Birkhoff Simplex.
#	
#

print("computing the birkhoff simplex:");
#Zeros:=[[1,3],[1,4],[1,5],[2,4],[2,5],[3,5]];
#birkhoff:=gensimplexB_n(2,5,Zeros):

		        #n, k , F
birkhoff:=gensimplexB_n(3,2,[]);
birkhoff:=convert(birkhoff, listlist);

#mysimplex:=create_random_simplex(3);

test_top_ehrhart_given_simplex_v2(birkhoff);






