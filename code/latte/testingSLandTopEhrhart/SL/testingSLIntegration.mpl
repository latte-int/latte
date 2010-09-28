read("testingSL_lib.mpl");	#load the testing functions.
with(CodeTools[Profiling]):
#
#  Sept. 9, 2010
#  Author: Brandon
#  Description: Uses the main SL testing lib to compute the volume of simplies and compare the answer's with latte's tirangulation method.
#

lattice_random_simplex:=proc(d,N) local R,U;
  R := rand(N):
  U:=proc()[seq(R(),i=1..d)] end proc:
  [ seq(U(), i=1..d+1) ];
end:


#Input:
#@parm: simplexDim: the abm. dim of the simplix
#@parm: numTests: how many simpleices you want to test at once
#@parm: degreeL: The power you want to take the linear form to.
test_sl_integration(5, 5, 3);

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
