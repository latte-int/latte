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


