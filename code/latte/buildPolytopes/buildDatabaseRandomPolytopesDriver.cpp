/*
 * buildDatabaseRandomPolytopesDriver.cpp
 *
 *  Created on: Nov 24, 2010
 *      Author: Brandon Dutra and Gregory Pinto
 *
 *      Computes H-rep latte files.
 */
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include "BuildRandomPolytope.h"


using namespace std;

/**
 * Builds a random polytope with the name "rp_dim_currentPolytopeCount in the directory fileBaseName
 * We add a slash after the fileBaseName/director/filepath.
 */

void buildDataBaseOfPolytopes(const string & fileBaseName,
		                      const int dim,
		                      const int numberTestCases)
{
	stringstream polytopeBaseFileName;
	int currentPolytopeCount = 1;

	while( currentPolytopeCount <= numberTestCases)
	{
		polytopeBaseFileName.str("");
		cout << "Going to make " << currentPolytopeCount << "th polytope for dim " << dim << " in location " << fileBaseName.c_str() << endl;
		BuildRandomPolytope rp; //random polytope.

		polytopeBaseFileName << fileBaseName << "/rp_" << dim << "_" << currentPolytopeCount;
		rp.setBaseFileName(polytopeBaseFileName.str());
		rp.makePoints(dim, dim*2+1);
		rp.buildPolymakeFile();
		rp.buildLatteVRepFile();

		if ( rp.getDim() != dim)
		{
			cout << "darn, rp dim=" << rp.getDim() << endl;
			cout << "files " << rp.getLatteVRepFile().c_str() << ", " << rp.getPolymakeFile().c_str() << endl;
			continue; //start over with this test case.
		}


		//rp.deletePolymakeFile(); //maybe we should also keep the polymake file in case you just want a V-rep or for some other reason.
		++currentPolytopeCount;
	}//while
}//buildDataBaseOfPolytopes



int main(int argc, char *argv[])
{
	//these 3 arrays should be parallel.
	vector<string> fileBaseNames;
	vector<int> dim;
	vector<int> numberTestCases;
	string currentDir;
	int currentDim;
	int currentNumberTestCase;
	int caseNumber = 1;


	cout << "Give me a list of polytope classes you want to build.\n"
		 << "Example: build class 1 > dir2 2 10\n"
		 << "         build class 2 > dir3 3 12\n"
		 << "         build class 3 > done\n"
		 << "This will find 10 random polytopes of dim 2 and save it in dir2, and save 12 random polytopes of dim 3 in dir3."
		 << " The directory should already exist." << endl;


	cout << "build class 1 > ";
	cin >> currentDir;
	while ( currentDir != "done")
	{
		cin >> currentDim >> currentNumberTestCase;

		fileBaseNames.push_back(currentDir);
		dim.push_back(currentDim);
		numberTestCases.push_back(currentNumberTestCase);

		cout << "build class " << ++caseNumber << " > ";
		cin >> currentDir;
	}//while dir.

	//now that we are done asking the user for what they want to build and where, do it.
	for(int i = 0; i < (int) fileBaseNames.size(); ++i)
		buildDataBaseOfPolytopes(fileBaseNames[i], dim[i], numberTestCases[i]);




	return 0;
}//main
