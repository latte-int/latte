/*
 * testVolumeForSLDriver.cpp
 *
 *  Created on: Sept 8, 2010
 *      Author: Brandon Dutra and Gregory Pinto
 *
 *  A maple script uses this exe to find the volume of many polytopes and output the number in an output file.
 *
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include "../valuation.h"


using namespace std;


int main(int argc, char *argv[])
{
	const char * arguments[6];
	ofstream ouput;
	int numberTests = 1;

	streambuf * cout_strbuf(cout.rdbuf()); //keep copy of the real cout.


	if (argc != 4)
	{
		cout << "Usage " << argv[0] << " inputLatteFile outputFile numTest" << endl;
		cout << "\n  numTests: > we will compute the volume "
			 << "of the polytope defined in files inputLatteFile1...inputLatteFile(numTests)" << endl;
		exit(1);
	}//if errors.

	if ( ! (numberTests = atoi(argv[argc-1])))
	{
		cout << "Sorry, " << argv[argc-1] << " is not an integer." << endl;
		exit(1);
	}//if not an inteter.


	ouput.open(argv[2]);
	if ( ! ouput.is_open() )
	{
		cout << "cannot open ouput file!!!" << endl;
		exit(1);
	}

	arguments[0] = (const char *) argv[0];
	arguments[1] = (const char *) "--valuation=volume";
	arguments[2] = (const char *) "--triangulate";
	arguments[3] = (const char *) "--redundancy-check=none";
	//arguments[4] = file name.


	for(int i = 1; i <= numberTests; ++i)
	{
		cout.rdbuf(cout_strbuf); //revert cout back to cout!
		cout << "Latte Valuation: going to compute volume of simplex " << i << " out of " << numberTests << endl;
		cout.rdbuf(cerr.rdbuf()); //change cout to cerr
		Valuation::ValuationContainer answer_times;

		stringstream ss;
		char fileName[100];
		ss << argv[1] << i;
		strcpy(fileName, ss.str().c_str());
		arguments[4] = (const char *)fileName;

		answer_times = Valuation::mainValuationDriver(arguments, 5);

		for(int k = 0; k < answer_times.answers.size(); ++k)
			if ( answer_times.answers[k].valuationType == Valuation::ValuationData::volumeTriangulation)//ValuationData
			{
				ouput << answer_times.answers[k].answer << endl;
				cout.rdbuf(cout_strbuf); //revert cout back to cout!
				cout << "Latte Volume: " << i << ": " << answer_times.answers[k].answer
				     << "\nLatte Time: " << answer_times.answers[k].timer << '\n' << endl;
				cout.rdbuf(cerr.rdbuf());
				break;
			}
	}//for every test file.
	cout.rdbuf(cout_strbuf); //revert cout back to cout!
	ouput.close();
	return 0;
}//main()


