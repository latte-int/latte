/*
 * runVolumeTestFromDatabaseDriver.cpp
 *
 *  Created on: Jan 25, 2011
 *      Author: bedutra
 */

/*
 * runIntegrationTestFromDatabaseDriver.cpp
 *
 *  Created on: Jan 5, 2011
 *      Author: bedutra
 */

#include <string>
#include <cstdio>
#include <vector>
#include <sstream>
#include <iostream>
#include <ctime>
#include "../valuation.h"
#include "../../sqlite/VolumeDB.h"



/**
 * Integrate each row in the toTest vector.
 * Makes a log file of the current row being tested.
 *
 * @alg. 1=lawrence, 2=triangulate
 * @dbFile* : sqlite database file with tables set up as in the create script
 * @toTest: sql query listing what rows in the integrate table we need to run.
 * @logFileName : print the current test info to the log file.
 */
void runTheTests(const vector<vector<string> > &toTest, int alg, const char*dbFile, string logFileName)
{
	ofstream log;
	for(vector<vector<string> >::const_iterator row = toTest.begin(); row != toTest.end(); ++row)
	{
		string latte       =  "."+(*row)[0];//look at ../filePath (assumes latte is in the form ./filePath).
		string lawrence    = (*row)[1];
		string triangulate = (*row)[2];
		string value       = (*row)[3];
		string rowid       = (*row)[4];

		if ( value != "NA")
		{
			if ( alg == 1 && atof(lawrence.c_str()) >= 0)
			{
				cout << "runTheTests::done with lawrence volume on " << latte.c_str() << endl;
				continue;
			}
			if ( alg == 2 && atof(triangulate.c_str()) >= 0)
			{
				cout << "runTheTests::done with triangulation volume on " << latte.c_str() << endl;
				continue;
			}
			if (atof(triangulate.c_str()) >= 0 && atof(lawrence.c_str()) >= 0)
			{
				cout << "runTheTests::done with volume on " << latte.c_str() << endl;
				continue;
			}
		}//check to see if the test is already done or not

		const char *argv[7];
		int index = 0;
		argv[index++] = "./runTheTests";
		argv[index++] = "--redundancy-check=none";
		argv[index++] =  "--valuation=volume";
		if (alg == 1) //if you change this, also need to change the other alg usage below.
			argv[index++] = "--lawrence";
		else if (alg == 2)
			argv[index++] = "--triangulate";
		else
		{
			cout << "runTheTests::unknown alg." << endl;
			exit(1);
		}
		if ( string::npos != latte.rfind(".vrep."))
			argv[index++] = "--vrep";
		else if (string::npos == latte.rfind(".hrep."))
		{
			cout << "unknown latte file data type: " << latte.c_str() << endl;
			exit(1);
		}//if not a vrep or hrep file, error.
		argv[index++] = latte.c_str();

		//open the log file and print what we are about to do.
		ofstream log(logFileName.c_str(), ios::app);
		time_t rawtime;
		struct tm * timeinfo;
		time ( &rawtime );
		timeinfo = localtime ( &rawtime );

		log << "\n" << asctime (timeinfo) << ": Testing rowid " << rowid.c_str();
		for(int i = 0; i < index; ++i)
			log << "\n\t" << argv[i];
		log.close();
		//RUN IT. finally.
		Valuation::ValuationContainer vc = Valuation::mainValuationDriver(argv, index);

		//ok, now save the results.

		RationalNTL theComputedIntegral;
		float theTotalTime;
		for( int i = 0; i < vc.answers.size(); ++i)
		{
			if ( vc.answers[i].valuationType == Valuation::ValuationData::volumeLawrence
					&& alg == 1)
			{
				theComputedIntegral = vc.answers[i].answer;
			}
			else if (vc.answers[i].valuationType == Valuation::ValuationData::volumeTriangulation
					&& alg == 2)
			{
				theComputedIntegral = vc.answers[i].answer;
			}

			if ( vc.answers[i].valuationType == Valuation::ValuationData::entireValuation)
			{
				theTotalTime = vc.answers[i].timer.get_seconds();
			}

		}//find the integral and time

		//print results to the log file
		log.open(logFileName.c_str(), ios::app);
		log << "\n\tvaluation type: " << alg
			<< "\n\tvaluation time: " << theTotalTime
			<< "\n\tvaluation ans: "  << theComputedIntegral;
		log.close();

		//finally, save it.
		VolumeDB db;
		db.open(dbFile);
		db.updateVolumeTimeAndValue((alg == 1 ? VolumeDB::Lawrence : VolumeDB::Triangulate)
							, theTotalTime, theComputedIntegral, value, rowid);
		db.close();

		if ( theTotalTime > 660 ) //660 sec= 10.5mins
		{
			cout << "Valuation took " << theTotalTime << " seconds, skipping the test of the test case" << endl;
			return;
		}//if it took too long, don't do the rest.
	}

}//runTheTests



void runIntegrationTest(char * dbFile, int dim, int vertex, int degree, bool useDual, int alg, int limit, string log)
{
	/*
	vector<vector<string> > toTest;
	IntegrationDB db;
	db.open(dbFile);
	if ( limit <= db.testCasesCompleted((alg == 2 ? IntegrationDB::Triangulate : IntegrationDB::Lawrence),dim, vertex, degree, useDual) )
	{
		cout << "Skipping " << dim << " " << vertex << " " << degree << " " << useDual << " " << alg << ":: test class already done" << endl;
		db.close();
		exit(1);
	}//if test already done
	if( db.canTestFinish((alg == 2 ? IntegrationDB::Triangulate : IntegrationDB::Lawrence),dim, vertex, degree, useDual, 600) )
	{
		toTest = db.getRowsToIntegrate(dim, vertex, degree, useDual, limit);
	}//if we think the test can finish in 600 sec, then do it.
	db.close();

	if (! toTest.size() )
	{
		cout << "test case is empty" << endl;
		db.close();
		exit(1);
	}
	runTheTests(toTest, alg, dbFile, log);
	*/
}


//db file, file name, use dual, alg, log
void runSpecficPolytopeTest(char * dbFile, char * polymakeFile, bool useDual, int alg,string log)
{
	vector<vector<string> > toTest;
	VolumeDB db;

	db.open(dbFile);

	toTest = db.getRowsToFindVolumeGivenSpecficFile(polymakeFile, useDual);

	assert(toTest.size() == 1);
	if ( (atof(toTest[0][1].c_str()) == -2.0 && alg == 1)
		|| (atof(toTest[0][2].c_str()) == -2.0 && alg == 2) )
	{
		cout << polymakeFile << " cannot finish" << endl;
		db.close();
		return;
	}

	//}
	//else
	//{
	//	cout << "runSpecficPolytopeTest:: skipping " << polymakeFile << " " << degree << " " << useDual << " " << alg << endl;
	//}
	db.close();

	runTheTests(toTest, alg, dbFile, log);

}



int main(int argc, char *argv[])
{


	if( argc != 6 )
	{
		cout << "error: usage: " << argv[0] << " database-file dim vertex dual[true|false] algorithm[1=lawrence,2=triangulate] " << endl;
		cout << "error: usage: " << argv[0] << " database-file specficFile fileName dual[true|fale] algorithm[1|2]" << endl;

		cout << "NOTE::we assume this file is one directory BELOW the polytope folders. If the database says there is a polytope in ./a/b, we will look at ../a/b for the file"
			 << "This is done because I want to run volume and integration tests at the same time using the same database." << endl;

		exit(1);
	}

	//	char * dbFile, int dim, int vertex, int degree, bool useDual, int alg, int limit)
	if (argc == 6 && strcmp(argv[2], "specficFile"))
	{
		cout << "sorry, volume for regular polytopes not finished yet." << endl;
		exit(1);
		//runIntegrationTest(argv[1], atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), string(argv[5]) == "true", atoi(argv[6]), 50, string(argv[0])+".log");
	}
	else if ( argc == 6 && !strcmp(argv[2], "specficFile") )
		                     //db file, file name, use dual, alg, log
		runSpecficPolytopeTest(argv[1], argv[3], string(argv[4]) == "true", atoi(argv[5]), string(argv[0])+".log");
	else
		cout << "unkown sequence of parameters" << endl;
	return 0;
}//main
