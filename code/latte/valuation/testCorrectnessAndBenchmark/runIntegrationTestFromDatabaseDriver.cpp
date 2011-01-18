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
#include "../../sqlite/IntegrationDB.h"



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
		string polynomial  = "--monomials="+(*row)[0];
		string latte       = (*row)[1];
		string lawrence    = (*row)[2];
		string triangulate = (*row)[3];
		string value       = (*row)[4];
		string rowid       = (*row)[5];

		if ( value != "NA")
		{
			if ( alg == 1 && atof(lawrence.c_str()) >= 0)
			{
				cout << "Skipping lawrence integration on " << latte.c_str() << " " << polynomial.c_str() << endl;
				continue;
			}
			if ( alg == 2 && atof(triangulate.c_str()) >= 0)
			{
				cout << "Skipping triangulation integration on " << latte.c_str() << " " << polynomial.c_str() << endl;
				continue;
			}
			if (atof(triangulate.c_str()) >= 0 && atof(lawrence.c_str()) >= 0)
			{
				cout << "Skipping integration on " << latte.c_str() << " " << polynomial.c_str() << endl;
				continue;
			}
		}//check to see if the test is already done or not

		const char *argv[6];
		argv[0] = "./runTheTests";
		argv[1] =  "--valuation=integrate";
		if (alg == 1) //if you cange this, allso need to chage the other alg get below.
			argv[2] = "--lawrence";
		else if (alg == 2)
			argv[2] = "--triangulate";
		else
		{
			cout << "runTheTests::unknown alg." << endl;
			exit(1);
		}
		argv[3] = polynomial.c_str();
		argv[4] = "--vrep";
		argv[5] = latte.c_str();

		//open the log file and print what we are about to do.
		ofstream log(logFileName.c_str(), ios::app);
		time_t rawtime;
		struct tm * timeinfo;
		time ( &rawtime );
		timeinfo = localtime ( &rawtime );

		log << "\n" << asctime (timeinfo) << ": Testing " << (alg == 1? "lawrence" : "triangulate") << endl;
		log << "\tLatte: " << latte.c_str()
			<< "\n\tPolynomial: " << polynomial.c_str()
			<< "\n\trowid: " << rowid.c_str();
		log.close();
		//RUN IT. finally.
		Valuation::ValuationContainer vc = Valuation::mainValuationDriver(argv, 6);

		//ok, now save the results.

		RationalNTL theComputedIntegral;
		float theTotalTime;
		for( int i = 0; i < vc.answers.size(); ++i)
		{
			if ( vc.answers[i].valuationType == Valuation::ValuationData::integrateLawrence
					&& alg == 1)
			{
				theComputedIntegral = vc.answers[i].answer;
			}
			else if (vc.answers[i].valuationType == Valuation::ValuationData::integrateTriangulation
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
		IntegrationDB db;
		db.open(dbFile);
		db.updateIntegrationTimeAndValue((alg == 1 ? IntegrationDB::Lawrence : IntegrationDB::Triangulate)
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
}


int main(int argc, char *argv[])
{
	//if (argc != 2)
	//{
	//	cout << "error. usage: " << argv[0] << "sqlite-file-name" << endl;
	//	exit(1);
	//}
	if( argc != 7 )
	{
		cout << "error: usage: " << argv[0] << " database-file dim vertex degree dual[true|false] algorithm[1=lawrence,2=triangulate] " << endl;
		exit(1);
	}

//	char * dbFile, int dim, int vertex, int degree, bool useDual, int alg, int limit)
	runIntegrationTest(argv[1], atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), string(argv[5]) == "true", atoi(argv[6]), 50, string(argv[0])+".log");
	return 0;
}//main
