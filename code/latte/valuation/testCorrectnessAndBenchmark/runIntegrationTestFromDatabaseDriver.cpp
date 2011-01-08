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
				continue;
			if ( alg == 2 && atof(triangulate.c_str()) >= 0)
				continue;
			if (atof(triangulate.c_str()) >= 0 && atof(lawrence.c_str()) >= 0)
				continue;
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

		log << asctime (timeinfo) << ": Testing " << (alg == 1? "lawrence" : "triangulate") << endl;
		log << "\tLatte: " << latte.c_str()
			<< "\n\tPolynomial: " << polynomial.c_str()
			<< "\n\trowid: " << rowid.c_str();
		log.close();
		//RUN IT. finally.
		Valuation::ValuationContainer vc = Valuation::mainValuationDriver(argv, 6);

		//ok, now save the results.
		IntegrationDB db;
		db.open(dbFile);
		for( int i = 0; i < vc.answers.size(); ++i)
		{
			stringstream sql;
			if ( vc.answers[i].valuationType != Valuation::ValuationData::entireValuation)
				continue;

			db.updateIntegrationTimeAndValue((alg == 1 ? IntegrationDB::Lawrence : IntegrationDB::Triangulate)
					, vc.answers[i].timer.get_seconds(), vc.answers[i].answer, value, rowid);
		}
		db.close();
	}

}//runTheTests



void runIntegrationTest(char * dbFile, int dim, int vertex, int degree, bool useDual, int alg, int limit, string log)
{
	vector<vector<string> > toTest;
	IntegrationDB db;
	db.open(dbFile);
	toTest = db.getRowsToIntegrate(dim, vertex, degree, useDual, limit);
	db.close();

	if (! toTest.size() )
	{
		cout << "test case is empty" << endl;
		exit(1);
	}
	runTheTests(toTest, alg, dbFile, log);
}

/*
void runIntegrationTest(char * dbFile)
{
	stringstream sql;
	string line;
	sql << "select p.filePath, t.latteFilePath, i.timeLawrence, i.timeTriangulate, i.integral, i.rowid"
		<< "\n from polynomial as p, polytope as t, integrate as i "
		<< "\n where p.rowid = i.polynomialID and t.rowid = i.polytopeID ";
		//and p.degree = ?
		//and t.dim = ?
		//and t.vertexCount = ?
		//and t.simple = ?
	    //and t.dual =?
		//order by t.latteFilePath, p.degree
	cout << "Please finish the sql selection statement" << endl;
	cout << sql.str().c_str() << endl;

	cout << " p.degree [>|<|=|<=,etc] [number] ";
	getline(cin, line);
	sql << " and p.degree " << line << " ";

	cout << " t.dim [>|<|=|<=,etc] [number] ";
	getline(cin, line);
	sql << " and t.dim " << line << " ";

	cout << " t.vertexCount [>|<|=|<=,etc] [number] ";
	getline(cin, line);
	sql << " and t.vertexCount " << line << " ";

	cout << " t.simple [>|<|=|<=,etc] [0|1] ";
	getline(cin, line);
	sql << " and t.simple " << line << " ";

	cout << " t.dual [is|is not] null :";
	getline(cin, line);
	sql << " and t.dual " << line << " null ";

	cout << " order by [t.latteFilePath, p.degree] ";
	getline(cin, line);
	sql << " order by " << line;

	vector<vector<string> > toTest;
	IntegrationDB db;
	db.open(dbFile);
	toTest = db.query(sql.str().c_str());
	db.close();

	if (! toTest.size() )
	{
		cout << "test case is empty" << endl;
		exit(1);
	}

	cout << " Use lawrence(1) or triangulation(2) or both(3)? ";
	int alg;
	cin >> alg;

	runTheTests(toTest, alg, dbFile);
}
*/

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
