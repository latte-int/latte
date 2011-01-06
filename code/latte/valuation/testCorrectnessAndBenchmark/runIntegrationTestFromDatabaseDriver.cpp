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
#include "../valuation.h"
#include "../../sqlite/IntegrationDB.h"



/**
 * @alg. 1=lawrence, 2=triangulate, 3 = both!
 */
void runTheTests(const vector<vector<string> > &toTest, int alg, const char*dbFile)
{
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
		if (alg == 1)
			argv[2] = "--lawrence";
		else if (alg == 2)
			argv[2] = "--triangulate";
		else
			argv[2] = "--all";
		argv[3] = polynomial.c_str();
		argv[4] = "--vrep";
		argv[5] = latte.c_str();

		//RUN IT. finally.
		Valuation::ValuationContainer vc = Valuation::mainValuationDriver(argv, 6);

		//ok, now save the results.
		IntegrationDB db;
		db.open(dbFile);
		for( int i = 0; i < vc.answers.size(); ++i)
		{
			stringstream sql;
			sql << "update integrate set ";
			if ( vc.answers[i].valuationType == Valuation::ValuationData::integrateLawrence && atof(lawrence.c_str()) < 0)
				sql << " timeLawrence = " << vc.answers[i].timer.get_seconds();
			else if (vc.answers[i].valuationType ==Valuation::ValuationData:: integrateTriangulation && atof(triangulate.c_str()) < 0)
				sql << " timeTriangulate = " << vc.answers[i].timer.get_seconds();
			else
				continue;

			if ( value == "NA")
			{
				sql << " , integral = ' " << vc.answers[i].answer << "' ";
			}
			else
			{
				RationalNTL previousValue(value), computedValue(vc.answers[i].answer);
				if ( previousValue != computedValue)
				{
					cout << "ERROR: the integrals differ." << endl;
					cout << "previousValue: " << previousValue << endl;
					cout << "computedValue: " << computedValue << endl;
					cout << "current sql stm:" << sql.str().c_str();
					cout << "rowid: " << rowid.c_str() << endl;
					exit(1);
				}
			}
			sql << " where rowid = " << rowid << endl;
			db.query(sql.str().c_str());
		}
		db.close();

	}

}//runTheTests



void runIntegrationTest(char * dbFile, int dim, int vertex, int degree, bool useDual, int alg, int limit)
{
	stringstream sql;
	if (useDual == false)
	{
		sql << "select p.filePath, t.latteFilePath, i.timeLawrence, i.timeTriangulate, i.integral, i.rowid"
			<< " from polynomial as p, polytope as t, integrate as i "
			<< " where p.rowid = i.polynomialID and t.rowid = i.polytopeID "
			<< " and p.degree = " << degree
			<< " and t.dim = " << dim
			<< " and t.vertexCount =" << vertex
			<< " and t.dual is null"
			<< " order by t.latteFilePath, p.degree";
	}
	else
	{
		sql << "select p.filePath, t.latteFilePath, i.timeLawrence, i.timeTriangulate, i.integral, i.rowid"
			<< " from polynomial as p, polytope as t, polytope as orgP, integrate as i "
			<< " where p.rowid = i.polynomialID and t.rowid = i.polytopeID "
			<< " and orgP.rowid = t.dual"
			<< " and p.degree = " << degree
			<< " and t.dim = " << dim
			<< " and orgP.vertexCount =" << vertex
			<< " and t.dual is not null"
			<< " order by t.latteFilePath, p.degree";
	}


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
	runTheTests(toTest, alg, dbFile);

}


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

int main(int argc, char *argv[])
{
	//if (argc != 2)
	//{
	//	cout << "error. usage: " << argv[0] << "sqlite-file-name" << endl;
	//	exit(1);
	//}
	if( argc != 7 )
	{
		cout << "erro: usage: " << argv[0] << " database dim vertex degree dual[true|false] algorithem[1=lawrence,2=triangulate] " << endl;
		exit(1);
	}

//	char * dbFile, int dim, int vertex, int degree, bool useDual, int alg, int limit)
	runIntegrationTest(argv[1], atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), string(argv[5]) == "true", atoi(argv[6]), 50);
	//runIntegrationTest(argv[1]);
	return 0;
}//main
