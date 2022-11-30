/*
 * fixDualPolytopesInDatabaseDriver.cpp
 *
 *  Created on: Jan 13, 2011
 *      Author: bedutra
 *
 * 1st goal: Fix the dual v-reps files
 * 2nd gola: replace the dual v-reps with h-reps.
 */

#include <string>
#include <cstdio>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include "../../buildPolytopes/BuildRandomPolytope.h"
#include "../../buildPolytopes/BuildRandomPolynomials.h"
#include "../../sqlite/IntegrationDB.h"

using namespace std;


void addDualIntegrationRowsFromVolume(char *dbFile, int dim, int vertexCount)
{
	stringstream sql;
	string countStr;
	int preExistingCount;
	IntegrationDB db;
	int degrees[] = { 1, 2, 5, 10, 20, 30, 40, 50 };
	const int degreeSize = 8;


	db.open(dbFile);

	for(int i = 0; i < degreeSize; ++i)
	{

		sql.str("");
		sql << "select distinct p.rowid from polynomial as p where p.dim =" << dim << " and degree = " << degrees[i];
		vector<vector<string> > polynomialID = db.query(sql.str().c_str());

		sql.str("");
		sql << "select count(*) from "
			<<	" integrate as v"
			<<  " join polytope as dp on dp.rowid = v.polytopeID"
			<<  " join polytope as p on p.rowid = dp.dual "
			<<  " join polynomial as mon on mon.rowid = v.polynomialID"
		    <<  " where dp.dim = " << dim
		    <<  " and p.vertexCount = " << vertexCount
		    <<  " and mon.degree = " << degrees[i];
		countStr = sql.str();

		preExistingCount = db.queryAsInteger(countStr.c_str());
		if ( preExistingCount >= 50)
			continue;

		sql.str("");
		sql << "select distinct v.polytopeID from "
			<<	"volume as v"
			<<  " join polytope as dp on dp.rowid = v.polytopeID"
			<<  " join polytope as p on p.rowid = dp.dual "
		    <<  " where dp.dim = " << dim
		    <<  " and p.vertexCount = " << vertexCount;
		vector<vector<string> > polytopeID = db.query(sql.str().c_str());

		assert(polynomialID.size() >= polytopeID.size());

		int newRows = 0;
		for(int j = 0; j < polytopeID.size() ; ++j)
		{
			if ( newRows + preExistingCount >= 50)
				break;//the j loop

			//check if this is missing.
			//note, I am adding 49-50 new rows (so I don't really need to check that this polynomial has not already been used.
			sql.str("");
			sql << "select count(*) from integrate where polynomialID = " << polynomialID[j][0] << " and polytopeID = " << polytopeID[j][0];
			if ( db.queryAsInteger(sql.str().c_str()) )
				continue;

			//insert a new row.
			sql.str("");
			sql << "insert into integrate ("
				<< "polynomialID, "
				<< "polytopeID, "
				<< "timeLawrence, "
				<< "timeTriangulate, "
				<< "integral) "
				<< "values ("
				<< polynomialID[j][0] << ","
				<< polytopeID[j][0] << ","
				<< "-1,"
				<< "-1,"
				<< "'NA'"
				<< ")";

			cout << "insert sql::" << sql.str().c_str() << endl;
			++newRows;
			db.query(sql.str().c_str());


		}//for j. loop over every (polynomial/polytope) pair and add it to the integration table if missing.

		assert(db.queryAsInteger(countStr.c_str()) >= 50);

	}



	db.close();
}//addDualIntegrationRowsFromVolume

/**
 * At first, I was making dual v-reps. This is a bad idea because the files become HUGE
 * Go back, center every polymake file,
 */
void makeDualHrepFiles(char *dbFile, int dim, int vertexCount)
{
	if (dim == 3 || dim == 4)
	{
		cout	<< "makeDualHrepFiles::sorry, this test case was done before this function was written. I will not overwrite the latte files for which you already tested!!!"
				<< endl;
		return;
	}//if dim 3 or 4.

	//first get the rows that we need to update.
	IntegrationDB db;
	vector<vector<string> > rows;
	stringstream sql;
	sql		<< "select p.rowid, p.latteFilePath, p.polymakeFilePath, org.polymakeFilePath as originalPolytope"
			<< " from polytope as p, polytope as org "
			<< " where p.dual = org.rowid and p.dual is not null "
			<< " and org.dim = " << dim << " and org.vertexCount = "<< vertexCount
			<< " and p.latteFilePath not like '.vrep.dual.latte' ";
	db.open(dbFile);
	rows = db.query(sql.str().c_str());
	db.close();


	if ( rows.size() == 0)
	{
		cout << "Nothing to fix for " << dim << " " << vertexCount << endl;
		return;
	}

	//now loop over everything and make sure it exist.
	for (int i = 0; i < (int) rows.size(); ++i)
	{
		string rowid = rows[i][0];
		string dualLatteFile = rows[i][1];
		string dualPolymakeFile = rows[i][2]; //this file should not exist.
		string polymakeFile = rows[i][3];

		string ext = ".polymake";
		string dualExt = ".dual.polymake";

		//get the base file name.
		string baseFileName = polymakeFile;
		assert(string::npos == baseFileName.rfind(dualExt)); //this better by the org. polymake file.
		baseFileName.replace(baseFileName.rfind(ext), ext.length(), ""); //replace the last use of .polymake with ""

		cout << "Going to center polytope, remake latte v-file, delete dual v-file, and make dual h-file for: " << baseFileName.c_str() << endl;

		//read the org. polymake file.
		BuildPolytope bp;
		bp.setBaseFileName(baseFileName);
		bp.setBuildPolymakeFile(true);

		if (bp.getLatteVRepDualFile() != dualLatteFile
				|| bp.getPolymakeDualFile() != dualPolymakeFile)
		{
			cout << "new and old dual polymake file names differ\n"
				 << " dual latte file expected: " << bp.getLatteVRepDualFile().c_str()
				 << " dual latte file received: " << dualLatteFile.c_str()
				 << " dual polymake file expected: " << bp.getPolymakeDualFile().c_str()
				 << " dual polymkae file received: " << dualPolymakeFile.c_str();
			exit(1);
		}//if error. the database should contain vrep.dual.latte files.

		//center the polytope and build the files.
		bp.centerPolytope();
		bp.buildLatteVRepFile(); //over write latte file.
		bp.buildLatteHRepDualFile(); //make a dual h-rep
		bp.setBuildLatteVRepDualFile(true); //set to true so we can delete it
		bp.deleteLatteVRepDualFile(); //delete it.

		sql.str("");
		sql << "update polytope set  latteFilePath= '" << bp.getLatteHRepDualFile() << "'"
				<< " where rowid = " << rowid;
		cout << sql.str().c_str() << endl;
		db.open(dbFile);
		db.query(sql.str().c_str());
		db.close();
		cout << endl;
	}//for i. for every row.

}//makeDualHrepFiles

/**
 * This function should never be used again. I have here 'just in case'
 * It's job: find references to dual v-rep latte files in the db and check that they exist
 * If they do, ok. If no, make them.
 */
void makeDualFiles(char *dbFile, int dim, int vertexCount)
{
	//first get the rows that we need to update.
	IntegrationDB db;
	vector<vector<string> > rows;
	stringstream sql;
	sql
			<< "select p.rowid, p.latteFilePath, p.polymakeFilePath, org.polymakeFilePath as originalPolytope,  org.simple, p.simple, p.vertexCount "
			<< " from polytope as p, polytope as org "
			<< " where p.dual = org.rowid and p.dual is not null "
			<< " and org.dim = " << dim << " and org.vertexCount = "
			<< vertexCount;
	db.open(dbFile);
	rows = db.query(sql.str().c_str());
	db.close();

	if (rows.size() == 0)
	{
		cout << "Nothing to do for dim " << dim << " and vertexCount" << endl;
		return;
	}

	//now loop over everything and make sure it exist.
	for (int i = 0; i < (int) rows.size(); ++i)
	{
		string rowid = rows[i][0];
		string dualLatteFile = rows[i][1];
		string dualPolymakeFile = rows[i][2]; //this file should not exist.
		string polymakeFile = rows[i][3];
		string orgSimple = rows[i][4];
		string dualSimple = rows[i][5];
		string vertexCountStr = rows[i][6]; //of the dual poly.

		string ext = ".polymake";
		string dualExt = ".dual.polymake";

		ifstream dualLatte(dualLatteFile.c_str());

		cout << "\n\nProcessing " << rowid.c_str() << endl;
		if (dualLatte.is_open() && dualSimple != "-1" && vertexCountStr != "-1")
		{
			cout << "Skipping dual polytope generation." << endl;
			cout << "  Latte file exist: " << dualLatteFile.c_str() << endl;
			cout << "  Polymake file exist: " << dualPolymakeFile.c_str()
					<< endl;
			dualLatte.close();
			continue;
		}//if both files already exist. or if the database does not have a record of its vertex/simple (again, really simplicial) count.

		//get the base file name.
		string baseFileName = polymakeFile;
		assert(string::npos == baseFileName.rfind(dualExt)); //this better by the org. polymake file.
		baseFileName.replace(baseFileName.rfind(ext), ext.length(), ""); //replace the last use of .polymake with ""

		cout << "Going to build polytope and latte files for: "
				<< baseFileName.c_str() << endl;

		//read the org. polymake file.
		BuildPolytope bp;
		bp.setBaseFileName(baseFileName);
		bp.setBuildPolymakeFile(true);


		if (bp.getLatteVRepDualFile() != dualLatteFile
				|| bp.getPolymakeDualFile() != dualPolymakeFile)
		{
			cout << "new and old file names differ\n" << " dual latte v rep: "
					<< bp.getLatteVRepDualFile().c_str()
					<< " dual polymake   : "
					<< bp.getPolymakeDualFile().c_str()
					<< " db d latte v rep: " << dualLatteFile.c_str()
					<< " db dual polymake: " << dualPolymakeFile.c_str();
			exit(1);
		}//if error. somehow my naming is off :(

		//make the dual latte files
		//bp.buildPolymakeDualFile(); //don't, it will take too much space.
		bp.buildLatteVRepDualFile();

		int dualVertexCount = bp.getVertexDualCount();
		int simplicial = bp.isDualSimplicial();

		sql.str("");
		//yes, I know, all sql should go through the IntegrationDB class, but this is only a fix.
		//update the vertex count
		sql << "update polytope set vertexCount = " << dualVertexCount
				<< ", simple = " << simplicial << " where rowid = " << rowid;
		cout << sql.str().c_str() << endl;
		db.open(dbFile);
		db.query(sql.str().c_str());
		db.close();
	}//for i. for every row.


}//makeDualFiles


int main(int argc, char *argv[])
{
	if (argc <= 3)
	{
		cout
				<< "Hello,\n I was making v reps for the org. poly and for its dual."
				<< "But this is a bad idea because after dilation, the dual v-rep has HUGE numbers that slow donw the \"vertex-to-tanget-cones\""
				<< " alg. that is done when latte reads in files."
				<< "This program will go back, center the latte files, remove the dual v-reps, and make dual h-reps, and update the database."
				<< endl;

		cout << "error. usage: " << argv[0]
				<< " sqlite-db-file dim vertex-count" << endl;

		cout << "Example: " << argv[0] << " file.sqlite3 polytope 15 30\n"
				<< "        \t will build dual polymake and dual latte files for (org.) polytopes of dim 15 and vertex count 30 (the dual does not have to have 30 vertices)"
				<< endl;
		exit(1);
	}

	//There was also a point in time where I was not making dual latte files at all,
	//this function goes back and makes them. This function should never be needed again
	//because the next "fix" function supersedes it.
	// db file,   dim       vertex-count
	//	makeDualFiles(argv[1], atoi(argv[2]),atoi(argv[3]));

	//makeDualHrepFiles(argv[1], atoi(argv[2]), atoi(argv[3]));



	//I'm not sure why, but the integration table seems to be missing integration test case (the polytope dual files exist, but their dual test rows do not).
	//However, the volume table has volumes test rows for these missing dual polytopes.
	//This function asses the missing integration test cases to the integration db table.
	addDualIntegrationRowsFromVolume(argv[1], atoi(argv[2]), atoi(argv[3]));
	return 0;
}//main


