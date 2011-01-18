/*
 * fixDualPolytopesInDatabaseDriver.cpp
 *
 *  Created on: Jan 13, 2011
 *      Author: bedutra
 *
 *  When making the database of polytopes, I was inserting file paths to dual polytoes
 *  but I was not making the files on disk because I did not know how.
 *
 *  Now I do. This program will look at all the dual polytopes in the database
 *  and make their latte file if it does not already exist.
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

// db file,   dim       vertex-count
void makeDualFiles(char *dbFile, int dim, int vertexCount)
{
	//first get the rows that we need to update.
	IntegrationDB db;
	vector<vector<string> > rows;
	stringstream sql;
	sql << "select p.rowid, p.latteFilePath, p.polymakeFilePath, org.polymakeFilePath as origionalPolytope,  org.simple, p.simple, p.vertexCount "
		<< " from polytope as p, polytope as org "
		<< " where p.dual = org.rowid and p.dual is not null "
		<< " and org.dim = " << dim
		<< " and org.vertexCount = " << vertexCount;
	db.open(dbFile);
	rows = db.query(sql.str().c_str());
	db.close();

	if ( rows.size() == 0)
	{
		cout << "Nothing to do for dim " << dim << " and vertexCount" << endl;
		return;
	}

	//now loop over everything and make sure it exist.
	for(int i = 0; i < (int) rows.size(); ++i)
	{
		string rowid            = rows[i][0];
		string dualLatteFile    = rows[i][1];
		string dualPolymakeFile = rows[i][2]; //this file should not exist.
		string polymakeFile     = rows[i][3];
		string orgSimple        = rows[i][4];
		string dualSimple       = rows[i][5];
		string vertexCountStr   = rows[i][6]; //of the dual poly.

		string ext = ".polymake";

		ifstream dualLatte(dualLatteFile.c_str());

		cout << "\n\nProcessing " << rowid.c_str() << endl;
		if ( dualLatte.is_open()
			&& dualSimple != "-1" && vertexCountStr != "-1")
		{
			cout << "Skipping dual polytope generation." << endl;
			cout << "  Latte file exist: " << dualLatteFile.c_str() << endl;
			cout << "  Polymake file exist: " << dualPolymakeFile.c_str() << endl;
			dualLatte.close();
			continue;
		}//if both files already exist. or if the database does not have a record of its vertex/simple (again, really simplicial) count.

		string baseFileName = polymakeFile;
		baseFileName.replace( baseFileName.rfind(ext), ext.length(),""); //replace the last use of .polymake with ""

		cout << "Going to build polytope and latte files for: " << baseFileName.c_str() << endl;

		BuildPolytope bp;
		bp.setBaseFileName(baseFileName);
		bp.setBuildPolymakeFile(true);

		if ( bp.getLatteVRepDualFile() != dualLatteFile
				|| bp.getPolymakeDualFile() != dualPolymakeFile)
		{
			cout << "new and old file names differ\n"
				 << " dual latte v rep: " << bp.getLatteVRepDualFile().c_str()
				 << " dual polymake   : " << bp.getPolymakeDualFile().c_str()
				 << " db d latte v rep: " << dualLatteFile.c_str()
				 << " db dual polymake: " << dualPolymakeFile.c_str();
			exit(1);
		}

		bp.buildPolymakeDualFile();
		bp.buildLatteVRepDualFile();

		int dualVertexCount = bp.getVertexDualCount();
		int simplicial = bp.isDualSimplicial();

		sql.str("");
		//yes, I know, all sql should go through the IntegrationDB class, but this is only a fix.
		sql << "update polytope set vertexCount = " << dualVertexCount
			<< ", simple = " << simplicial
			<< " where rowid = " << rowid;
		cout << sql.str().c_str() << endl;
		db.open(dbFile);
		db.query(sql.str().c_str());
		db.close();
	}//for i. for every row.


}//makeDualFiles



int main(int argc, char *argv[])
{
	if (argc <= 3 )
	{
		cout << "Hello,\n I was incorrectly making dual polytope files."
			 << " The database contains file paths to dual polymake and dual "
			 << " vrep files. This program will go back and physically build the "
			 << " latte files that the database says exist. Dual polymake files should not exist (waste of space)\n"
			 << endl;

		cout << "error. usage: " << argv[0] << " sqlite-db-file dim vertex-count" << endl;

		cout << "Example: " << argv[0] << " file.sqlite3 polytope 15 30\n"
			 << "        \t will build dual polymake and dual latte files for (org.) polytopes of dim 15 and vertex count 30 (the dual does not have to have 30 vertices)" << endl;
		exit(1);
	}
					// db file,   dim       vertex-count
	makeDualFiles(argv[1], atoi(argv[2]),atoi(argv[3]));

	return 0;
}//main


