/*
 * buildDatabaseRandomPolytopesDriver.cpp
 *
 *  Created on: Nov 24, 2010
 *      Author: Brandon Dutra and Gregory Pinto
 *
 *      Computes H-rep latte files.
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


struct BuildClass
{
	int dim;
	int degree; //of polynomial
	int vertex; //number of vertices in polytope.
	int number; //how many of this
};


		// db file,  count         polymake file, degree
void buildSpecficPolytopeAndTestCase(char * dbFile, int count, char * polymakeFile, int degree)
{
	string baseFileName(polymakeFile);
	string fileEnding(".polymake");
	baseFileName.replace(baseFileName.find(fileEnding),fileEnding.length(),"");
	//now baseFileName does not contain the .polymake file extention.

	BuildPolytope newPolytope;
	newPolytope.setBaseFileName(baseFileName);
	newPolytope.setBuildPolymakeFile(true);
	newPolytope.findVertices();
	newPolytope.findFacets();

	//first insert the polytopes into the polytope table.
	IntegrationDB db;
	db.open(dbFile);
	bool doesPolytopeExist =db.doesPolytopeExist(newPolytope.getPolymakeFile().c_str());
	db.close(); //close because next part can take a long time.
	if ( doesPolytopeExist )
	{
		cout << newPolytope.getPolymakeFile().c_str() << " already exist in database, skipping" << endl;
	}
	else
	{
		int dualVertexCount, buildVertexCount, buildDim;
		int simple, dualSimple;

		//make the polymake file and latte file.
		cout << "**" << newPolytope.getPolymakeFile().c_str() << " degree " << degree << endl;
		newPolytope.buildLatteVRepFile();
		newPolytope.buildPolymakeDualFile();
		newPolytope.buildLatteVRepDualFile();


		//get info about this polytope.
		buildVertexCount = newPolytope.getVertexCount();
		buildDim         = newPolytope.getDim();
		dualVertexCount  = newPolytope.getVertexDualCount();
		simple           = -1;//save time, don't compute this. newPolytope.isSimplicial(); //again, historical naming error.
		dualSimple       = -1;//newPolytope.isDualSimplicial();

		//print to screen --debugging.
		cout << newPolytope.getLatteVRepFile().c_str() << " dim: " << buildDim << "\tvertex " << buildVertexCount << "\tsimple " << simple << endl;
		cout << "  polymake: " << newPolytope.getPolymakeFile().c_str() << endl;
		cout << "  " << newPolytope.getLatteVRepDualFile().c_str() << " dim: " << buildDim << "\tvertex " << dualVertexCount << "\tsimple " << dualSimple << endl;
		cout << "  polymake: " << newPolytope.getPolymakeDualFile().c_str() << endl;

		//save the polytopes in the database.
		db.open(dbFile);
		int dualID;
		dualID = db.insertPolytope(buildDim, buildVertexCount, simple    , -1    , newPolytope.getLatteVRepFile().c_str()    , newPolytope.getPolymakeFile().c_str());
		         db.insertPolytope(buildDim, dualVertexCount, dualSimple, dualID, newPolytope.getLatteVRepDualFile().c_str(), newPolytope.getPolymakeDualFile().c_str());
		db.close();
		newPolytope.deletePolymakeDualFile();
	}//insert this polytope into the database.

	//now that the polytope is in the database, lets add testcases for the given polynomial degree.
	//(we assume the db already has enough polynomials...that is, we do not make any now).


	db.open(dbFile);
	db.insertSpecficPolytopeIntegrationTest(newPolytope.getPolymakeFile(), degree, count);
	db.insertSpecficPolytopeIntegrationTest(newPolytope.getPolymakeDualFile(), degree, count);
	db.close();

}//buildSpecficPolytopeAndTestCase


void buildIntegrationTest(char *dbFile, int count, int dim, int vertexCount, int degree)
{
/*
	cout << "Enter list of test classes to build (dim vertex-count degree count\n Example:"<<endl;
	cout << "> 2 3 4 5 " << endl;
	cout << "> -1"<< endl;
	cout << "will build 5 test cases of dim-2 polytopes which have 3 vertices on polynomials of degree 4" << endl;

	//I don't really need to do use the BuildClass or the get-then-build style, but I'll keep things similar to how stuff is done in the other functions.
	int dim;
	vector<BuildClass> toBuildList;
	cout << "> ";
	cin >> dim;
	while(dim > 0)
	{
		BuildClass b;
		b.dim = dim;
		cin >> b.vertex >> b.degree >> b.number;
		toBuildList.push_back(b);

		cout << "> ";
		cin >> dim;
	}//while there is more testing classes to build.
*/
	vector<BuildClass> toBuildList;
	BuildClass b;
	b.dim = dim;
	b.vertex = vertexCount;
	b.degree = degree;
	b.number = count;
	toBuildList.push_back(b);

	for(int i = 0; i < (int) toBuildList.size(); ++i)
	{
		cout << "Going to build " << count << " many tests using: "
			 << "\n dim " << dim
			 << "\n vertexCount " << vertexCount
			 << "\n degree " << degree
			 << endl;
		IntegrationDB db;
		db.open(dbFile);
		db.insertIntegrationTest(toBuildList[i].dim, toBuildList[i].degree, toBuildList[i].vertex, toBuildList[i].number);
	}//for i
}//buildIntegrationTest

//build polytopes and insert them with unused polynomials for the class.
void buildPolytopes(char *dbFile, int count, int dim, int vertexCount)
{
/*
	cout << "Enter list of polytope classes to build (vertex-count dim count)\n Example:" << endl;
	cout << "> 5 2 10" << endl;
	cout << "> -1" << endl;
	cout << "Will build 10 polytopes in folder dim2 of dim-2 with 5 vertices" << endl;

	int vertex;
	vector<BuildClass> toBuildList;
	cout << "> ";
	cin >> vertex;
	while ( vertex >0)
	{
		BuildClass b;
		b.vertex = vertex;
		cin >> b.dim >> b.number;
		toBuildList.push_back(b);

		cout << "> ";
		cin >> vertex;
	}//while there are more build classes
*/
	vector<BuildClass> toBuildList;
	BuildClass b;
	b.dim = dim;
	b.vertex = vertexCount;
	b.number = count;
	toBuildList.push_back(b);

	//now build everything!
	for(int i = 0; i < (int) toBuildList.size(); ++i)
	{
		for(; toBuildList[i].number > 0; --toBuildList[i].number)
		{
			//make the name
			stringstream baseFileName;
			baseFileName << "dim" << toBuildList[i].dim <<  "/ext" << toBuildList[i].vertex << "/polytope" << toBuildList[i].dim << "_vertex" << toBuildList[i].vertex << "_num" << toBuildList[i].number;

			//try to build the polytope
			bool correctFlag = false;
			int additionalPoints = 0;
			while (! correctFlag)
			{
				BuildRandomPolytope newPolytope;
				int buildVertexCount, buildDim;

				newPolytope.setBaseFileName(baseFileName.str());

				//check that this polytope does not already esist
				IntegrationDB db;
				db.open(dbFile);
				if ( db.doesPolytopeExist(newPolytope.getPolymakeFile().c_str()))
				{
					cout << newPolytope.getPolymakeFile().c_str() << " already exist in database, skipping" << endl;
					correctFlag = true;
					continue;//for the while loop.
				}
				db.close(); //close it because the next part can take a really long time!

				//make the polymake file and latte file.
				cout << "**dim = " << toBuildList[i].dim << " vertex=" << max(toBuildList[i].vertex, toBuildList[i].vertex + additionalPoints) <<endl;
				newPolytope.makePoints(toBuildList[i].dim, max(toBuildList[i].vertex, toBuildList[i].vertex + additionalPoints), 50, 0.5);
				newPolytope.buildLatteVRepFile();//also makes a polymake file

				//get info about this polytope.
				buildVertexCount = newPolytope.getVertexCount();
				buildDim         = newPolytope.getDim();

				if ( buildVertexCount != toBuildList[i].vertex || buildDim != toBuildList[i].dim)
				{
					cout << "Deleting " << baseFileName.str().c_str() << " because vertex count=" << buildVertexCount << "and dim=" << buildDim << " additionalPoints=" << additionalPoints << endl;
					newPolytope.deletePolymakeFile();
					newPolytope.deleteLatteVRepFile();

					//if the vertex count is wrong, adjust the number of points used.
					additionalPoints += toBuildList[i].vertex - buildVertexCount; //will add/sub points if the current polytope's vertex count is too low/high.
					continue; //for the while loop.
				}//if polytope is not built to the requirements.

				//now make the latte file and the dual polymake and dual latte file.
				cout << "dual not being made";
				newPolytope.buildPolymakeDualFile();
				newPolytope.buildLatteVRepDualFile();

				//collect more statistics
				int dualVertexCount;
				int simple, dualSimple;
				dualVertexCount = newPolytope.getVertexDualCount();
				simple          = newPolytope.isSimplicial();
				dualSimple      = newPolytope.isDualSimplicial();

				//print to screen --debugging.
				cout << newPolytope.getLatteVRepFile().c_str() << " dim: " << buildDim << "\tvertex " << buildVertexCount << "\tsimple " << simple << endl;
				cout << "  polymake: " << newPolytope.getPolymakeFile().c_str() << endl;
				cout << "  " << newPolytope.getLatteVRepDualFile().c_str() << " dim: " << buildDim << "\tvertex " << dualVertexCount << "\tsimple " << dualSimple << endl;
				cout << "  polymake: " << newPolytope.getPolymakeDualFile().c_str() << endl;
				correctFlag = true;

				//save the polytopes in the database.
				db.open(dbFile);
				int dualID;
				dualID = db.insertPolytope(toBuildList[i].dim, toBuildList[i].vertex, simple    , -1    , newPolytope.getLatteVRepFile().c_str()    , newPolytope.getPolymakeFile().c_str());
				         db.insertPolytope(toBuildList[i].dim, dualVertexCount      , dualSimple, dualID, newPolytope.getLatteVRepDualFile().c_str(), newPolytope.getPolymakeDualFile().c_str());

				newPolytope.deletePolymakeDualFile();
				//yes, we are done! Let the deconstructor clean things up.
			}//while the current polytope is not of the correct dim or vertex-count, try again.
		}//for every polytope
	}//for i. for every build class.
}//buildPolytopes

//saves the polynomials in ./polynomials/dimN/
void buildPolynomials(char *dbFile, int count, int dim, int degree)
{
/*	cout << "Enter list of polynomial classes to build (degree dim count)\n Example:"<<endl;
	cout << "> 4 10 50" << endl;
	cout << "> -1" <<endl;
	cout << "Will build 50 degree 4, dim-10 polynomials in folder polynomials/dim10" << endl;

	int degree;
	vector<BuildClass> toBuildList;
	cout << "> ";
	cin >> degree;
	while(degree > 0)
	{
		BuildClass b;
		b.degree = degree;
		cin >> b.dim >> b.number;
		toBuildList.push_back(b);

		cout << "> ";
		cin >> degree;
	}//while another case to bild.
*/
	vector<BuildClass> toBuildList;
	BuildClass b;
	b.dim = dim;
	b.degree = degree;
	b.number = count;
	toBuildList.push_back(b);

	//connect to the db.
	IntegrationDB db;
	db.open(dbFile);

	//for each class of polynomials, build them!
	for(int i = 0; i < (int)toBuildList.size(); ++i)
	{
		for(; toBuildList[i].number > 0; --toBuildList[i].number)
		{
			stringstream fileName;
			fileName << "polynomials/dim" << toBuildList[i].dim  << "/deg" << toBuildList[i].degree << "/poly" << toBuildList[i].dim << "_deg" << toBuildList[i].degree << "_num" << toBuildList[i].number << ".polynomial";

			//insert the name in the db.
			int newRowid;
			try {newRowid = db.insertPolynomial(toBuildList[i].dim, toBuildList[i].degree, fileName.str().c_str());}
			catch (exception &e)
			{
				cout << e.what() << endl;//maybe the file name already exist?
				continue;
			}//catch

			//make the file.
			ofstream file(fileName.str().c_str());
			if ( file.is_open())
			{
				file << makeRandomPolynomial(toBuildList[i].dim, toBuildList[i].degree, 1).c_str() << endl;
				file.close();
			}//we can make the file, do it.
			else
			{
				db.deletePolynomial(newRowid);
				cout << fileName.str().c_str() << " cannot be written to. removed from db." << endl;
				continue;
			}//else, cannot write file, so delete it from the db.

			cout << fileName.str().c_str() << endl;
		}//for each polynomial number

	}//for i. for each build class
}//buildPolynomials

int main(int argc, char *argv[])
{
	if (argc <= 2 )
	{
		cout << "Hello,\n This program can"
			 << "\n\t1) Make polynomials in a set folder (see below), and update a sqlite database"
			 << "\n\t2) Make polytopes (with polymake) in a set folder (see below),, and update a sqlite database"
			 << "\n\t3) Insert test cases into the sqlite database" << endl;

		cout << "error. usage: " << argv[0] << " sqlite-db-file [polytope | polynomial | test | specficPolytoepAndTest ] [more parameters]" << endl;
		cout << "  polytope parameters:   count dim vertex-count " << endl;
		cout << "  polynomial parameters: count dim degree" << endl;
		cout << "  test parameters:       count dim vertex-count degree\n" << endl;
		cout << "  specficPolytoepAndTest parameters : count polytope-file-path degree" << endl;

		cout << "Example: " << argv[0] << " file.sqlite3 polytope 50 4 8\n"
			 << "\t will make 50 dim-4 polytopes with 8 vertices in folder dim4/ext4/" << endl;
		cout << "Example: " << argv[0] << " file.sqlite3 polyomial 50 4 8\n"
			 << "\t will make 50 dim-4 degree-8 polynomials in the folder polynomials/dim4/deg8" << endl;
		cout << "Example: " << argv[0] << " file.sqlite3 test 50 4 8 10"
			 << "\t will make 50 test covering dim-4 polytopes with 8 vertices on polynomials of dim-4 and degree 10" << endl;
		cout << "Example: " << argv[0] << " file.sqlite3 specficPolytoepAndTest 50 dir/Birkhoff.polymake 9\n"
			 << "\t will make 50 polynomial-degree 9 tests for this polytope, and will make the latte and dual latte vrep files. The polymake file should already exist, and the db should have polynomials already.\n"
			 << "\t The polymake file must end in .polymake and is assumed to NOT to be dual." << endl;


		exit(1);
	}

	if (strcmp(argv[2], "polynomial") == 0)
					// db file,   count,       dim              degree
		buildPolynomials(argv[1], atoi(argv[3]),atoi(argv[4]), atoi(argv[5]));
	else if (strcmp(argv[2], "polytope") == 0)
		            // db file,   count,       dim              vertex-count
		buildPolytopes(argv[1], atoi(argv[3]),atoi(argv[4]), atoi(argv[5]));
	else if (strcmp(argv[2], "test") == 0)
					       // db file,   count,       dim              vertex-count, degree
		buildIntegrationTest(argv[1], atoi(argv[3]),atoi(argv[4]), atoi(argv[5]), atoi(argv[6]));
	else if ( strcmp(argv[2], "specficPolytoepAndTest") == 0 )
									// db file,  count         polymake file, degree
		buildSpecficPolytopeAndTestCase(argv[1], atoi(argv[3]), argv[4], atoi(argv[5]));
	else
	{
		cout << "unknown option" << endl;
	}

	return 0;
}//main


