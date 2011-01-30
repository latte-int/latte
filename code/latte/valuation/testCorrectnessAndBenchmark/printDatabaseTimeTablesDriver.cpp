/*
 * printDatabaseTimeTablesDriver.cpp
 *
 *  Created on: Jan 14, 2011
 *      Author: bedutra
 *
 *
 */

#include <string>
#include <cstdio>
#include <cassert>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include "../../sqlite/IntegrationDB.h"


using namespace std;


void printHTMLtable(char * dbFile, int dim)
{

	vector<vector<IntegrationPrintData> > results, resultsDual;
	string fasterStyle, slowerStyle;
	string timeT, timeL;

	fasterStyle = " fasterTime";
	slowerStyle = " slowerTime";

	IntegrationDB db;
	db.open(dbFile);
	resultsDual = db.getStatisticsByDim(dim, true);
	results     = db.getStatisticsByDim(dim, false);
	db.close();
	assert(results.size() == resultsDual.size());

	cout << "<html><head><title>Integration results</title>" << endl;

	cout << "<style type=\"text/css\">" << endl;
	cout << ".time  {background-color: #ffffff}" << endl;
	cout << ".dualTime {background-color: #ffffff}" << endl;
	cout << ".fasterTime {color:green; font-weight:bold; text-decoration: underline;}" << endl;
	cout << ".slowerTime {color:red;}" << endl;
	cout << ".range {background-color: #EBEBEB}" << endl;
	cout << ".dualRange {background-color: #EBEBEB}" << endl;
	cout << ".countNum {background-color: #C4C4C4}" << endl;
	cout << ".dualCountNum {background-color: #C4C4C4}" << endl;

	cout << "</style></head><body>" << endl;
	cout << "Dim " << dim << ":: row number is vertex count, and column is polynomial degree<br />\n";
	cout << "<table border=\"1\">" << endl;

	//print 1st header row.
	cout << "\n<tr><td>x</td>";
	for(int i = 0; i < (int) results[0].size(); ++i)
		cout << "\n\t<td colspan=\"2\">" << results[0][i].degree << "</td>";
	cout << "\n</tr>" << endl;

	//print 2nd header row.
	cout << "\n<tr><td>x</td>";
	for(int i = 0; i < (int) results[0].size(); ++i)
			cout << "\n\t<td>Law.</td><td>Tri.</td>";
	cout << "\n</tr>" << endl;

	//print other rows.
	for(int i = 0; i < (int) results.size(); ++i)
	{
		cout << "\n<tr><!--first row-->";

		//print 1st col.
		cout << "\n\t<td rowspan=\"6\">" << results[i][0].vertexCount << "</td>";

		for(int j = 0; j < (int) results[i].size(); ++j)
		{
			if (results[i][j].avgLawrenceTime < results[i][j].avgTriangulationTime)
			{
				timeL = fasterStyle;
				timeT = slowerStyle;
			}
			else if (results[i][j].avgLawrenceTime > results[i][j].avgTriangulationTime)
			{
				timeL = slowerStyle;
				timeT = fasterStyle;
			}
			else
			{
				timeL = fasterStyle;
				timeT = fasterStyle;
			}//if tie.

			cout << "\n\t<td class=\"time" << timeL.c_str() << "\">" << results[i][j].avgLawrenceTime << "</td><td class=\"time" << timeT.c_str() << "\">" << results[i][j].avgTriangulationTime << "</td>";
		}
		cout << "\n</tr>";

		cout << "\n<tr><!--min max row-->";
		for(int j = 0; j < (int) results[i].size(); ++j)
			cout << "\n\t<td class=\"range\">[" << results[i][j].minLawrenceTime  << ", " << results[i][j].maxLawrenceTime << "]</td><td class=\"range\">[" << results[i][j].minTriangulationTime << ", " << results[i][j].maxTriangulationTime << "]</td>";
		cout << "\n</tr>";


		cout << "\n<tr><!--count row-->";
		for(int j = 0; j < (int) results[i].size(); ++j)
			cout << "\n\t<td class=\"countNum\">[" << results[i][j].totalFinishedLawrenceTestCases  << "/" << results[i][j].totalTestCases << "</td><td class=\"countNum\">" << results[i][j].totalFinishedTriangulationTestCases << "/" << results[i][j].totalTestCases << "</td>";
		cout << "\n</tr>";

		//dual rows.
		cout << "\n<tr><!--dual avg row-->";
		for(int j = 0; j < (int) resultsDual[i].size(); ++j)
		{
			if (resultsDual[i][j].avgLawrenceTime < resultsDual[i][j].avgTriangulationTime)
			{
				timeL = fasterStyle;
				timeT = slowerStyle;
			}
			else if (resultsDual[i][j].avgLawrenceTime > resultsDual[i][j].avgTriangulationTime)
			{
				timeL = slowerStyle;
				timeT = fasterStyle;
			}
			else
			{
				timeL = fasterStyle;
				timeT = fasterStyle;
			}//if tie.

			cout << "\n\t<td class=\"dualTime" << timeL.c_str() << "\">" << resultsDual[i][j].avgLawrenceTime << "</td><td class=\"dualTime" << timeT.c_str() << "\">" << resultsDual[i][j].avgTriangulationTime << "</td>";
		}
		cout << "\n</tr>";

		cout << "\n<tr><!--dual min max row-->";
		for(int j = 0; j < (int) resultsDual[i].size(); ++j)
			cout << "\n\t<td class=\"dualRange\">[" << resultsDual[i][j].minLawrenceTime  << ", " << resultsDual[i][j].maxLawrenceTime << "]</td><td class=\"dualRange\">[" << resultsDual[i][j].minTriangulationTime << ", " << resultsDual[i][j].maxTriangulationTime << "]</td>";
		cout << "\n</tr>";


		cout << "\n<tr><!--dual count row-->";
		for(int j = 0; j < (int) resultsDual[i].size(); ++j)
			cout << "\n\t<td class=\"dualCountNum\">" << resultsDual[i][j].totalFinishedLawrenceTestCases  << "/" << resultsDual[i][j].totalTestCases << "</td><td class=\"dualCountNum\">" << resultsDual[i][j].totalFinishedTriangulationTestCases << "/" << resultsDual[i][j].totalTestCases << "</td>";
		cout << "\n</tr>";



	}//for each row.


	cout << "</table>" << endl;
	cout << "</body></head>" << endl;
}//printHTMLtable()


void printHTMLtablePerPolymakeFile(const char *dbFile)
{

	vector<vector<IntegrationPrintData> > results, resultsDual;
	vector<vector<string> > allPolymakeFiles;
	string fasterStyle, slowerStyle;
	string timeT, timeL;

	fasterStyle = " fasterTime";
	slowerStyle = " slowerTime";

	IntegrationDB db;
	db.open(dbFile);
	allPolymakeFiles = db.getAllPolymakeFiles();
	resultsDual = db.getStatisticsByFile(allPolymakeFiles, true);
	results     = db.getStatisticsByFile(allPolymakeFiles, false);
	db.close();
	assert(results.size() == resultsDual.size());

	cout << "<html><head><title>Integration results</title>" << endl;

	cout << "<style type=\"text/css\">" << endl;
	cout << ".time  {background-color: #ffffff}" << endl;
	cout << ".dualTime {background-color: #ffffff}" << endl;
	cout << ".fasterTime {color:green; font-weight:bold; text-decoration: underline;}" << endl;
	cout << ".slowerTime {color:red;}" << endl;
	cout << ".range {background-color: #EBEBEB}" << endl;
	cout << ".dualRange {background-color: #EBEBEB}" << endl;
	cout << ".countNum {background-color: #C4C4C4}" << endl;
	cout << ".dualCountNum {background-color: #C4C4C4}" << endl;

	cout << "</style></head><body>" << endl;
	//cout << " " << ) << ":: row number is vertex count, and column is polynomial degree<br />\n";
	cout << "<table border=\"1\">" << endl;

	//print 1st header row.
	cout << "\n<tr><td>x</td>";
	for(int i = 0; i < (int) results[0].size(); ++i)
		cout << "\n\t<td colspan=\"2\">" << results[0][i].degree << "</td>";
	cout << "\n</tr>" << endl;

	//print 2nd header row.
	cout << "\n<tr><td>x</td>";
	for(int i = 0; i < (int) results[0].size(); ++i)
			cout << "\n\t<td>Law.</td><td>Tri.</td>";
	cout << "\n</tr>" << endl;

	//print other rows.
	for(int i = 0; i < (int) results.size(); ++i)
	{
		if ( results[i][0].avgLawrenceTime == 0 && results[i][0].avgTriangulationTime == 0)
			continue; //don' print empty rows.
		cout << "\n<tr><!--first row-->";

		//print 1st col.
		cout << "\n\t<td rowspan=\"6\">" << allPolymakeFiles[i][0].c_str() << "<br>dim: " << results[i][0].dim << " vertex: " << results[i][0].vertexCount << "</td>";

		for(int j = 0; j < (int) results[i].size(); ++j)
		{
			if (results[i][j].avgLawrenceTime < results[i][j].avgTriangulationTime)
			{
				timeL = fasterStyle;
				timeT = slowerStyle;
			}
			else if (results[i][j].avgLawrenceTime > results[i][j].avgTriangulationTime)
			{
				timeL = slowerStyle;
				timeT = fasterStyle;
			}
			else
			{
				timeL = fasterStyle;
				timeT = fasterStyle;
			}//if tie.

			cout << "\n\t<td class=\"time" << timeL.c_str() << "\">" << results[i][j].avgLawrenceTime << "</td><td class=\"time" << timeT.c_str() << "\">" << results[i][j].avgTriangulationTime << "</td>";
		}
		cout << "\n</tr>";

		cout << "\n<tr><!--min max row-->";
		for(int j = 0; j < (int) results[i].size(); ++j)
			cout << "\n\t<td class=\"range\">[" << results[i][j].minLawrenceTime  << ", " << results[i][j].maxLawrenceTime << "]</td><td class=\"range\">[" << results[i][j].minTriangulationTime << ", " << results[i][j].maxTriangulationTime << "]</td>";
		cout << "\n</tr>";


		cout << "\n<tr><!--count row-->";
		for(int j = 0; j < (int) results[i].size(); ++j)
			cout << "\n\t<td class=\"countNum\">[" << results[i][j].totalFinishedLawrenceTestCases  << "/" << results[i][j].totalTestCases << "</td><td class=\"countNum\">" << results[i][j].totalFinishedTriangulationTestCases << "/" << results[i][j].totalTestCases << "</td>";
		cout << "\n</tr>";

		//dual rows.
		cout << "\n<tr><!--dual avg row-->";
		for(int j = 0; j < (int) resultsDual[i].size(); ++j)
		{
			if (resultsDual[i][j].avgLawrenceTime < resultsDual[i][j].avgTriangulationTime)
			{
				timeL = fasterStyle;
				timeT = slowerStyle;
			}
			else if (resultsDual[i][j].avgLawrenceTime > resultsDual[i][j].avgTriangulationTime)
			{
				timeL = slowerStyle;
				timeT = fasterStyle;
			}
			else
			{
				timeL = fasterStyle;
				timeT = fasterStyle;
			}//if tie.

			cout << "\n\t<td class=\"dualTime" << timeL.c_str() << "\">" << resultsDual[i][j].avgLawrenceTime << "</td><td class=\"dualTime" << timeT.c_str() << "\">" << resultsDual[i][j].avgTriangulationTime << "</td>";
		}
		cout << "\n</tr>";

		cout << "\n<tr><!--dual min max row-->";
		for(int j = 0; j < (int) resultsDual[i].size(); ++j)
			cout << "\n\t<td class=\"dualRange\">[" << resultsDual[i][j].minLawrenceTime  << ", " << resultsDual[i][j].maxLawrenceTime << "]</td><td class=\"dualRange\">[" << resultsDual[i][j].minTriangulationTime << ", " << resultsDual[i][j].maxTriangulationTime << "]</td>";
		cout << "\n</tr>";


		cout << "\n<tr><!--dual count row-->";
		for(int j = 0; j < (int) resultsDual[i].size(); ++j)
			cout << "\n\t<td class=\"dualCountNum\">" << resultsDual[i][j].totalFinishedLawrenceTestCases  << "/" << resultsDual[i][j].totalTestCases << "</td><td class=\"dualCountNum\">" << resultsDual[i][j].totalFinishedTriangulationTestCases << "/" << resultsDual[i][j].totalTestCases << "</td>";
		cout << "\n</tr>";



	}//for each row.


	cout << "</table>" << endl;
	cout << "</body></head>" << endl;
}//printHTMLtablePerPolymakeFile



int main(int argc, char *argv[])
{
	if (argc < 2 )
	{
		cout << "Hello,\n "
			 << "This program will print time tables for easy viewing"
			 << endl;

		cout << "error. usage: " << argv[0] << " sqlite-db-file dim (optional)" << endl;
		exit(1);
	}
	if ( argc == 3)
					// db file,   dim
		printHTMLtable(argv[1], atoi(argv[2]));
	else if ( argc == 2)
		printHTMLtablePerPolymakeFile(argv[1]);

	return 0;
}//main



