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
#include <iomanip>
#include "../../sqlite/IntegrationDB.h"
#include "../../sqlite/VolumeDB.h"


using namespace std;
void printLatexTableIntegration(char * dbFile, int dim)
{

	vector<vector<ValuationDBStatistics> > results, resultsDual;
	string fasterStyle, slowerStyle, sameStyle;
	string timeT, timeL;

	fasterStyle = " \\timeFaster";
	slowerStyle = " \\timeSlower";
	sameStyle   = " \\timeBothZero";

	IntegrationDB db;
	db.open(dbFile);
	resultsDual = db.getStatisticsByDim(dim, true);
	results     = db.getStatisticsByDim(dim, false);
	db.close();
	assert(results.size() == resultsDual.size());

	cout << fixed;

	for(int i = 0; i < (int) results.size(); ++i)
	{
		for(int j = 0; j < (int) results[i].size(); ++j)
		{
			if ( results[i][j].totalFinishedLawrenceTestCases < 50)
				results[i][j].avgLawrenceTime = 0;
			if ( results[i][j].totalFinishedTriangulationTestCases < 50 )
				results[i][j].avgTriangulationTime = 0;
		}
		for(int j = 0; j < (int) resultsDual[i].size(); ++j)
		{
			if ( resultsDual[i][j].totalFinishedLawrenceTestCases < 50)
				resultsDual[i][j].avgLawrenceTime = 0;
			if ( resultsDual[i][j].totalFinishedTriangulationTestCases < 50 )
				resultsDual[i][j].avgTriangulationTime = 0;
		}
	}//for i. clear any incomplete test cases.


	cout << "\\begin{table}\n";
	cout << "\\caption{Triangulation vs Lawrence Integration on Random Polytopes in Dimension " << dim << "}\n";
	cout << "\\label{tabel:lawrence-random-integration-dim" << dim <<  "}\n";
	cout << "\\begin{tabular}{l";
	for(int i = 0; i < results.size(); ++i)
		cout << 'r';
	cout << "}\n";
	cout << "\\toprule \n";
	//print 1st header row.
	cout << "& \\multicolumn{" << results.size() << "}{c}{Vertex Count}\\\\ \n";
	cout << " \\cmidrule(c){2-" << results.size() + 1 << "} \n";

	for (int i = 0; i < (int) results.size(); ++i)
		cout << " & \\multicolumn{2}{c}{" << results[i][0].vertexCount << "} ";
	cout << " \\\\ \n";
	for (int i = 0; i < (int) results.size(); ++i)
		cout << "\\cmidrule(r){" << 2*i+1+1 << "-" << 2*i+2+1<< "}";
	cout << "\n";

	//print 2nd header row.
	cout << "Degree";
	for(int i = 0; i < (int) results[0].size(); ++i)
			cout << " & Law. & Tri. ";
	cout << "\\\\ \n";

	//print other rows.
	for(int i = 0; i < (int) results[0].size(); ++i)
	{
		//print 1st col.
		cout << results[0][i].degree;

		for(int j = 0; j < (int) results.size(); ++j)
		{
			if (results[j][i].avgLawrenceTime < results[j][i].avgTriangulationTime)
			{
				timeL = fasterStyle;
				timeT = slowerStyle;
			}
			else if (results[j][i].avgLawrenceTime > results[j][i].avgTriangulationTime)
			{
				timeL = slowerStyle;
				timeT = fasterStyle;
			}
			else
			{
				timeL = sameStyle;
				timeT = sameStyle;
			}//if tie.

			cout << " & " << timeL.c_str() << "{ " << setw(6) << setprecision(2) << results[j][i].avgLawrenceTime << "}"
				 << " & " << timeT.c_str() << "{ " << setw(6) << setprecision(2) << results[j][i].avgTriangulationTime << "}";
		}
		cout << "\\\\ \n";

		//dual rows.
		for(int j = 0; j < (int) resultsDual.size(); ++j)
		{
			if (resultsDual[j][i].avgLawrenceTime < resultsDual[j][i].avgTriangulationTime)
			{
				timeL = fasterStyle;
				timeT = slowerStyle;
			}
			else if (resultsDual[j][i].avgLawrenceTime > resultsDual[j][i].avgTriangulationTime)
			{
				timeL = slowerStyle;
				timeT = fasterStyle;
			}
			else
			{
				timeL = sameStyle;
				timeT = sameStyle;
			}//if tie.

			cout << " & " << timeL.c_str() << "{ " << setw(6) << setprecision(2) << resultsDual[j][i].avgLawrenceTime << "}"
				 << " & " << timeT.c_str() << "{ " << setw(6) << setprecision(2) << resultsDual[j][i].avgTriangulationTime << "}";
		}
		cout << "\\\\ \n";
		cout << "\\hline \n";
	}//for each row.

	cout << "\\bottomrule \n";
	cout << "\\end{tabular} \n";
	cout << "\\end{table} \n";
}//printLatexIntegrationTable()


void printLatexTable_core(const vector<vector<ValuationDBStatistics> > & results,
		const vector<vector<ValuationDBStatistics> >  &resultsDual,
		const string &fasterStyle,
		const string &slowerStyle,
		const string &sameStyle,
		const string &zeroStyle)
{
	string timeT, timeL;

	for(int i = 0; i < (int) results.size(); ++i)
	{
		//print 1st col.
		cout << results[i][0].vertexCount;

		for(int j = 0; j < (int) results[i].size(); ++j)
		{
			if ( results[i][j].avgLawrenceTime == 0 )
			{
				if ( results[i][j].avgTriangulationTime != 0)
				{
					timeT = fasterStyle;
					timeL = zeroStyle;
				}
				else
				{
					timeT = zeroStyle;
					timeL = zeroStyle;
				}
			}
			else if ( results[i][j].avgTriangulationTime == 0)
			{
				timeT = zeroStyle;
				timeL = fasterStyle; //not zero.
			}
			else if (results[i][j].avgLawrenceTime < results[i][j].avgTriangulationTime)
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
				timeL = sameStyle;
				timeT = sameStyle;
			}//tie.

			cout << " & " << timeL.c_str() << "{ " << setw(6) << setprecision(2) << results[i][j].avgLawrenceTime << "}"
				 << " & " << timeT.c_str() << "{ " << setw(6) << setprecision(2) << results[i][j].avgTriangulationTime << "}";
		}
		cout << "\\\\ \n";

		//dual rows.
		for(int j = 0; j < (int) resultsDual[i].size(); ++j)
		{
			if ( resultsDual[i][j].avgLawrenceTime == 0 )
			{
				if ( resultsDual[i][j].avgTriangulationTime != 0)
				{
					timeT = fasterStyle;
					timeL = zeroStyle;
				}
				else
				{
					timeT = zeroStyle;
					timeL = zeroStyle;
				}
			}
			else if ( resultsDual[i][j].avgTriangulationTime == 0)
			{
				timeT = zeroStyle;
				timeL = fasterStyle; //not zero.
			}
			else if (resultsDual[i][j].avgLawrenceTime < resultsDual[i][j].avgTriangulationTime)
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
				timeL = sameStyle;
				timeT = sameStyle;
			}//tie.

			cout << " & " << timeL.c_str() << "{ " << setw(6) << setprecision(2) << resultsDual[i][j].avgLawrenceTime << "}"
				 << " & " << timeT.c_str() << "{ " << setw(6) << setprecision(2) << resultsDual[i][j].avgTriangulationTime << "}";
		}
		cout << "\\\\ \n";
		cout << "\\hline \n";
	}//for each row.

}//printLatexTable_core

void printLatexTableIntegration_sideways(char * dbFile, int dim)
{

	vector<vector<ValuationDBStatistics> > results, resultsDual;
	string fasterStyle, slowerStyle, sameStyle, zeroStyle;


	fasterStyle = " \\timeFaster";
	slowerStyle = " \\timeSlower";
	sameStyle   = " \\timeTie";
	zeroStyle   = " \\timeNotComputed";

	IntegrationDB db;
	db.open(dbFile);
	resultsDual = db.getStatisticsByDim(dim, true);
	results     = db.getStatisticsByDim(dim, false);
	db.close();
	assert(results.size() == resultsDual.size());

	cout << fixed;

	for(int i = 0; i < (int) results.size(); ++i)
	{
		for(int j = 0; j < (int) results[i].size(); ++j)
		{
			if ( results[i][j].totalFinishedLawrenceTestCases < 50)
				results[i][j].avgLawrenceTime = 0;
			if ( results[i][j].totalFinishedTriangulationTestCases < 50 )
				results[i][j].avgTriangulationTime = 0;
		}
		for(int j = 0; j < (int) resultsDual[i].size(); ++j)
		{
			if ( resultsDual[i][j].totalFinishedLawrenceTestCases < 50)
				resultsDual[i][j].avgLawrenceTime = 0;
			if ( resultsDual[i][j].totalFinishedTriangulationTestCases < 50 )
				resultsDual[i][j].avgTriangulationTime = 0;
		}
	}//for i. clear any incomplete test cases.


	cout << "\\begin{sidewaystable}\n";
	cout << "\\caption{Triangulation vs Lawrence Integration on Random Polytopes in Dimension " << dim << "}\n";
	cout << "\\label{tabel:lawrence-random-integration-dim" << dim <<  "}\n";
	cout << "\\begin{tabular}{l";
	for(int i= 0; i < results[0].size(); ++i)
		cout << "rr";
	cout << "}\n";
	cout << "\\toprule \n";
	//print 1st header row.
	cout << "& \\multicolumn{" << results[0].size()*2 << "}{c}{Monomial Degree}\\\\ \n";

	for (int i = 0; i < (int) results[0].size(); ++i)
		cout << " & \\multicolumn{2}{c}{" << results[0][i].degree << "} ";
	cout << " \\\\ \n";
	for (int i = 0; i < (int) results[0].size(); ++i)
		cout << "\\cmidrule(r){" << 2*i+1+1 << "-" << 2*i+2+1<< "}";
	cout << "\n";

	//print 2nd header row.
	cout << "Vert.";
	for(int i = 0; i < (int) results[0].size(); ++i)
			cout << " & Law. & Tri. ";
	cout << "\\\\ \n";
	cout << "\\hline \n";
	//print other rows.

	printLatexTable_core(results, resultsDual, fasterStyle, slowerStyle, sameStyle, zeroStyle);


	cout << "\\bottomrule \n";
	cout << "\\end{tabular} \n";
	cout << "\\end{sidewaystable} \n";
}//printLatexIntegrationTable()


void printHTMLtableIntegration(char * dbFile, int dim)
{

	vector<vector<ValuationDBStatistics> > results, resultsDual;
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

		cout << "\n<tr><!--min max row and manual stop MS-->";
		for(int j = 0; j < (int) results[i].size(); ++j)
			cout << "\n\t<td class=\"range\">[" << results[i][j].minLawrenceTime  << ", " << results[i][j].maxLawrenceTime << "]" << (results[i][j].manuallyLimitedLawrence ? "MS":"") << "</td><td class=\"range\">[" << results[i][j].minTriangulationTime << ", " << results[i][j].maxTriangulationTime << "]" << (results[i][j].manuallyLimitedTriangulation ? "MS":"") << "</td>";
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
			cout << "\n\t<td class=\"dualRange\">[" << resultsDual[i][j].minLawrenceTime  << ", " << resultsDual[i][j].maxLawrenceTime << "]" << (resultsDual[i][j].manuallyLimitedLawrence ? "MS":"") << "</td><td class=\"dualRange\">[" << resultsDual[i][j].minTriangulationTime << ", " << resultsDual[i][j].maxTriangulationTime << "]" << (resultsDual[i][j].manuallyLimitedTriangulation? "MS":"") << "</td>";
		cout << "\n</tr>";


		cout << "\n<tr><!--dual count row-->";
		for(int j = 0; j < (int) resultsDual[i].size(); ++j)
			cout << "\n\t<td class=\"dualCountNum\">" << resultsDual[i][j].totalFinishedLawrenceTestCases  << "/" << resultsDual[i][j].totalTestCases << "</td><td class=\"dualCountNum\">" << resultsDual[i][j].totalFinishedTriangulationTestCases << "/" << resultsDual[i][j].totalTestCases << "</td>";
		cout << "\n</tr>";



	}//for each row.


	cout << "</table>" << endl;
	cout << "</body></head>" << endl;
}//printHTMLtable()


void printLatexTableVolume(const char *dbFile)
{
	cerr << "going to print latex volume table.\n";
	vector<vector<ValuationDBStatistics> > results, resultsDual;
	string fasterStyle, slowerStyle, sameStyle, zeroStyle;


	fasterStyle = " \\timeFaster";
	slowerStyle = " \\timeSlower";
	sameStyle   = " \\timeTie";
	zeroStyle   = " \\timeNotComputed";

	VolumeDB db;
	db.open(dbFile);
	resultsDual = db.getStatistics(true);
	results     = db.getStatistics(false);
	db.close();

	//check data.
	assert(results.size() == resultsDual.size());
	for(int i = 0; i < results.size(); ++i)
	{
		int dim = results[i][0].dim;

		for(int k = 0; k < results[i].size(); ++k)
		{
			assert(results[i][k].dim == dim);
			assert(results[i][k].vertexCount - dim == results[0][k].vertexCount - results[0][k].dim);
		}
	}

	cout << fixed;

	for(int i = 0; i < (int) results.size(); ++i)
	{
		for(int j = 0; j < (int) results[i].size(); ++j)
		{
			if ( results[i][j].totalFinishedLawrenceTestCases < 50)
				results[i][j].avgLawrenceTime = 0;
			if ( results[i][j].totalFinishedTriangulationTestCases < 50 )
				results[i][j].avgTriangulationTime = 0;
		}
		for(int j = 0; j < (int) resultsDual[i].size(); ++j)
		{
			if ( resultsDual[i][j].totalFinishedLawrenceTestCases < 50)
				resultsDual[i][j].avgLawrenceTime = 0;
			if ( resultsDual[i][j].totalFinishedTriangulationTestCases < 50 )
				resultsDual[i][j].avgTriangulationTime = 0;
		}
	}//for i. clear any incomplete test cases.


	cout << "\\begin{table}\n";
	cout << "\\caption{Triangulation vs Lawrence Volume on Random Polytopes}\n";
	cout << "\\label{tabel:lawrence-random-volume}\n";
	cout << "\\begin{tabular}{l";
	for(int i= 0; i < results[0].size(); ++i)
		cout << "rr";
	cout << "}\n";
	cout << "\\toprule \n";
	//print 1st header row.
	cout << "& \\multicolumn{" << results[0].size()*2 << "}{c}{Primal Vertex Count}\\\\ \n";

	for (int i = 0; i < (int) results[0].size(); ++i)
		cout << " & \\multicolumn{2}{c}{Dim+" << results[0][i].vertexCount - results[0][i].dim << "} ";
	cout << " \\\\ \n";
	for (int i = 0; i < (int) results[0].size(); ++i)
		cout << "\\cmidrule(r){" << 2*i+1+1 << "-" << 2*i+2+1<< "}";
	cout << "\n";

	//print 2nd header row.
	cout << "Dim.";
	for(int i = 0; i < (int) results[0].size(); ++i)
			cout << " & Law. & Tri. ";
	cout << "\\\\ \n";
	cout << "\\hline \n";
	//print other rows.

	printLatexTable_core(results, resultsDual, fasterStyle, slowerStyle, sameStyle, zeroStyle);


	cout << "\\bottomrule \n";
	cout << "\\end{tabular} \n";
	cout << "\\end{table} \n";
}//printLatexTableVolume

void printHTMLtableVolume(const char *dbFile)
{
	vector<vector<ValuationDBStatistics> > results, resultsDual;
	string fasterStyle, slowerStyle;
	string timeT, timeL;

	fasterStyle = " fasterTime";
	slowerStyle = " slowerTime";

	VolumeDB db;
	db.open(dbFile);
	resultsDual = db.getStatistics(true);
	results     = db.getStatistics(false);
	db.close();

	//check data.
	assert(results.size() == resultsDual.size());
	for(int i = 0; i < results.size(); ++i)
	{
		int dim = results[i][0].dim;

		for(int k = 0; k < results[i].size(); ++k)
		{
			assert(results[i][k].dim == dim);
			assert(results[i][k].vertexCount - dim == results[0][k].vertexCount - results[0][k].dim);
		}
	}

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
	cout << "row number is dim, and column is (vertexCount - dim)<br />\n";
	cout << "<table border=\"1\">" << endl;

	//print 1st header row.
	cout << "\n<tr><td>x</td>";
	for(int i = 0; i < (int) results[0].size(); ++i)
		cout << "\n\t<td colspan=\"2\">dim+1+" << results[0][i].vertexCount  - results[0][i].dim - 1<< "</td>";
	cout << "\n</tr>" << endl;

	//print 2nd header row.
	cout << "\n<tr><td>x</td>";
	for(int i = 0; i < (int) results[0].size(); ++i)
			cout << "\n\t<td>Law.</td><td>Tri.</td>";
	cout << "\n</tr>" << endl;

	//print other rows.
	for(int i = 0; i < (int) results.size(); ++i)
	{
		cout << "\n<tr><!--start of dim " << results[i][0].dim << "-->";

		//print 1st col.
		cout << "\n\t<td rowspan=\"6\">" << results[i][0].dim << "</td>";

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

			cout << "\n\t<td class=\"time" << timeL.c_str() << "\">" << results[i][j].avgLawrenceTime << "</td><!--dim:" << results[i][j].dim << ", vc:" << results[i][j].vertexCount << ", dual:" << results[i][j].useDual<< "--><td class=\"time" << timeT.c_str() << "\">" << results[i][j].avgTriangulationTime << "</td>";
		}
		cout << "\n</tr>";

		cout << "\n<tr><!--min max row-->";
		for(int j = 0; j < (int) results[i].size(); ++j)
			cout << "\n\t<td class=\"range\">[" << results[i][j].minLawrenceTime  << ", " << results[i][j].maxLawrenceTime << "]" << (results[i][j].manuallyLimitedLawrence? "MS":"") << "</td><td class=\"range\">[" << results[i][j].minTriangulationTime << ", " << results[i][j].maxTriangulationTime << "]" << (results[i][j].manuallyLimitedTriangulation? "MS":"") << "</td>";
		cout << "\n</tr>";


		cout << "\n<tr><!--count row-->";
		for(int j = 0; j < (int) results[i].size(); ++j)
			cout << "\n\t<td class=\"countNum\">" << results[i][j].totalFinishedLawrenceTestCases  << "/" << results[i][j].totalTestCases << "</td><td class=\"countNum\">" << results[i][j].totalFinishedTriangulationTestCases << "/" << results[i][j].totalTestCases << "</td>";
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
			cout << "\n\t<td class=\"dualRange\">[" << resultsDual[i][j].minLawrenceTime  << ", " << resultsDual[i][j].maxLawrenceTime << "]" << (resultsDual[i][j].manuallyLimitedLawrence? "MS":"") << "</td><td class=\"dualRange\">[" << resultsDual[i][j].minTriangulationTime << ", " << resultsDual[i][j].maxTriangulationTime << "]" << (resultsDual[i][j].manuallyLimitedTriangulation? "MS":"") << "</td>";
		cout << "\n</tr>";


		cout << "\n<tr><!--dual count row-->";
		for(int j = 0; j < (int) resultsDual[i].size(); ++j)
			cout << "\n\t<td class=\"dualCountNum\">" << resultsDual[i][j].totalFinishedLawrenceTestCases  << "/" << resultsDual[i][j].totalTestCases << "</td><td class=\"dualCountNum\">" << resultsDual[i][j].totalFinishedTriangulationTestCases << "/" << resultsDual[i][j].totalTestCases << "</td>";
		cout << "\n</tr>";



	}//for each row.

	cout << "</table>" << endl;
	cout << "</body></head>" << endl;

}//printHTMLtableVolume



void printMatlabTableVolume(const char * dbFile, const char * mFile)
{
	vector<vector<ValuationDBStatistics> > results, resultsDual;
	string timeT, timeL;

	VolumeDB db;
	db.open(dbFile);
	resultsDual = db.getStatistics(true);
	results     = db.getStatistics(false);
	db.close();

	//check data.
	assert(results.size() == resultsDual.size());
	for(int i = 0; i < results.size(); ++i)
	{
		int dim = results[i][0].dim;

		for(int k = 0; k < results[i].size(); ++k)
		{
			assert(results[i][k].dim == dim);
			assert(results[i][k].vertexCount - dim == results[0][k].vertexCount - results[0][k].dim);
		}
	}

	//open the file.
	ofstream file(mFile);

	//make dim vector
	file << "dimArray = [";
	for(int i = 0; i < (int) results.size(); ++i)
		file << results[i][0].dim << ' ';
	file << "]\n\n";

	//make point vector
	file << "pointArray = [";
	for(int i = 0; i < (int) results[0].size(); ++i)
		file << results[0][i].vertexCount - results[0][i].dim << ' ';
	file << "]\n\n";

	//make primal volume matrix lawrence.
	file << "primalVolumeLawrence = [";
	for(int i = 0; i < (int) results.size(); ++i)
	{
		for(int j = 0; j < (int) results[i].size(); ++j)
			file << results[i][j].avgLawrenceTime << ' ';
		file << ";\n";
	}
	file << "]\n\n";

	//make primal volume matrix triang.
	file << "primalVolumeTriangulation = [";
	for(int i = 0; i < (int) results.size(); ++i)
	{
		for(int j = 0; j < (int) results[i].size(); ++j)
			file << results[i][j].avgTriangulationTime << ' ';
		file << ";\n";
	}
	file << "]\n\n";

	//make dual volume matrix lawrence.
	file << "dualVolumeLawrence = [";
	for(int i = 0; i < (int) resultsDual.size(); ++i)
	{
		for(int j = 0; j < (int) resultsDual[i].size(); ++j)
			file << resultsDual[i][j].avgLawrenceTime << ' ';
		file << ";\n";
	}
	file << "]\n\n";

	//make dual volume matrix triang.
	file << "dualVolumeTriangulation = [";
	for(int i = 0; i < (int) resultsDual.size(); ++i)
	{
		for(int j = 0; j < (int) resultsDual[i].size(); ++j)
			file << resultsDual[i][j].avgTriangulationTime << ' ';
		file << ";\n";
	}
	file << "]\n\n";

	file.close();

}//printMatlabTavleVolume.


void printLatexTableIntegrationPerPolymakeFile(const char *dbFile)
{

	vector<vector<ValuationDBStatistics> > results, resultsDual;
	vector<vector<string> > allPolymakeFiles;
	string fasterStyle, slowerStyle, sameStyle, zeroStyle;
	string timeT, timeL;

	fasterStyle = " \\timeFaster";
	slowerStyle = " \\timeSlower";
	sameStyle   = " \\timeTie";
	zeroStyle   = " \\timeNotComputed";

	IntegrationDB db;
	db.open(dbFile);
	allPolymakeFiles = db.getAllPolymakeFiles();
	resultsDual = db.getStatisticsByFile(allPolymakeFiles, true);
	results     = db.getStatisticsByFile(allPolymakeFiles, false);
	db.close();
	assert(results.size() == resultsDual.size());

	cout << fixed;

	for(int i = 0; i < (int) results.size(); ++i)
	{
		for(int j = 0; j < (int) results[i].size(); ++j)
		{
			if ( results[i][j].totalFinishedLawrenceTestCases < 50)
				results[i][j].avgLawrenceTime = 0;
			if ( results[i][j].totalFinishedTriangulationTestCases < 50 )
				results[i][j].avgTriangulationTime = 0;
		}
		for(int j = 0; j < (int) resultsDual[i].size(); ++j)
		{
			if ( resultsDual[i][j].totalFinishedLawrenceTestCases < 50)
				resultsDual[i][j].avgLawrenceTime = 0;
			if ( resultsDual[i][j].totalFinishedTriangulationTestCases < 50 )
				resultsDual[i][j].avgTriangulationTime = 0;
		}
	}//for i. clear any incomplete test cases.


	cout << "\\begin{table}\n";
	cout << "\\caption{Triangulation vs Lawrence Integration on the Ziegler Database}\n";
	cout << "\\label{tabel:lawrence-ziegler-integration}\n";
	cout << "\\begin{tabular}{l";
	for(int i= 0; i < results[0].size(); ++i)
		cout << "rr";
	cout << "}\n";
	cout << "\\toprule \n";
	//print 1st header row.
	cout << "& \\multicolumn{" << results[0].size()*2 << "}{c}{Monomial Degree}\\\\ \n";

	for (int i = 0; i < (int) results[0].size(); ++i)
		cout << " & \\multicolumn{2}{c}{" << results[0][i].degree << "} ";
	cout << " \\\\ \n";
	for (int i = 0; i < (int) results[0].size(); ++i)
		cout << "\\cmidrule(r){" << 2*i+1+1 << "-" << 2*i+2+1<< "}";
	cout << "\n";

	//print 2nd header row.
	cout << "Vert.";
	for(int i = 0; i < (int) results[0].size(); ++i)
			cout << " & Law. & Tri. ";
	cout << "\\\\ \n";
	cout << "\\hline \n";
	//print other rows.
	for(int i = 0; i < (int) results.size(); ++i)
	{
		//print 1st col.
		cout << allPolymakeFiles[i][0].c_str();

		for(int j = 0; j < (int) results[i].size(); ++j)
		{
			if ( results[i][j].avgLawrenceTime == 0 )
			{
				if ( results[i][j].avgTriangulationTime != 0)
				{
					timeT = fasterStyle;
					timeL = zeroStyle;
				}
				else
				{
					timeT = zeroStyle;
					timeL = zeroStyle;
				}
			}
			else if ( results[i][j].avgTriangulationTime == 0)
			{
				timeT = zeroStyle;
				timeL = fasterStyle; //not zero.
			}
			else if (results[i][j].avgLawrenceTime < results[i][j].avgTriangulationTime)
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
				timeL = sameStyle;
				timeT = sameStyle;
			}//tie.

			cout << " & " << timeL.c_str() << "{ " << setw(6) << setprecision(2) << results[i][j].avgLawrenceTime << "}"
				 << " & " << timeT.c_str() << "{ " << setw(6) << setprecision(2) << results[i][j].avgTriangulationTime << "}";
		}
		cout << "\\\\ \n";

		//dual rows.
		for(int j = 0; j < (int) resultsDual[i].size(); ++j)
		{
			if ( resultsDual[i][j].avgLawrenceTime == 0 )
			{
				if ( resultsDual[i][j].avgTriangulationTime != 0)
				{
					timeT = fasterStyle;
					timeL = zeroStyle;
				}
				else
				{
					timeT = zeroStyle;
					timeL = zeroStyle;
				}
			}
			else if ( resultsDual[i][j].avgTriangulationTime == 0)
			{
				timeT = zeroStyle;
				timeL = fasterStyle; //not zero.
			}
			else if (resultsDual[i][j].avgLawrenceTime < resultsDual[i][j].avgTriangulationTime)
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
				timeL = sameStyle;
				timeT = sameStyle;
			}//tie.

			cout << " & " << timeL.c_str() << "{ " << setw(6) << setprecision(2) << resultsDual[i][j].avgLawrenceTime << "}"
				 << " & " << timeT.c_str() << "{ " << setw(6) << setprecision(2) << resultsDual[i][j].avgTriangulationTime << "}";
		}
		cout << "\\\\ \n";
		cout << "\\hline \n";
	}//for each row.

	cout << "\\bottomrule \n";
	cout << "\\end{tabular} \n";
	cout << "\\end{table} \n";
}//printLatexTablePerPolymakeFile


void printHTMLtableIntegrationPerPolymakeFile(const char *dbFile)
{

	vector<vector<ValuationDBStatistics> > results, resultsDual;
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
			cout << "\n\t<td class=\"range\">[" << results[i][j].minLawrenceTime  << ", " << results[i][j].maxLawrenceTime << "]" << (results[i][j].manuallyLimitedLawrence? "MS":"") << "</td><td class=\"range\">[" << results[i][j].minTriangulationTime << ", " << results[i][j].maxTriangulationTime << "]" << (results[i][j].manuallyLimitedTriangulation? "MS":"") << "</td>";
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
			cout << "\n\t<td class=\"dualRange\">[" << resultsDual[i][j].minLawrenceTime  << ", " << resultsDual[i][j].maxLawrenceTime << "]" << (resultsDual[i][j].manuallyLimitedLawrence? "MS":"") << "</td><td class=\"dualRange\">[" << resultsDual[i][j].minTriangulationTime << ", " << resultsDual[i][j].maxTriangulationTime << "]" << (resultsDual[i][j].manuallyLimitedTriangulation? "MS":"") << "</td>";
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
	if (argc < 3 )
	{
		cerr << "Hello,\n "
			 << "This program will print time tables for easy viewing"
			 << endl;

		cerr << "error. usage:\n";
		cerr << "ex: \"sqlite-db volume [html | latex]\" prints the volume tables" << endl;
		//cerr << "ex: \"sqlite-db volume \" prints the volume table per file in latex." << endl;
		//cerr << "ex: \"sqlite-db volume matlab\" prints the volume mesh to a matlab file " << argv[0] << "-volume.m " << endl;
		cerr << "ex: \"sqlite-db integration\" prints the integration db on the specfic file db" << endl;
		cerr << "ex: \"sqlite-db integration [html | latex ] 5\" prints the integration table for dim 5 only" << endl;
		exit(1);
	}

	for(int i = 0; i < argc; ++i)
		cerr << "arg " << i << " = " << argv[i] << endl;

	if ( strcmp(argv[2], "volume") == 0)
	{
		if ( 0 == strcmp(argv[3], "html"))
			printHTMLtableVolume(argv[1]); //exe db volume
		else if ( 0 == strcmp(argv[3], "latex"))
			printLatexTableVolume(argv[1]); //exe db volume
		else
			cerr << "volume table style " << argv[3] << "not known" << endl;
	}
	else if (strcmp(argv[2], "integration") == 0)
	{
		if ( argc == 5)
		{
			if ( 0 == strcmp(argv[3], "html"))
				printHTMLtableIntegration(argv[1], atoi(argv[4])); //exe db integration html dim
			else if ( 0 == strcmp(argv[3], "latex"))
				printLatexTableIntegration_sideways(argv[1], atoi(argv[4])); //exe db integration latex dim
		}
		else
			if ( 0 == strcmp(argv[3], "html"))
				printHTMLtableIntegrationPerPolymakeFile(argv[1]); //exe db integration
			else if ( 0 == strcmp(argv[3], "latex"))
				printLatexTableIntegrationPerPolymakeFile(argv[1]); //exe db integration
			else
				cerr << "table style " << argv[3] << " not known" << endl;
	}//print an integration table.
	else
		cerr << "valuation " << argv[2] << "is not known" << endl;



	return 0;
}//main



