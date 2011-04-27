
#include "IntegrationDB.h"
#include <sstream>
#include <cstdlib> //for atoi
#include <cassert>
#include <iostream>
IntegrationDB::IntegrationDB()
{
}


/**
 * @parm seconds: time limit how how each polytope should take.
 * @parm dim: of polytope
 * @parm vertexCount of the polytope (or it's dual if useDual is ture)
 * @parm degree of the polynomial
 * @parm useDual: true if we are testing the dual polytopes
 * @return bool: true if the test case "dim-vertecCount-degree-dual" should finish with secondsLimit
 */
bool IntegrationDB::canTestFinish(AlgorithemUsed alg, int dim, int vertexCount, int degree, bool useDual, int secondsLimit)
{
	stringstream sql;
	vector<vector<string> > result;
	//if we have data for this test case, great, find the average.

	if ( useDual == true)
	{
		sql << "select count(*)"
			<< " from polynomial as p, integrate as i, polytope as dualP, polytope as orgP "
			<< " where i.polytopeID = dualP.rowid and i.polynomialID = p.rowid " //joint dual-integrate-polynomial
			<< " and dualP.dual = orgP.rowid and dualP.dual is not null" //join the dual and org. polytope.
			<< " and orgP.dim = " << dim
			<< " and orgP.vertexCount = " << vertexCount
			<< " and p.degree <= " << degree
			<< " and " << (alg == Lawrence ? " i.timeLawrence " : " i.timeTriangulate ") << " = -2";//-2 = I manually said skip this test case

	}
	else
	{
		sql << "select count(*)"
			<< " from polynomial as p, integrate as i, polytope as t "
			<< " where i.polytopeID = t.rowid and i.polynomialID = p.rowid "
			<< " and t.dim = " << dim
			<< " and t.vertexCount = " << vertexCount
			<< " and p.degree <= " << degree
			<< " and " << (alg == Lawrence ? " i.timeLawrence " : " i.timeTriangulate ") << " = -2";
	}//else not using the dual.

	if ( queryAsInteger(sql.str().c_str()) )
		return false; //skip this test case.

	sql.str("");
	if ( useDual == true)
	{
		sql << "select avg(" << (alg == Lawrence ? " i.timeLawrence " : " i.timeTriangulate ") << "), count(*)"
			<< " from polynomial as p, integrate as i, polytope as dualP, polytope as orgP "
			<< " where i.polytopeID = dualP.rowid and i.polynomialID = p.rowid " //joint dual-integrate-polynomial
			<< " and dualP.dual = orgP.rowid and dualP.dual is not null" //join the dual and org. polytope.
			<< " and orgP.dim = " << dim
			<< " and orgP.vertexCount = " << vertexCount
			<< " and p.degree = " << degree
			<< " and " << (alg == Lawrence ? " i.timeLawrence " : " i.timeTriangulate ") << " >= 0";

	}
	else
	{
		sql << "select avg(" << (alg == Lawrence ? " i.timeLawrence " : " i.timeTriangulate ") << "), count(*)"
			<< " from polynomial as p, integrate as i, polytope as t "
			<< " where i.polytopeID = t.rowid and i.polynomialID = p.rowid "
			<< " and t.dim = " << dim
			<< " and t.vertexCount = " << vertexCount
			<< " and p.degree = " << degree
			<< " and " << (alg == Lawrence ? " i.timeLawrence " : " i.timeTriangulate ") << " >= 0";
	}//else not using the dual.

	result = query(sql.str().c_str());
	if ( result[0][0] != "NULL" && atoi(result[0][1].c_str()) >= 1 )
		return (atof(result[0][0].c_str()) < secondsLimit); //if there are at least 1 test case

	//darn, if we are here, then there where no (or less than 3) test cases of this class.
	//now look at the previous tests and see if they did not finish.
	if (vertexCount <= dim + 3)
		return true; //always do the basic cases!

	//find the max time for all vertexCount and polynomial degree less than the current setting for fix dim.
	sql.str("");
	if ( useDual == true)
	{
		sql << "select max(" << (alg == Lawrence ? " i.timeLawrence " : " i.timeTriangulate ") << ")"
			<< " from polynomial as p, integrate as i, polytope as dualP, polytope as orgP "
			<< " where i.polytopeID = dualP.rowid and i.polynomialID = p.rowid " //joint dual-integrate-polynomial
			<< " and dualP.dual = orgP.rowid and dualP.dual is not null" //join the dual and org. polytope.
			<< " and orgP.dim = " << dim
			<< " and orgP.vertexCount = " << vertexCount //the difference here with the above is the
			<< " and p.degree <= " << degree              // <=  on degree.
			<< " and " << (alg == Lawrence ? " i.timeLawrence " : " i.timeTriangulate ") << " >= 0";

	}
	else
	{
		sql << "select max(" << (alg == Lawrence ? " i.timeLawrence " : " i.timeTriangulate ") << ")"
			<< " from polynomial as p, integrate as i, polytope as t "
			<< " where i.polytopeID = t.rowid and i.polynomialID = p.rowid "
			<< " and t.dim = " << dim
			<< " and t.vertexCount = " << vertexCount
			<< " and p.degree <= " << degree
			<< " and " << (alg == Lawrence ? " i.timeLawrence " : " i.timeTriangulate ") << " >= 0";
	}//else not using the dual.
	result = query(sql.str().c_str());
	if ( result[0][0] == "NULL")
		return true; //we have no data to say one way or the other.


	return (atof(result[0][0].c_str()) < secondsLimit);

}//canTestFinish


/**
 * @polymakeFile a .polymakefile, not a .dual.polymake file
 * @degree: of the polynomial to test
 * @return bool. true = we thing we can finish this file on the current degree.
 */
bool IntegrationDB::canSpecficFileFinish(AlgorithemUsed alg, const char *polymakeFile, int degree, int useDual, int secondsLimit)
{
	stringstream sql;
	vector<vector<string> > result;

	//fist check that I did not manually want to skip this test case.
	if ( useDual == true)
	{
		sql << "select count(*)"
			<< " from polynomial as p, integrate as i, polytope as dualP, polytope as orgP "
			<< " where i.polytopeID = dualP.rowid and i.polynomialID = p.rowid " //joint dual-integrate-polynomial
			<< " and dualP.dual = orgP.rowid and dualP.dual is not null" //join the dual and org. polytope.
			<< " and orgP.polymakeFilePath = '" << polymakeFile << "'"
			<< " and p.degree <= " << degree
			<< " and " << (alg == Lawrence ? " i.timeLawrence " : " i.timeTriangulate ") << " = -2"; // -2 means "skip this"
	}
	else
	{
		sql << "select count(*)"
			<< " from polynomial as p, integrate as i, polytope as t "
			<< " where i.polytopeID = t.rowid and i.polynomialID = p.rowid "
			<< " and t.polymakeFilePath = '" << polymakeFile << "'"
			<< " and p.degree <= " << degree
			<< " and " << (alg == Lawrence ? " i.timeLawrence " : " i.timeTriangulate ") << " = -2";
	}//else not using the dual.

	//cout << "check -2: " << sql.str().c_str() << "\n" << endl;
	result = query(sql.str().c_str());
	//cout << "check-2 ans: " << result[0][0].c_str() << '\n' << endl;
	if ( atoi(result[0][0].c_str()) )
		return false; //skip this test case.

	sql.str("");
	//find the average on the set degree.
	if ( useDual == true)
	{
		sql << "select avg(" << (alg == Lawrence ? " i.timeLawrence " : " i.timeTriangulate ") << ")"
			<< " from polynomial as p, integrate as i, polytope as dualP, polytope as orgP "
			<< " where i.polytopeID = dualP.rowid and i.polynomialID = p.rowid " //joint dual-integrate-polynomial
			<< " and dualP.dual = orgP.rowid and dualP.dual is not null" //join the dual and org. polytope.
			<< " and orgP.polymakeFilePath = '" << polymakeFile << "'"
			<< " and p.degree = " << degree // <=  on degree.
			<< " and " << (alg == Lawrence ? " i.timeLawrence " : " i.timeTriangulate ") << " >= 0";
	}
	else
	{
		sql << "select avg(" << (alg == Lawrence ? " i.timeLawrence " : " i.timeTriangulate ") << ")"
			<< " from polynomial as p, integrate as i, polytope as t "
			<< " where i.polytopeID = t.rowid and i.polynomialID = p.rowid "
			<< " and t.polymakeFilePath = '" << polymakeFile << "'"
			<< " and p.degree = " << degree
			<< " and " << (alg == Lawrence ? " i.timeLawrence " : " i.timeTriangulate ") << " >= 0";
	}//else not using the dual.

	//cout << "get avg: " << sql.str().c_str() << '\n' << endl;
	result = query(sql.str().c_str());
	//cout << "get avg: ans " << result[0][0].c_str() << endl;
	/*
	if ( !strcmp(polymakeFile, "./Various/3simp3simp.polymake") && degree == 50)
	{
		cout << "BEFORE file: " << polymakeFile
			 << "\ndegree: " << degree
			 << "\ndual: " << useDual
			 << "\ntime: " << (result[0][0].c_str())
			 << "\n" << (atof(result[0][0].c_str()) < secondsLimit) << endl;

		if ( result[0][0] != "NULL")
			exit(1);
	}*/


	if ( result[0][0] != "NULL")
		return (atof(result[0][0].c_str()) < secondsLimit);

	//if we got here, we never saw this file with this degree before.
	//now look at the previous degrees and see what was the last degree it finished at.
	if (degree <= 3)
		return true; //always do the basic cases!

	//find the max of any degree done.
	sql.str("");
	if ( useDual == true)
	{
		sql << "select max(" << (alg == Lawrence ? " i.timeLawrence " : " i.timeTriangulate ") << ")"
			<< " from polynomial as p, integrate as i, polytope as dualP, polytope as orgP "
			<< " where i.polytopeID = dualP.rowid and i.polynomialID = p.rowid " //joint dual-integrate-polynomial
			<< " and dualP.dual = orgP.rowid and dualP.dual is not null" //join the dual and org. polytope.
			<< " and orgP.polymakeFilePath = '" << polymakeFile << "'"
			<< " and p.degree <= " << degree // <=  on degree.
			<< " and " << (alg == Lawrence ? " i.timeLawrence " : " i.timeTriangulate ") << " >= 0";
	}
	else
	{
		sql << "select max(" << (alg == Lawrence ? " i.timeLawrence " : " i.timeTriangulate ") << ")"
			<< " from polynomial as p, integrate as i, polytope as t "
			<< " where i.polytopeID = t.rowid and i.polynomialID = p.rowid "
			<< " and t.polymakeFilePath = '" << polymakeFile << "'"
			<< " and p.degree <= " << degree
			<< " and " << (alg == Lawrence ? " i.timeLawrence " : " i.timeTriangulate ") << " >= 0";
	}//else not using the dual.


	//cout << "avg anything: " << sql.str().c_str() << '\n' << endl;
	result = query(sql.str().c_str());

	/*
	if ( !strcmp(polymakeFile, "./Various/3simp3simp.polymake") && degree == 50)
	{
		cout << "AFTER file: " << polymakeFile
			 << "\ndegree: " << degree
			 << "\ndual: " << useDual
			 << "\ntime: " << (result[0][0].c_str())
			 << "\n"<< (atof(result[0][0].c_str()) < secondsLimit) << endl;
		exit(1);
	}*/
	if ( result[0][0] == "NULL")
		return true; //we have no data to say one way or the other.

	return (atof(result[0][0].c_str()) < secondsLimit);
}//canSpecficTestFinish

void IntegrationDB::deletePolynomial(int id)
{
	stringstream sql;
	sql << "delete from polynomial where rowid = " << id;
	query(sql.str().c_str());
}
/**
 * @parm p: polynomial filename to search by
 * @return int: rowid of matching filename or 0.
 * We assume all the rowid's are positive ( you never manually insert a non-positive!)
 */
int IntegrationDB::doesPolynomialExist(const char * p)
{
	stringstream sql;
	sql << "select rowid from polynomial where filePath = '" << p << "' limit 1";
	vector<vector<string> > result = query(sql.str().c_str());

	if (result.size())
	{
		return atoi(result[0][0].c_str());
	}
	return 0;
}//doesPolynomialExist

/**
 * Searches for polytope rowid by polymake file path.
 */
int IntegrationDB::doesPolytopeExist(const char *polymake)
{
	stringstream sql;
	sql << "select rowid from polytope where polymakeFilePath = '" << polymake << "'";

	vector<vector<string> > result = query(sql.str().c_str());

	if (result.size())
	{
		return atoi(result[0][0].c_str());
	}
	return 0;
}//doesPolytopeExist

/**
 * Return all the non-dual polymake files.
 */
vector<vector<string> > IntegrationDB::getAllPolymakeFiles()
{
	stringstream sql;
	sql << "select distinct polymakeFilePath from polytope where dual is null order by dim, vertexCount asc";
	return query(sql.str().c_str());
}//getAllPolymakeFiles

/*
 * Returns true if this case contained a "-2" as the time value for a test.
 */
bool IntegrationDB::getLimit(AlgorithemUsed alg, int dim, int vertexCount, int degree, bool useDual)
{
	stringstream sql;
	string strAlg;

	if (alg == Lawrence)
		strAlg = "timeLawrence";
	else
		strAlg = "timeTriangulate";


	if (useDual == true)
	{
		sql << "select count(i." << strAlg << ")"
			<< " from integrate as i"
			<< " join polytope as dualP on dualP.rowid = i.polytopeID"
			<< " join polytope as orgP on orgP.rowid = dualP.dual"
			<< " join polynomial as p on p.rowid = i.polynomialID"
			<< " where dualP.dual is not null " //and with orgp
			<< " and orgP.dim = " << dim
			<< " and orgP.vertexCount = " << vertexCount
			<< " and p.degree = " << degree
			<< " and i." << strAlg << " = -2 ";
	}//if dual
	else
	{

		sql << "select avg(i." << strAlg << ") "
			<< " from integrate as i"
			<< " join polynomial as p on p.rowid = i.polynomialID"
			<< " join polytope as t on t.rowid = i.polytopeID"
			<< " where t.dual is null"
			<< " and t.dim = " << dim
			<< " and t.vertexCount = " << vertexCount
			<< " and p.degree = " << degree
			<< " and i." << strAlg << " = -2 ";
	}//regular

	return (queryAsInteger(sql.str().c_str()) > 0);
}

/*
 * Returns true if this test case as a -2 as a time for a test result.
 */
bool IntegrationDB::getLimitByFile(AlgorithemUsed alg, const string &polymakeFile, int degree, bool useDual)
{
	stringstream sql;
	string strAlg;

	strAlg = (alg == Lawrence ? "timeLawrence" : "timeTriangulate");

	if (useDual == true)
	{
		sql << "select count(i." << strAlg << ") "
			<< " from polynomial as p, polytope as dualP, polytope as orgP, integrate as i "
			<< " where i.polynomialID = p.rowid and i.polytopeID = dualP.rowid " //join with p, i, dualP
			<< " and dualP.dual is not null and dualP.dual = orgP.rowid" //and with orgp
			<< " and orgP.polymakeFilePath = '" << polymakeFile << "'"
			<< " and p.degree = " << degree
			<< " and i." << strAlg << " = -2 ";
	}//if dual
	else
	{

		sql << "select count(i." << strAlg << ") "
			<< " from polynomial as p, polytope as t, integrate as i "
			<< " where i.polynomialID = p.rowid and i.polytopeID = t.rowid " //join with p, i, dualP
			<< " and t.dual is null"
			<< " and t.polymakeFilePath = '" << polymakeFile << "'"
			<< " and p.degree = " << degree
			<< " and i." << strAlg << " = -2 ";
	}//regular

	return (queryAsInteger(sql.str().c_str()) > 0);
}


/**
 * Returns the number of polynomials in the polynomials table with this dim and degree
 */
int IntegrationDB::getNumberPolynomials(int dim, int degree)
{
	stringstream sql;
	sql << "select count(*) from polynomial where dim = " << dim << " and degree = " << degree;
	return queryAsInteger(sql.str().c_str());
}//getNumberPolynomials

/**
 * Return the number of polytopes with set dim and vertexCount and dual value.
 */
int IntegrationDB::getNumberPolytopes(int dim, int vertexCount, bool useDuals)
{
	stringstream sql;
	if (useDuals == false)
		sql << "select count(*) from polytope where dim = " << dim
			<< " and vertexCount = " << vertexCount
			<< " and dual is  null";
	else
		sql << "select count(*) from polytope as orgP, polytope as dualP"
			<< " where orgP.dim = " << dim
			<< "  and orgP.rowid = dualP.dual"
			<< "  and orgP.vertexCount = " << vertexCount;

	cout << sql.str().c_str() << endl;
	//when the column for dual is null, then the polytope is not a dual.
	//when the column for dual is NOT null, then the polytope is a dual.
	return queryAsInteger(sql.str().c_str());
}//getNumberPolytopes

/**
 * Returns the number of integration tests that have this dim, vertex count, and polynomial degree (and is or isn't a dual)
 */
int IntegrationDB::getNumberIntegrationTest(int dim, int vertexCount, int degree, bool useDuals)
{
	stringstream sql; //we are going to join the tables.
	if ( useDuals == false)
	{
		sql << "select count(*) from polynomial as p, polytope as t, integrate as i "
			<< "where i.polynomialID = p.rowid and i.polytopeID = t.rowid " //join
			<< " and t.dim = " << dim << " and t.vertexCount = " << vertexCount
			<< " and p.degree = " << degree
			<< " and t.dual is null ";
	}
	else
	{
		sql << "select count(*) from polynomial as p, polytope as orgP, polytope as dualP, integrate as i "
				<< "where i.polynomialID = p.rowid and i.polytopeID = dualP.rowid " //join
				<< " and dualP.dual = orgP.rowid "
				<< " and orgP.dim = " << dim << " and orgP.vertexCount = " << vertexCount
				<< " and p.degree = " << degree
				<< " and dualP.dual is not null ";
	}
	return queryAsInteger(sql.str().c_str());
}//getNumberIntegrationTest


int IntegrationDB::getNumberIntegrationTest(const string &polymakeFile, int degree, bool useDual)
{
	stringstream sql; //we are going to join the tables.
	if ( useDual == false)
	{
		sql << "select count(*) from polynomial as p, polytope as t, integrate as i "
			<< "where i.polynomialID = p.rowid and i.polytopeID = t.rowid " //join
			<< " and t.polymakeFilePath = '" << polymakeFile << "'"
			<< " and p.degree = " << degree
			<< " and t.dual is null ";
	}
	else
	{
		sql << "select count(*) from polynomial as p, polytope as orgP, polytope as dualP, integrate as i "
				<< "where i.polynomialID = p.rowid and i.polytopeID = dualP.rowid " //join
				<< " and dualP.dual = orgP.rowid "
				<< " and orgP.polymakeFilePath = '" << polymakeFile << "'"
				<< " and p.degree = " << degree
				<< " and dualP.dual is not null ";
	}
	return queryAsInteger(sql.str().c_str());
}//getNumberIntegrationTest

/**
 * @polytopeID: rowid of a polytope.
 * @degree: of a polynomial
 * @return: number of integration tests in the integrate table with this polytope and a polynomial of this degree.
 */
int IntegrationDB::getNumberIntegrationTest(int polytopeID, int degree)
{
	stringstream sql;
	sql << "select count(*) from integrate as i, polytope as t, polynomial as p"
		<< " where i.polytopeID = t.rowid and i.polynomialID = p.rowid " //join
		<< " and t.rowid = " << polytopeID
		<< " and p.degree = " << degree;

	cout << sql.str().c_str() << endl;
	return queryAsInteger(sql.str().c_str());
}//getNumberIntegrationTest


/**
 * Gets all the integrate rows that have a set dim, vertex count, degree, and dual values.
 */
vector<vector<string> > IntegrationDB::getRowsToIntegrate(int dim, int vertex, int degree, bool useDual, int limit)
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
			<< " order by t.latteFilePath, p.degree"
			<< " limit " << limit;
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
			<< " order by t.latteFilePath, p.degree"
			<< " limit " << limit;
	}
	return query(sql.str().c_str());
}//getRowsToIntegrate

/**
 * Given a org. polymake file, find all the test (or it's dual) with a set degree
 */
vector<vector<string> > IntegrationDB::getRowsToIntegrateGivenSpecficFile(char *polymakeFile, int degree, bool useDual, int limit)
{
	stringstream sql;
	if (useDual == false)
	{
		sql << "select p.filePath, t.latteFilePath, i.timeLawrence, i.timeTriangulate, i.integral, i.rowid"
			<< " from polynomial as p, polytope as t, integrate as i "
			<< " where p.rowid = i.polynomialID and t.rowid = i.polytopeID "
			<< " and p.degree = " << degree
			<< " and t.dual is null"
			<< " and t.polymakeFilePath = '" << polymakeFile << "'"
			<< " limit " << limit;
	}
	else
	{
		sql << "select p.filePath, dualP.latteFilePath, i.timeLawrence, i.timeTriangulate, i.integral, i.rowid"
			<< " from polynomial as p, polytope as dualP, polytope as orgP, integrate as i "
			<< " where p.rowid = i.polynomialID and dualP.rowid = i.polytopeID "
			<< " and orgP.rowid = dualP.dual"
			<< " and p.degree = " << degree
			<< " and dualP.dual is not null"
			<< " and orgP.polymakeFilePath = '" << polymakeFile << "'"
			<< " limit " << limit;
	}
	//cout << "getRowsToIntegrateGivenSpecficFile:: " << sql.str().c_str() << endl;
	return query(sql.str().c_str());

}//getRowsToIntegrateGivenSpecficFile

/**
 * Given a set dim, will find what polynomial degrees and vertex counts exist in the db
 * and give the results/stats in a 2d matrix.
 *
 * Assumes the db only contains a few different polynomial degrees and vertex counts.
 * (That is, we didn't blindly insert non-dual polynomials of any degree, only of degree 2,3,4,5,15,20,30,...etc
 * If we did, then we would get a new column in the matrix for each different degree!!!)
 *
 * @parm: dim. dimension of polytope.
 * @useDua: tue if we want to use he dual polytopes.
 * @return: a 2d matrix of polytope test stats.
 * answer[i][k] contains information on polytope vertex-count case i and
 *   polynomial degree case k.
 */
vector<vector<ValuationDBStatistics> >  IntegrationDB::getStatisticsByDim(int dim, bool useDual)
{
	vector<vector<ValuationDBStatistics> > ans;
	vector<vector<string> > polynomialDegrees;
	vector<vector<string> > vertexCounts;
	stringstream sql;

	sql << "select distinct p.degree from polynomial as p where p.dim = " << dim;
	polynomialDegrees = query(sql.str().c_str());


	sql.str("");
	sql << "select distinct t.vertexCount from polytope as t "
		<< " where t.dual is null"
		<< " and t.dim = " << dim;
	vertexCounts = query(sql.str().c_str());

	ans.resize(vertexCounts.size()); //make room for each row.
	//now loop over every vertexCount and polynomial degree and get
	//the statistics for the class dim-vertecCount-degree-dual
	cerr << "table is " << vertexCounts.size() << " by " << polynomialDegrees.size();
	for(int row = 0; row < (int)vertexCounts.size(); ++row)
	{

		for(int col = 0; col < (int)polynomialDegrees.size(); ++col)
		{
			cerr << "row, col=" << row << ", " << col << endl;
			ans[row].push_back(getStatisticsByDimVertexDegree(dim, atoi(vertexCounts[row][0].c_str()), atoi(polynomialDegrees[col][0].c_str()), useDual));
		}//for
	}//for row. vertexCounts

	return ans;
}//getResultsByDim


vector<vector<ValuationDBStatistics> > IntegrationDB::getStatisticsByFile(const vector<vector<string> > &polymakeFile, bool useDual)
{
	vector<vector<ValuationDBStatistics> > ans;
	vector<vector<string> > polynomialDegrees;
	vector<vector<string> > vertexCounts;
	stringstream sql;

	if ( polymakeFile.size() == 0)
		return ans; //return it empty

	//get list of different degrees.
	sql << "select distinct p.degree from polynomial as p";
	polynomialDegrees = query(sql.str().c_str());


	ans.resize(polymakeFile.size()); //make room for each row.

	//now loop over every file and polynomial degree and get
	//the statistics for the class polymakefile-degree-dual
	cerr << "table is " << vertexCounts.size() << " by " << polynomialDegrees.size();
	for(int row = 0; row < (int) polymakeFile.size(); ++row)
	{

		for(int col = 0; col < (int)polynomialDegrees.size(); ++col)
		{
			cerr << "row, col=" << row << ", " << col << endl;
			ans[row].push_back(getStatisticsByFileDegree(polymakeFile[row][0], atoi(polynomialDegrees[col][0].c_str()), useDual));
		}//for
	}//for row. vertexCounts

	return ans;
}//getStatisticsByFile

/**
 * Given the (dual) polytope dim, vertex count, and polynomial degree,
 * Will find and return basic statistics (avg time, min/max, totals, etc) about this test class.
 */
ValuationDBStatistics IntegrationDB::getStatisticsByDimVertexDegree(int dim, int vertexCount, int degree, bool useDual)
{
	ValuationDBStatistics vdbs;
	vector<double> avgMinMaxCountLawrence, avgMinMaxCountTriangulate;

	//save how this function was called.
	vdbs.dim = dim;
	vdbs.vertexCount = vertexCount;
	vdbs.degree = degree;
	vdbs.useDual = useDual;

	//get avg, min, man, and number finished
	avgMinMaxCountLawrence    = getStatisticsAvgMinMaxCount(Lawrence, dim, vertexCount, degree, useDual);
	avgMinMaxCountTriangulate = getStatisticsAvgMinMaxCount(Triangulate, dim, vertexCount, degree, useDual);

	vdbs.avgTriangulationTime = avgMinMaxCountTriangulate[0];
	vdbs.avgLawrenceTime      = avgMinMaxCountLawrence[0];

	vdbs.minTriangulationTime = avgMinMaxCountTriangulate[1];
	vdbs.minLawrenceTime      = avgMinMaxCountLawrence[1];

	vdbs.maxTriangulationTime = avgMinMaxCountTriangulate[2];
	vdbs.maxLawrenceTime      = avgMinMaxCountLawrence[2];

	vdbs.totalFinishedTriangulationTestCases = avgMinMaxCountTriangulate[3];
	vdbs.totalFinishedLawrenceTestCases      = avgMinMaxCountLawrence[3];

	vdbs.totalTestCases = getNumberIntegrationTest(dim, vertexCount, degree, useDual);

	vdbs.manuallyLimitedLawrence = getLimit(Lawrence, dim, vertexCount, degree, useDual);
	vdbs.manuallyLimitedTriangulation = getLimit(Triangulate, dim, vertexCount, degree, useDual);

	return vdbs;
}//getStatisticsByDimVertexDegree


ValuationDBStatistics IntegrationDB::getStatisticsByFileDegree(const string & polymakeFile, int degree, bool useDual)
{
	ValuationDBStatistics vdbs;
	vector<double> avgMinMaxCountLawrence, avgMinMaxCountTriangulate;

	//save how this function was called.
	stringstream d, v;
	d << "select dim from polytope where polymakeFilePath = '" << polymakeFile << "'";
	v << "select vertexCount from polytope where polymakeFilePath = '" << polymakeFile << "'";
	vdbs.dim = queryAsInteger(d.str().c_str());
	vdbs.vertexCount = queryAsInteger(v.str().c_str());
	vdbs.degree = degree;
	vdbs.useDual = useDual;

	//get avg, min, man, and number finished
	avgMinMaxCountLawrence    = getStatisticsAvgMinMaxCount(Lawrence, polymakeFile, degree, useDual);
	avgMinMaxCountTriangulate = getStatisticsAvgMinMaxCount(Triangulate, polymakeFile, degree, useDual);

	vdbs.avgTriangulationTime = avgMinMaxCountTriangulate[0];
	vdbs.avgLawrenceTime      = avgMinMaxCountLawrence[0];

	vdbs.minTriangulationTime = avgMinMaxCountTriangulate[1];
	vdbs.minLawrenceTime      = avgMinMaxCountLawrence[1];

	vdbs.maxTriangulationTime = avgMinMaxCountTriangulate[2];
	vdbs.maxLawrenceTime      = avgMinMaxCountLawrence[2];

	vdbs.totalFinishedTriangulationTestCases = avgMinMaxCountTriangulate[3];
	vdbs.totalFinishedLawrenceTestCases      = avgMinMaxCountLawrence[3];

	vdbs.totalTestCases = getNumberIntegrationTest(polymakeFile, degree, useDual);

	vdbs.manuallyLimitedLawrence = getLimitByFile(Lawrence, polymakeFile, degree, useDual);
	vdbs.manuallyLimitedTriangulation = getLimitByFile(Triangulate, polymakeFile, degree, useDual);

	return vdbs;
}///getStatisticsByFileDegree


vector<double> IntegrationDB::getStatisticsAvgMinMaxCount(AlgorithemUsed alg, int dim, int vertexCount, int degree, bool useDual)
{
	vector<double> ans;
	vector<vector<string> > strAns;
	stringstream sql;
	string strAlg;

	strAlg = (alg == Lawrence ? "timeLawrence" : "timeTriangulate");

	if (useDual == true)
	{
		sql << "select avg(i." << strAlg << "), min(i." << strAlg << "), max(i." << strAlg << "), count(*) "
			<< " from integrate as i"
			<< " join polytope as dualP on dualP.rowid = i.polytopeID"
			<< " join polytope as orgP on orgP.rowid = dualP.dual"
			<< " join polynomial as p on p.rowid = i.polynomialID"
			<< " where dualP.dual is not null " //and with orgp
			<< " and orgP.dim = " << dim
			<< " and orgP.vertexCount = " << vertexCount
			<< " and p.degree = " << degree
			<< " and i." << strAlg << " >= 0 ";
	}//if dual
	else
	{

		sql << "select avg(i." << strAlg << "), min(i." << strAlg << "), max(i." << strAlg << "), count(*) "
			<< " from integrate as i"
			<< " join polynomial as p on p.rowid = i.polynomialID"
			<< " join polytope as t on t.rowid = i.polytopeID"
			<< " where t.dual is null"
			<< " and t.dim = " << dim
			<< " and t.vertexCount = " << vertexCount
			<< " and p.degree = " << degree
			<< " and i." << strAlg << " >= 0 ";
	}//regular

	//get the data and save it
	strAns = query(sql.str().c_str());

	ans.push_back(atof(strAns[0][0].c_str())); //avg
	ans.push_back(atof(strAns[0][1].c_str())); //min
	ans.push_back(atof(strAns[0][2].c_str())); //max
	ans.push_back(atof(strAns[0][3].c_str())); //count

	return ans;
}//	getStatisticsAvgMinMaxCount

vector<double> IntegrationDB::getStatisticsAvgMinMaxCount(AlgorithemUsed alg, const string &polymakeFile, int degree, bool useDual)
{
	vector<double> ans;
	vector<vector<string> > strAns;
	stringstream sql;
	string strAlg;

	strAlg = (alg == Lawrence ? "timeLawrence" : "timeTriangulate");

	if (useDual == true)
	{
		sql << "select avg(i." << strAlg << "), min(i." << strAlg << "), max(i." << strAlg << "), count(*) "
			<< " from polynomial as p, polytope as dualP, polytope as orgP, integrate as i "
			<< " where i.polynomialID = p.rowid and i.polytopeID = dualP.rowid " //join with p, i, dualP
			<< " and dualP.dual is not null and dualP.dual = orgP.rowid" //and with orgp
			<< " and orgP.polymakeFilePath = '" << polymakeFile << "'"
			<< " and p.degree = " << degree
			<< " and i." << strAlg << " >= 0 ";
	}//if dual
	else
	{

		sql << "select avg(i." << strAlg << "), min(i." << strAlg << "), max(i." << strAlg << "), count(*) "
			<< " from polynomial as p, polytope as t, integrate as i "
			<< " where i.polynomialID = p.rowid and i.polytopeID = t.rowid " //join with p, i, dualP
			<< " and t.dual is null"
			<< " and t.polymakeFilePath = '" << polymakeFile << "'"
			<< " and p.degree = " << degree
			<< " and i." << strAlg << " >= 0 ";
	}//regular

	//get the data and save it
	strAns = query(sql.str().c_str());

	ans.push_back(atof(strAns[0][0].c_str())); //avg
	ans.push_back(atof(strAns[0][1].c_str())); //min
	ans.push_back(atof(strAns[0][2].c_str())); //max
	ans.push_back(atof(strAns[0][3].c_str())); //count

	return ans;
}//getStatisticsAvgMinMaxCount


/**
 * Within all polynomials of set dim and degree, find my the ones that are NOT being used in a test with a polytope of set vertexCount and dual values.
 */
vector<vector<string> > IntegrationDB::getUnusedPolynomials(int dim, int degree, int vertexCount, bool useDual)
{
	stringstream sql;
	//GOAL: 1) Find all the polynomials that are of a set degree and dim.
	//      2) Find all the tests where the polynomials are of a set degree dim, the polytope is a set dim, and vertexCount (or came from a polytope of a set vertexCount if useDual=ture),
	//      3) Return the polynomials that are in the set 1 but NOT set 2.
	if ( useDual == false)
	{//for every polynomial of set dim and degree ans ask the question
		//Is this polynomial already being used in a test with a non-dual polytope of dim 10 and set vertexCount?
		//if not, then return it.
		sql << "select distinct rowid from polynomial "
				<< " where rowid not in"
				<< " ( select p.rowid from polynomial as p, integrate as i, polytope as t "
				<< "     where i.polytopeID = t.rowid and i.polynomialID = p.rowid " //join
				<< "       and t.dim = " << dim
				<< "       and t.vertexCount = " << vertexCount
				<< "       and p.degree = " << degree
				<< "       and t.dual is null "
				<< " ) "
				<< " and dim = " << dim << " and degree = " << degree;
	}
	else
	{//for every polynomial of set dim and degree ask the question:
		//Is this polynomial already being used in a test with a DUAL polytope of dim 10 where the org. polytope had set vertexCount?
		//if not, then return it.
		sql << "select distinct rowid from polynomial where rowid not in"
			<< "  (select p.rowid from polynomial as p, integrate as i, polytope as dualPolytope, polytope as orgPolytope "
			<< "     where i.polytopeID = dualPolytope.rowid and i.polynomialID = p.rowid "
			<< " 	   and dualPolytope.dual is not null "				//make sure we have dual
			<< "       and dualPolytope.dual = orgPolytope.rowid "		//find the org. polytope
			<< "       and orgPolytope.dim = " << dim					//look at the dim of the test case
			<< "       and orgPolytope.vertexCount = " << vertexCount	//look at the vertex count of the org polytope
			<< "       and p.degree = " << degree						//polynomial degree.
			<< "  ) "
			<< " and dim = " << dim << " and degree = " << degree;
	}
	return query(sql.str().c_str());
}//getUnusedPolynomials


vector<vector<string> > IntegrationDB::getUnusedPolynomials(int dim, int degree, int polytopeID)
{
	stringstream sql;
	sql << "select distinct rowid from polynomial where rowid not in "
			<< " ( select p.rowid from polynomial as p, integrate as i, polytope as t "
			<< "     where i.polytopeID = t.rowid and i.polynomialID = p.rowid " //join
			<< "       and t.rowid = " << polytopeID
			<< "       and p.degree = " << degree
			<< " ) "
			<< " and dim = " << dim << " and degree = " << degree;
	return query(sql.str().c_str());
}//getUnusedPolynomials


/**
 * Given a set dim, degree, vertexCount and dual values, find me all polytopes that are not used in the integrate table.
 */
vector<vector<string> > IntegrationDB::getUnusedPolytopes(int dim, int degree, int vertexCount, bool useDual)
{
	stringstream sql;
	//GOAL: 1) Find all the polytopes that are of a set dim and vertex count (or game from a polytope of a set vertex count if useDual = true).
	//      2) Find all the tests where the polynomials are of a set degree dim, the polytope is a set dim, and vertexCount (or came from a polytope of a set vertexCount if useDual=ture),
	//      3) Return the polytopes that are in the set 1 but NOT set 2.
	if (useDual == false)
	{//for every non-dual polytope of set dim and vertexCount ask the question:
		//Is this polytope already being used in a test with a polynomial of set degree and vertexCount?
		//if not, then return it
		sql << "select distinct rowid from polytope where rowid not in"
			<< " ( select t.rowid from polynomial as p, integrate as i, polytope as t "
			<< "     where i.polytopeID = t.rowid and i.polynomialID = p.rowid " //join
			<< "       and t.dim = " << dim
			<< "       and t.vertexCount = " << vertexCount
			<< "       and p.degree = " << degree
			<< "       and t.dual is  null "
			<< " ) "
			<< " and dim = " << dim << " and vertexCount = " << vertexCount
			<< " and dual is  null";
	}
	else
	{//for every dual polytope of set dim and that came from an org. polytope of set vertexCount ask the question:
		//Is this polytope already being used in a test with a polynomial of set degree and vertexCount?
		//if not, then return it
		sql << "select distinct dualP.rowid from polytope as dualP, polytope as orgP "
			<< " where dualP.rowid not in "
			<< " ( select dualPoly.rowid from polynomial as p, integrate as i, polytope as dualPoly, polytope as orgPoly "
			<< "     where i.polytopeID = dualPoly.rowid and i.polynomialID = p.rowid " //join
			<< "       and dualPoly.dual is not null"
			<< "       and dualPoly.dual = orgPoly.rowid"
			<< "       and orgPoly.dim = " << dim
			<< "       and orgPoly.vertexCount = " << vertexCount
			<< "       and p.degree = " << degree
			<< " ) "
			<< " and dualP.dim = " << dim
			<< " and dualP.dual is not null"
			<< " and dualP.dual = orgP.rowid" //join the orgPoly and dualPoly.
			<< " and orgP.vertexCount = " << vertexCount;
	}
	return query(sql.str().c_str());
}//getUnusedPolytopes


/**
 * No longer used....to delte.
 *
 *	Inserts a polynomial, a polytope, and the dual polytope in the polynomial/polytope tables.
 * Also, inserts two integration tests (one for the org. polytope and one for the dual polytope).
 *
 * @parm polynomialPath:	file name to the latte-style polynomial file
 * @parm dim: dim of both polytopes and the polynomial.
 * @parm degree: of the polynomial
 * @parm (dual)polytopePath: file name of the (dual) polytope
 * @parm (dual)polymakePath: file name of the (dual) polymake file.
 * @parm (dual)vertexCount: number of vertices in the (dual) polytope
 * @parm (dual)simple:	is the (dual) polytope simple? Yes=true=1.
 */
void IntegrationDB::insertEmptyIntegrationTest(
			const char* polynomialPath, int dim, int degree,
			const char* polytopePath, const char* polymakePath, int vertexCount, int simple,
			const char* dualPolytopePath, const char* dualPolymakePath, int dualVertexCount, int dualSimple)
{
	int polynomialID, polytopeID, dualPolytopeID;
	
	//insert the dual polytope first (and get the rowid) and then insert the org. polytope.
	dualPolytopeID = insertPolytope(dim, dualVertexCount, dualSimple, -1            , dualPolytopePath, dualPolymakePath);
	polytopeID     = insertPolytope(dim, vertexCount    , simple    , dualPolytopeID, polytopePath, polymakePath);

	//now insert the polynomial.
	polynomialID   = insertPolynomial(dim, degree, polynomialPath);
	
	//now insert two test cases: org. and dual polytope.
	insertIntegrationTest(polynomialID, polytopeID);
	insertIntegrationTest(polynomialID, dualPolytopeID);
}//insertEmptyIntegrationTest

/**
 *	Inserts 1 empty row in the integrate table.
 *
 * @parm polynomialID: rowid of a polynomial from the polynomial table.
 * @parm polytopeID: rowid of a polytope from the polytope table.
 * @return rowid of the new just-inserted integrate row.
 */
int IntegrationDB::insertIntegrationTest(int polynomialID, int polytopeID)
{
	stringstream sql;
	sql << "insert into integrate (polynomialID, polytopeID, timeLawrence, timeTriangulate, integral) " 
	    << "values (" << polynomialID << ", " << polytopeID << ", " << "-1, -1, 'NA')";
	query(sql.str().c_str());
	return last_insert_rowid();
}

/**
 * Pre conditions: the polynomial and polytope tables have enough rows of a
 * set dim/degree/vertex count to make 'count' many integration tests
 * and dual integration tests.
 *
 * If there is enough, we will then pick unused polytopes and polynomials and
 * insert there combination into the integrate table.
 */
void IntegrationDB::insertIntegrationTest(int dim, int degree, int vertexCount, int count)
{
	//find how many of each we have.
	int numPolynomials          = getNumberPolynomials(dim, degree);
	//cout << "numPolynomials" << numPolynomials << endl;
	int numPolytopes            = getNumberPolytopes(dim, vertexCount, false);
	//cout << "numPolytope" << numPolytopes <<endl;
	int numDualPolytopes        = getNumberPolytopes(dim, vertexCount, true);
	//cout << "numDualPolytope" << numDualPolytopes <<endl;
	int numIntegrationTests     = getNumberIntegrationTest(dim, vertexCount, degree, false);
	//cout << "numIntegrationTests" << numIntegrationTests << endl;
	int numDualIntegrationTests = getNumberIntegrationTest(dim, vertexCount, degree, true);
	//cout << "numDualIntegrationTests" << numDualIntegrationTests << endl;


	//check if we really can make 'count' many integration tests and dual integration tets.
	if (numPolynomials < count)
		throw SqliteDBexception("insertIntegrationTest::Not enough polynomials exist");
	if (numPolynomials < count)
		throw SqliteDBexception("insertIntegrationTest::Not enough polytopes exist");
	if (numDualPolytopes < count)
		throw SqliteDBexception("insertIntegrationTest::not enough dual polytopes exist");

	//ok, now we know we can make count many integration tests! so lets do it.

	if ( count - numIntegrationTests  <= 0)
	{
		cout << "There already exist " << numIntegrationTests << " tests; in fact, there might be more." << endl;
		return;
	}

	if (numIntegrationTests < count)
		makeMoreIntegrationTests(dim, degree, vertexCount, false, count, numIntegrationTests);
	if (numDualIntegrationTests < count)
		makeMoreIntegrationTests(dim, degree, vertexCount, true, count, numIntegrationTests);
}//insertIntegrationTest

void  IntegrationDB::insertSpecficPolytopeIntegrationTest(string polymakeFile, int degree, int count)
{
	stringstream sql;
	sql << "select dim from polytope where polymakeFilePath = '" << polymakeFile << "'";

	int dim = queryAsInteger(sql.str().c_str());
	int numPolynomials = getNumberPolynomials(dim, degree);
	int rowid = doesPolytopeExist(polymakeFile.c_str());
	int numIntegrationTests = getNumberIntegrationTest(rowid, degree);

	if ( rowid <= 0)
		throw SqliteDBexception("insertSpecficPolytopeIntegrationTest::polytope does not exist");
	if (count > numPolynomials)
		throw SqliteDBexception("insertSpecficPolytopeIntegrationTest::Not enough polynomials exist");
	if ( count - numIntegrationTests <= 0)
	{
		cout << "polytope dim: " << dim
			 << "\n number of polynomials: " << numPolynomials
			 << "\n rowid of polytope (" << polymakeFile.c_str() << "): " << rowid
			 << "\n number of current integration tests: " << numIntegrationTests << endl;
		cout << "There already exist " << numIntegrationTests << " tests." << endl;

		return;
	}

	makeMoreIntegrationTests(rowid, dim, degree, count, numIntegrationTests);

}//insertSpecficPolytopeIntegrationTest



/**
 * Inserts 1 row in the polynomial table.
 *
 * @parm dim: dim. of the polynomial (number of variables).
 * @parm degree: of the polynomial
 * @parm filePath: to the latte-style polynomial.
 * @return rowid of the inserted row.
 */
int IntegrationDB::insertPolynomial(int dim, int degree, const char*filePath) throw(SqliteDBexception)
{
	if ( doesPolynomialExist(filePath))
		throw SqliteDBexception(string("insertPolynomial::Polynomial ")+filePath+" already exist");

	stringstream sql;
	sql << "insert into polynomial (dim, degree, filepath) values (" << dim << ", " << degree << ", '" << filePath << "')";
	query(sql.str().c_str());	
	return last_insert_rowid();
}//insertPolynomial


/**
 * Inserts 1 row in the polytope table.
 *
 * @parm dim: dim. of polytope = dime of polynomial.
 * @parm vertexCount: what do you think?
 * @parm simple: true if the polytope is simple
 * @parm dualRowID: if positive, it is a "pointer" to the polytope table for the dual polytope. The dual does not point back to its "parent"
 * @parm latteFilePath: file name of the latte file.
 * @parm polymakeFilePath: polymake file path. could be null.
 */
int IntegrationDB::insertPolytope(int dim, int vertexCount, int simple, int dualRowID, const char* latteFilePath, const char* polymakeFilePath)
{
	stringstream sql;
	sql << "insert into polytope (dim, vertexCount, simple, latteFilePath, polymakeFilePath, dual) values (" 
	    << dim << ", " << vertexCount << ", " << simple  << ", '" << latteFilePath << "', '" << polymakeFilePath
	    << "', ";
	if (dualRowID > 0 )
		sql << dualRowID;
	else
		sql << "NULL";
	sql <<")";
	query(sql.str().c_str());
	return last_insert_rowid();
}//insertPolytope

//no longer used.
//to delete one day....
int IntegrationDB::insertPolytopeAndPickIntegrationTest(int dim, int vertexCount,     int simple    , const char * latteFile    , const char * polymakeFile
												               , int dualVertexCount, int dualSimple, const char * dualLatteFile, const char * dualPolymakeFile)
{
	stringstream sql;
	sql << "select * from polytope where polymakeFilePath = '" << latteFile << "' or polymakeFilePath = '" << dualLatteFile << "'";
	if ( query(sql.str().c_str()).size() )
	{
		throw SqliteDBexception("insertPolytopeAndPickIntegrationTest::Database already contains those polymake files!!!");
	}//if error.

	sql.str("");//clear query string.

	//now find 1 polynomial that have not been used on a polytope of this dim and vertex count.
	sql << "select rowid from polynomial where rowid not in "
			<<"(select p.rowid from polynomial as p, integrate as i, polytope as t "
			<<         "where i.polytopeID = t.rowid and i.polynomialID = p.rowid " // join the tables
			<<         "   and t.dim = " << dim << " and t.vertexCount = " << vertexCount  //restrict to the current dim/vertexcount
			<<") limit 1";
	//get this unused rowid.
	vector<vector<string> > result;
	result = query(sql.str().c_str());
	if ( ! result.size())
		throw SqliteDBexception("insertPolytopeAndPickIntegrationTest::There are not enough unique unused polynomials");
	int polynomialID = atoi(result[0][0].c_str());
	assert(polynomialID > 0);

	//insert two polytopes.
	int polytopeID     = insertPolytope(dim, vertexCount    , simple    , -1        , latteFile    , polymakeFile);
	int dualPolytopeID = insertPolytope(dim, dualVertexCount, dualSimple, polytopeID, dualLatteFile, dualPolymakeFile);

	//insert two tests with the same polynomial.
	insertIntegrationTest(polynomialID, polytopeID);
	insertIntegrationTest(polynomialID, dualPolytopeID);
}//insertPolytopeAndPickIntegrationTest


/**
 * Returns the number of completed test in the integrate table for the specfic test class
 */
int IntegrationDB::testCasesCompleted(AlgorithemUsed alg, int dim, int vertexCount, int degree, bool useDual)
{
	stringstream sql;
	if ( useDual == true)
	{
		sql << "select count(*) from polynomial as p, integrate as i, polytope as dualP, polytope as orgP "
			<< " where i.polynomialID = p.rowid and i.polytopeID = dualP.rowid "
			<< " and dualP.dual is not null and dualP.dual = orgP.rowid "
			<< " and orgP.vertexCount = " << vertexCount
			<< " and orgP.dim = " << dim
			<< " and p.degree = " << degree;
	}
	else
	{
		sql << "select count(*) from polynomial as p, integrate as i, polytope as t "
			<< " where i.polynomialID = p.rowid and i.polytopeID = t.rowid "
			<< " and t.vertexCount = " << vertexCount
			<< " and t.dim = " << dim
			<< " and p.degree = " << degree
			<< " and t.dual is null ";
	}//else not dual polytopes.

	sql << " and " << (alg == Lawrence  ? " i.timeLawrence " : " i.timeTriangulate ") << " >= 0";
	return queryAsInteger(sql.str().c_str());
}//isTestCaseFinished

/**
 *
 * Private function.
 *
 * Get a list of unused polynomials and polytopes, and add them to the integrate table.
 */
void IntegrationDB::makeMoreIntegrationTests(int dim, int degree, int vertexCount, bool useDual, int requestedCount, int existingCount)
{
	int newRows = requestedCount - existingCount;
	//unusedPolynomials and unusedPolytopes has 1 column of rowid's that do not already exist in the integrate table
	//again, unusedPolynomials[i] is a vector with 1 element.
	vector<vector<string> > unusedPolynomials = getUnusedPolynomials(dim, degree, vertexCount, useDual);
	vector<vector<string> > unusedPolytopes   = getUnusedPolytopes(dim, degree, vertexCount, useDual);

	if ( unusedPolynomials.size() < newRows || unusedPolytopes.size() < newRows)
		throw SqliteDBexception("makeMoreIntegrationTests: there are not enough free polynomials or polytopes"); //I think this should never be true...

	for(int i = 0; i < newRows; ++i)
	{
		insertIntegrationTest(atoi(unusedPolynomials[i][0].c_str()), atoi(unusedPolytopes[i][0].c_str()));
	}//for i
}//makeMoreIntegrationTests

/**
 * Private function.
 * get a list of unused polynomials for this set polytope and add them to the integrate table.
 */
void IntegrationDB::makeMoreIntegrationTests(int polytopeID, int dim, int degree, int requestedCount, int existingCount)
{
	int newRows = requestedCount - existingCount;
	vector<vector<string> > unusedPolynomials = getUnusedPolynomials(dim, degree, polytopeID);

	if ( unusedPolynomials.size() < newRows)
		throw SqliteDBexception("makeMoreIntegrationTests: there are not enough free polynomials"); //I think this should never be true...

	for(int i = 0; i < newRows; ++i)
	{
		insertIntegrationTest(atoi(unusedPolynomials[i][0].c_str()), polytopeID);
	}//for i
}//makeMoreIntegrationTests

/**
 * Updates the integral table with time and valuation results.
 *
 * @alg: what column should be updated in the integrate table?
 * @time: time from calling the mainValuationDriver
 * @currentValue: new value for the integral column. Checks to see if this is different from any previous values.
 * @previousValue: current value in the integral column, or the string "NA"
 * @rowid: which row are we updating?
 */
void IntegrationDB::updateIntegrationTimeAndValue(AlgorithemUsed alg, double time, RationalNTL computedValue, string previousValueStr,  string rowid)
{
	stringstream sql;

	sql << "update integrate set ";
	if ( alg == Lawrence)
		sql << " timeLawrence = " << time;
	else //Triangulate};
		sql << " timeTriangulate = " << time;


	if ( previousValueStr == "NA")
	{
		sql << " , integral = ' " << computedValue << "' ";
	}
	else
	{
		RationalNTL previousValue(previousValueStr);
		if ( previousValue != computedValue)
		{
			throw SqliteDBexception(string("updateIntegrationTimeAndValue::The integrals differ")
					+"\n\tpreviousValue: " + previousValueStr
					+"\n\tcomputedValue: " + computedValue.str()
					+"\n\tcurrent sql stm:" + sql.str()
					+"\n\trowid: " + rowid);
		}
	}
	sql << " , flagValue = '" << -1 << "', flagType = '" << LAWRECE_INTEGRATE_VERSION << "'";
	sql << " where rowid = " << rowid << endl;
	query(sql.str().c_str());

}//updateIntegrationTimeAndValue
