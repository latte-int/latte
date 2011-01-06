
#include "IntegrationDB.h"
#include <sstream>
#include <cstdlib> //for atoi
#include <cassert>
#include <iostream>
IntegrationDB::IntegrationDB()
{
}

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
 * Returns the number of polynomials in the polynomials table with this dim and degree
 */
int IntegrationDB::getNumberPolynomials(int dim, int degree)
{
	stringstream sql;
	sql << "select count(*) from polynomial where dim = " << dim << " and degree = " << degree;
	return queryAsInteger(sql.str().c_str());
}//getNumberPolynomials


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


//select rowid from polynomial where rowid not in
//	(select p.rowid from polynomial as p, integrate as i, polytope as t
	//where i.polytopeID = t.rowid and i.polynomialID = p.rowid
	 //  and t.dim = 5 and p.degree = 32 and t.vertexCount = 10
	//) and dim = 5 and degree = 32
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
			const char* polytopePath, const char* polymakePath, int vertexCount, bool simple,
			const char* dualPolytopePath, const char* dualPolymakePath, int dualVertexCount, bool dualSimple)
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

void IntegrationDB::insertIntegrationTest(int dim, int degree, int vertexCount, int count)
{
	int numPolynomials          = getNumberPolynomials(dim, degree);
	cout << "numPolynomials" << numPolynomials << endl;
	int numPolytopes            = getNumberPolytopes(dim, vertexCount, false);
	cout << "numPolytope" << numPolytopes <<endl;
	int numDualPolytopes        = getNumberPolytopes(dim, vertexCount, true);
	cout << "numDualPolT" << numDualPolytopes <<endl;
	int numIntegrationTests     = getNumberIntegrationTest(dim, vertexCount, degree, false);
	cout << "numIntegrationTests" << numIntegrationTests << endl;
	int numDualIntegrationTests = getNumberIntegrationTest(dim, vertexCount, degree, true);
	cout << "numDualIntegrationTests" << numDualIntegrationTests << endl;


	if (numPolynomials < count)
		throw SqliteDBexception("insertIntegrationTest::Not enough polynomials exist");
	if (numPolynomials < count)
		throw SqliteDBexception("insertIntegrationTest::Not enough polytopes exist");
	if (numDualPolytopes < count)
		throw SqliteDBexception("insertIntegrationTest::not enough dual polytopes exist");
	//ok, now we know we can make count many integration tests! so lets do it...

	if ( count - numIntegrationTests  <= 0)
	{
		cout << "There already exist " << numIntegrationTests << " tests, but added any more" << endl;
		return;
	}

	if (numIntegrationTests < count)
		makeMoreIntegrationTests(dim, degree, vertexCount, false, count, numIntegrationTests);
	if (numDualIntegrationTests < count)
		makeMoreIntegrationTests(dim, degree, vertexCount, true, count, numIntegrationTests);
}//insertIntegrationTest

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
int IntegrationDB::insertPolytope(int dim, int vertexCount, bool simple, int dualRowID, const char* latteFilePath, const char* polymakeFilePath)
{
	stringstream sql;
	sql << "insert into polytope (dim, vertexCount, simple, latteFilePath, polymakeFilePath, dual) values (" 
	    << dim << ", " << vertexCount << ", " << (simple ? 1 : 0) << ", '" << latteFilePath << "', '" << polymakeFilePath 
	    << "', ";
	if (dualRowID > 0 )
		sql << dualRowID;
	else
		sql << "NULL";
	sql <<")";
	query(sql.str().c_str());
	return last_insert_rowid();
}//insertPolytope


int IntegrationDB::insertPolytopeAndPickIntegrationTest(int dim, int vertexCount,     bool simple    , const char * latteFile    , const char * polymakeFile
												, int dualVertexCount, bool dualSimple, const char * dualLatteFile, const char * dualPolymakeFile)
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

										

void IntegrationDB::makeMoreIntegrationTests(int dim, int degree, int vertexCount, bool useDual, int requestedCount, int existingCount)
{
	int newRows = requestedCount - existingCount;
	//unusedPoly* has 1 column of rowid's that do not already exist in the integrate table
	vector<vector<string> > unusedPolynomials = getUnusedPolynomials(dim, degree, vertexCount, useDual);
	vector<vector<string> > unusedPolytopes   = getUnusedPolytopes(dim, degree, vertexCount, useDual);

	if ( unusedPolynomials.size() < newRows || unusedPolytopes.size() < newRows)
		throw SqliteDBexception("makeMoreIntegrationTests: there are not enough free polynomials or polytopes"); //I think this should never be true...

	for(int i = 0; i < newRows; ++i)
	{
		insertIntegrationTest(atoi(unusedPolynomials[i][0].c_str()), atoi(unusedPolytopes[i][0].c_str()));
	}//for i
}//makeMoreIntegrationTests

