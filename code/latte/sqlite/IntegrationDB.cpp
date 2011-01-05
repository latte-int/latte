
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
	sql << "select count(*) from polytope where dim = " << dim
		<< " and vertexCount = " << vertexCount
		<< " and dual is " << (useDuals ? "not":"") << " null";
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
	sql << "select count(*) from polynomial as p, polytope as t, integrate as i "
		<< "where i.polynomialID = p.rowid and i.polytopeID = t.rowid " //join
		<< " and t.dim = " << dim << " and t.vertexCount = " << vertexCount
		<< " and p.degree = " << degree
		<< " and t.dual is " << (useDuals ? "not":"") << " null ";
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
	sql << "select distinct rowid from polynomial where rowid not in"
		<< " ( select p.rowid from polynomial as p, integrate as i, polytope as t "
		<< "     where i.polytopeID = t.rowid and i.polynomialID = p.rowid " //join
	    << "       and t.dim = " << dim << "and t.vertexCount = " << vertexCount
	    << "       and p.degree = " << degree
	    << "       and t.dual is " << (useDual ? "not":"") << " null "
	    << " ) "
		<< " and dim = " << dim << " and degree = " << degree;
	return query(sql.str().c_str());
}//getUnusedPolynomials

vector<vector<string> > IntegrationDB::getUnusedPolytopes(int dim, int degree, int vertexCount, bool useDual)
{
	stringstream sql;
	sql << "select distinct rowid from polytope where rowid not in"
		<< " ( select t.rowid from polynomial as p, integrate as i, polytope as t "
		<< "     where i.polytopeID = t.rowid and i.polynomialID = p.rowid " //join
	    << "       and t.dim = " << dim << "and t.vertexCount = " << vertexCount
	    << "       and p.degree = " << degree
	    << "       and t.dual is " << (useDual ? "not":"") << " null "
	    << " ) "
		<< " and dim = " << dim << " and vertexCount = " << vertexCount
		<< " and dual is " << (useDual ? "not":"")  << " null";
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

