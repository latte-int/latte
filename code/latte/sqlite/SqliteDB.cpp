

#include "SqliteDB.h"
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <unistd.h> //for sleep()


//A B C D E F G H I J K L M N O P Q R S T U V W X Y Z
SqliteDBexception::SqliteDBexception() throw() {}
SqliteDBexception::SqliteDBexception(const char* s) throw() {msg = s;}
SqliteDBexception::SqliteDBexception (const string s) throw() {msg = s;}
SqliteDBexception::~SqliteDBexception() throw() {}
const char* SqliteDBexception::what() const  throw() {return msg.c_str();}



SqliteDB::SqliteDB(): dbPtr(NULL)
{}


SqliteDB::~SqliteDB()
{
	close();
}

void SqliteDB::close()
{
	sqlite3_close(dbPtr); 
}


int SqliteDB::last_insert_rowid()
{
	return sqlite3_last_insert_rowid(dbPtr);
}

/**
 *	@throw SqliteDBexception if cannot open db.
 */
void SqliteDB::open(const char *filename)
{
	if( sqlite3_open(filename, &dbPtr))
	{
		throw SqliteDBexception("open::Cannot open database file");
	}//true = could not open db. false = no errors
}

void SqliteDB::open(const string &filename)
{
	open(filename.c_str());
}

 //this is used for debugging.
 //print the resutls from query()
void SqliteDB::printResults(const vector<vector<string> > &v)
{

	for(vector<vector<string> >::const_iterator i = v.begin(); i != v.end(); ++i)
	{
		for(vector<string>::const_iterator j = (*i).begin(); j != (*i).end(); ++j)
			cout <<  *j << "\t";
		cout << endl;
			
	}
}

/**
* Returns the results of an SQL query as an array of strings.
* Null values are returned as the string "NULL"
* @parm query: sql string to run
* @return 2d array of results w/o any table heading.
*/
vector<vector<string> > SqliteDB::query(const char* query)
{
	sqlite3_stmt *statement;
	vector<vector<string> > results;
	int queryErrorCode;

	queryErrorCode =sqlite3_prepare_v2(dbPtr, query, -1, &statement, 0);
	while ( queryErrorCode == SQLITE_BUSY)
	{
		cout << "Database is busy, sleeping for 5sec." << endl;
		sleep(5);
		cout << "Re-running sql:" << query << endl;
		queryErrorCode = sqlite3_prepare_v2(dbPtr, query, -1, &statement, 0);
	}//while db is busy.


	if(queryErrorCode == SQLITE_OK)
	{
		int cols = sqlite3_column_count(statement);
		
		//The commented-out line can be used to get the header row (the names of the columns)
		//const char *sqlite3_column_name(sqlite3_stmt*, int N); gets header infor for nth col(starting from zero???).
		int result = 0;
		while(true)
		{
			result = sqlite3_step(statement);
	             
			if(result == SQLITE_ROW)
			{
				vector<string> values;

				for(int col = 0; col < cols; col++)
				{
					char *resultsPtr = (char*)sqlite3_column_text(statement, col);
					if ( resultsPtr)
						values.push_back(resultsPtr);
					else
						values.push_back("NULL");
				}//for every col in the current row.
				results.push_back(values);
			}//if there are results
			else
			{
				break;  
			}
		}//while there are more rows
	        
		sqlite3_finalize(statement);//free up the memory.
	}//if no errors
	else
	{
		string error = sqlite3_errmsg(dbPtr);
		stringstream s;
		s << "SqliteDB ERROR:: ." << query << ". " << error << endl;
		cout << s.str().c_str() << endl;
		throw SqliteDBexception(s.str());
	}//else there was errors running the sql statment
	     
	return results; 
}//query


vector<vector<string> > SqliteDB::query(const string &q)
{
	return query(q.c_str());
}


/**
 * For queries that return 1 int (example: select count(*) form...where...)
 */
int SqliteDB::queryAsInteger(const char* q)
{
	vector<vector<string> > result = query(q);
	if ( result.size() != 1 && result[0].size() !=1 )
		throw SqliteDBexception("queryAsInteger::Query contains more than 1 result");

	return atoi(result[0][0].c_str());
}

/**
 * For queries that return 1 floating-point number (example: select avg(*) form...where...)
 */
double SqliteDB::queryAsFloat(const char* q)
{
	vector<vector<string> > result = query(q);
	if ( result.size() != 1 && result[0].size() !=1 )
		throw SqliteDBexception("queryAsFloat::Query contains more than 1 result");

	return atof(result[0][0].c_str());
}



/**
 * For queries that return 1 column of numerical data.
 */
vector<double> SqliteDB::queryAsFloatArray(const char *q)
{
	vector<vector<string> > result = query(q);
	vector<double> ans;

	ans.resize(result.size());

	for(int i = 0; i < result.size(); ++i)
	{
		if ( result[i].size() != 1)
			throw SqliteDBexception("queryAsFloatArray::Query contains more than 1 column");

		if( result[i][0].compare("NULL") == 0)
			ans[i] = 0.0;
		else
			ans[i] = atof(result[i][0].c_str());
	}//for i.

	return ans;
}//queryAsFloatArray



