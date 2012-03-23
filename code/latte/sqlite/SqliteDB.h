#ifndef SQLITEDB_H
#define SQLITEDB_H

/** class SqliteDB is a basic class for using a sqlite database.
 *  In applications, you would want to subclass SqliteDB and write new query functions specific to the application.
 *  This class makes no assumptions about the form of the database. 
 *
 *  If you are new to SQL, do not worry. It is easy, there are only about 5 commands you really need to know.
 *  For a good quck introduction look at http://www.w3schools.com/sql/default.asp
 *
 *  There are 3 differences you need to know. SQL is a language (just like c++), and SQLite is
 *  just one (almost incomplete) implementation of the language and so has its own differences
 *  (just like how g++ and other compilers differ). Finally, we need an API to run sql.
 *  So there are 3 things to learn: basic sql, sql that sqlite uses, and the SQLite API.
 *
 *  The sqlite lib is passed to g++ with "-lsqlite"
 *
 *  This is my first time using sqlite in c++, so if you know of a better way to do things, go for it and let me know.
 *  ~Brandon
 */

#include <string>
#include <vector>
#include <sqlite3.h>
#include <exception>

using namespace std;

class SqliteDBexception: public exception
{
private:
	string msg;
public:
	SqliteDBexception () throw();
	SqliteDBexception (const char*) throw();
	SqliteDBexception (const string) throw();
    virtual ~SqliteDBexception() throw();
	virtual const char* what() const throw();
};

class SqliteDB
{
public: 
	//A B C D E F G H I J K L M N O P Q R S T U V W X Y Z
	SqliteDB(); 
	~SqliteDB();
	
	void close(); 


	int last_insert_rowid(void);	     //row id of the last inserted row that used this database pointer (dbPtr).
	void open(const char *filename);	//opens the database connection.
	void open(const string &filename);

	void printResults(const vector<vector<string> > &r); //this is used for debugging.
	
	vector<vector<string> > query(const char* query); //basic query function. returns results as 2d array of strings.
	vector<vector<string> > query(const string &query);

	double queryAsFloat(const char* q);
	vector<double> queryAsFloatArray(const char *q);
	int queryAsInteger(const char* query); //for querys that should only return 1 integer (ex: select count(*) from...)
 		
private: 

	sqlite3 *dbPtr; //database pointer 


};



#endif 



