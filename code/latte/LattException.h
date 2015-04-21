/*
 * LattException.h
 *
 *  Created on: Jul 19, 2011
 *      Author: bedutra
 */

#ifndef LATTEXCEPTION_H_
#define LATTEXCEPTION_H_

#include <exception>
#include <string>
#include <iostream>
#include <sstream>

using namespace std;

/**
 * Interface Examples:
 *   All the THROW_LATTE* mactos are exceptions. They differ in what message .what() returns.
 *
 *   THROW_LATTE(LattException::bug_Unknown); --In catch statement, if printing .what(), will print standard error message (error code, line and file numbers)
 *   THROW_LATTE(LattException::bug_Unknown, bool); --If bool is ...
 *                                                  --  1, then .what() will return the standard error message
 *                                                  --  0, and LATTEXCEPTION_PRINTSTATUS=0, then .what() will return an empty string
 *                                                  --  0, and LATTEXCEPTION_PRINTSTATUS=1, then .what() will return the standard error message. Good for debugging.
 *   THROW_LATTE(LattException::bug_Unknown, bool, msg); The value of bool plays the same job as above. msg is a string that gives more detail on the error.
 */
#define LATTEXCEPTION_PRINTSTATUS 0


//Neat trick to override a macro function. Can use THROW_LATTE(code, bool) or THROW_LATTE(code)
#define LATTE_GET_MACRO(_1,_2, NAME,...) NAME

//Should ONLY use THROW_LATTE(...) and NEVER use THROW_LATTE1 and THROW_LATTE2
#define THROW_LATTE(...) LATTE_GET_MACRO(__VA_ARGS__, THROW_LATTE2, THROW_LATTE1)(__VA_ARGS__)
#define THROW_LATTE1( type ) throw LattException( (type), __FILE__, __LINE__, 1, "" )
#define THROW_LATTE2( type , printBool) throw LattException( (type), __FILE__, __LINE__, (printBool), "" )

#define THROW_LATTE_MSG( type, printBool, msg) throw LattException( (type), __FILE__, __LINE__, (printBool), (msg) )


class LattException : public exception
{
public:
	enum UserError
	{
		none,  						//no error

		//user errors (ue)
		ue_UnknownCommandLineOption,   //bad option
		ue_BadCommandLineOption,		//bad combination of options.
		ue_BadFileOption,				//bad file keyword or style
		ue_BadCommandLineOptionCount,	//unexpected command count
		ue_HelpMenuDisplayed,			//not really an error.
		ue_FileNameMissing,				//missing file name
		ue_BadPolynomialLinFormInput,   //incorrect polynomial or linear form input

		//polyhedra errors (pe)
		pe_RationalPolytope, 			//expecting integer-vertex polytope.
		pe_Unbounded,					//expecting bounded polyheda
		pe_UnexpectedRepresentation,

		//file  error (fe)
		fe_Open,						//cannot open file or does not exist.
		fe_Parse,						//parse error

		//integration errors (ie)
		ie_BadIntegrandFormat,			//something's wrong with the integrand.
		ie_UnexpectedIntegrationOption,	//maybe wrong algorithm used

		//division errors
		de_divisionByZero, 				//divided by zero...perturbation is not

		//unknown (bug)
		bug_Unknown,					//something really bad happened.
		bug_NotImplementedHere,			//this case does not exist is the function that throw it.
	};



	virtual ~LattException() throw();
	LattException(UserError ue, const char * file, const int line, const bool printException, const char * message = "");

	virtual const char* what() const throw();

private:
	string printErrorMessages() const;

	UserError userError;
	string msg; //error message from programmer.
	int lineNumber;
	string fileName;
	bool printStatus;
};

#endif /* LATTEXCEPTION_H_ */
