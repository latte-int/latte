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

/// FIXME: Don't include the semicolon. Such macros should always be
/// usable in the same way as ordinary statements. --mkoeppe
#define THROW_LATTE( type ) throw LattException( (type), __FILE__, __LINE__, "" );
#define THROW_LATTE_MSG( type, msg ) throw LattException( (type), __FILE__, __LINE__, (msg) );

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

		//polyhedra errors
		pe_RationalPolytope, 			//expecting integer-vertex polytope.
		pe_Unbounded,					//expecting bounded polyheda
		pe_UnexpectedRepresentation,

		//file  error.
		fe_Open,						//cannot open file or does not exist.
		fe_Parse,						//parse error

		//integration errors.
		ie_badIntegrandFormat,			//something's wrong with the integrand.

		//unknown
		bug_Unknown,					//somthing really bad happened.
		bug_NotImplementedHere,			//this case does not exist is the function that throw it.
	};



	virtual ~LattException() throw();
	LattException(UserError ue, const char * file, const int line, const char * message = "");

	virtual const char* what() const throw();

private:
	void printErrorMessages() const;

	UserError userError;

	string msg; //error message from programmer.
	int lineNumber;
	string fileName;
};

#endif /* LATTEXCEPTION_H_ */
