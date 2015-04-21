/*
 * LattException.cpp
 *
 *  Created on: Jul 19, 2011
 *      Author: bedutra
 */

#include "LattException.h"


LattException::~LattException() throw() {}

LattException::LattException(UserError ue, const char * file, const int line, const bool b, const char * message)
{
	userError = ue;
	fileName = file;
	lineNumber = line;
	msg = message;
	printStatus = b;
}



string  LattException::printErrorMessages() const
{
	stringstream out;
	out << "\nLatte Exception"
		 << "\n  Error code : " << userError << ". ";
		
	switch( userError)
	{
		case none: 								out << ""; break;
		case ue_UnknownCommandLineOption: 
		case ue_BadCommandLineOption: 			out << "(Unknown command line options)"; break;
		case ue_BadFileOption: 					out << "(Wrong file keyword or style)"; break;
		case ue_BadCommandLineOptionCount: 		out << "(Unexpected command count)"; break;
		case ue_HelpMenuDisplayed: 				out << "(Help menu displayed)";break;
		case ue_FileNameMissing: 				out << "(Missing file name)";break;
		case ue_BadPolynomialLinFormInput: 		out << "(Incorrect polynomial or linear form input)";break;
		case pe_RationalPolytope: 				out << "(Expecting integer-vertex polytope)";break;
		case pe_Unbounded: 						out << "(Expecting bounded polyhedra)";break;
		case pe_UnexpectedRepresentation: 		out << "(Error in polyhedra representation)";break;
		case fe_Open:
		case fe_Parse: 							out << "(Cannot read file correctly)";break;
		case ie_BadIntegrandFormat: 			out << "(Wrong integrand)";break;
		case ie_UnexpectedIntegrationOption: 	out << "(Wrong integration options)";break;
		case de_divisionByZero: 				out << "(Divided by zero, perturbation failed)";break;
		case bug_Unknown: 						out << "(Unknown bug)";break;
		case bug_NotImplementedHere: 			out << "(Feature not yet implemented)";break;

	}
	if ( msg != "")
		out << "\n  Message    : " << msg.c_str();
	out << "\n  File       : " << fileName.c_str()
		 << "\n  Line number: " << lineNumber << endl;
	return out.str();
}

const char* LattException::what() const throw()
{
	string s = "";
	if ( LATTEXCEPTION_PRINTSTATUS || printStatus == 1 )
	s = printErrorMessages();

	return s.c_str();
}
