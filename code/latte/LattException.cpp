/*
 * LattException.cpp
 *
 *  Created on: Jul 19, 2011
 *      Author: bedutra
 */

#include "LattException.h"


LattException::~LattException() throw() {}

LattException::LattException(UserError ue, const char * file, const int line, const char * message)
{
	userError = ue;
	fileName = file;
	lineNumber = line;
	msg = message;
}



void LattException::printErrorMessages() const
{
	cout << "\nLatte message    : "
		 << "\n  error code: " << userError
		 << "\n  Optional message: " << msg.c_str()
		 << "\n  file: " << fileName.c_str()
		 << "\n  line number: " << lineNumber << endl;
}

const char* LattException::what() const throw()
{
	printErrorMessages();
	return exception::what();
}
