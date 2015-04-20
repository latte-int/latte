/*
 * Driver.cpp
 *
 *  Created on: Jun 24, 2010
 *      Author: Brandon Dutra and Gregory Pinto
 *
 *
 *
 *  type --help to print the help menu.
 */
#include <iostream>
#include "valuation.h"


int main(int argc, char *argv[])
{
	try {
		Valuation::mainValuationDriver((const char **) argv, argc);
	} catch ( LattException & e)
	{
		std::cout << e.what() << std::endl;
		return 1;
	}


	return 0;
}
//main()
