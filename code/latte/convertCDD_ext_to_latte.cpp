
/*
 * Converts a CDD style H-rep to a latte style H-rep.
 * Why this is needed: cdd's H-rep allow fractional coeff., whereas latte does not!
 */

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cassert>
#include <vector>
#include "rational.h"
#include <NTL/ZZ.h>


using namespace std;

int main( int argc, const char* argv[] )

{
	if ( argc != 3)
	{
		cout << "usage: " << argv[0] << "inputFileName.ine outputFileName" << endl;
		exit(1);
	}
	cout << "not implemented yet" << endl;
	exit(1);

	ifstream input(argv[1]);
	ofstream output(argv[2]);
	string line;
	string temp;
	long numRows, ambDim;


	while ( getline(input, line))
		if ( line == "begin")
			break;
//	getline(input, line);

	input >> numRows >> ambDim;
	input.ignore(250, '\n');
	output << numRows << " " << ambDim << endl;
//	hRep.resize(ambDim);

	for(int i = 0; i < numRows; ++i )
	{
		for(int k  = 0; k < ambDim; ++k)
		{

			input >> temp;
			output << temp << " ";
			//hRep[k] = RationalNTL(temp);
			//RationalNTL product;
			//lcmDenom = (hRep[k].getDenominator() * lcmDenom)/ (GCD(lcmDenom, hRep[k].getDenominator()));
		}
		//for(int k = 0; k < ambDim; ++k)
		//	output << (hRep[k] * lcmDenom) << " " ;
		output << endl;

	}//for each row.
	input >> temp;
	assert(temp == "end");

	input.close();
	output.close();
	return 0;
}
