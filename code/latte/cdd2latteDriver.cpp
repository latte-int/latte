/*
 * cdd2latteDriver.cpp
 *
 *  Created on: Feb 2, 2011
 *      Author: bedutra
 *
 *  usage: given a .ine or .ext file from cdd, we will make the correspoinding latte file.
 *  V-reps are not dilated, but H-reps are.
 *
 *  I needed to convert a cdd h-rep file to latte because latte cannot currently handel rational h-reps.
 */

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cassert>
#include <vector>
#include "gmp.h"
#include <gmpxx.h>
#include <cstdio>

using namespace std;

vector<vector<mpq_class> > readCDDfacets(const char * input)
{
	ifstream in;
	string line;
	mpz_class dim, count;
	int i, j;
	mpq_class term;
	vector<vector<mpq_class> > facets;

	in.open(input);
	if (!in.is_open())
	{
		cout << "Cannot open file " << input << "." << endl;
		exit(1);
	}

	for (getline(in, line, '\n'); line != "begin"; getline(in, line, '\n'))
		;

	in >> count >> dim;
	in >> line; //read in the "integer"/"rationial" text.


	for (i = 0; i < count; ++i)
	{
		vector<mpq_class> oneFacet;

		for (j = 0; j < dim; ++j)
		{
			in >> term;
			oneFacet.push_back(term);
		}//for j
		facets.push_back(oneFacet);//hyperplane has rational coeffs.
	}//while

	in >> line; //read in "end"
	assert ( line == "end");

	return facets;
}//readCDDfacets

/**
 * @parm facets: each row will be dilated.
 * @parm outpfile: file name.
 */
void writeLatteFacets(const char *output, vector<vector<mpq_class> > &list)
{
	ofstream out;
	for (int i = 0; i < (int) list.size(); ++i)
	{
		mpz_t currentLCM;
		mpz_init_set_si(currentLCM, 1);

		//find the lcm of al the denominator for this row only.
		for (int k = 0; k < list[i].size(); ++k)
		{
			if (list[i][k] == mpz_class(0))
				continue;

			mpz_lcm(currentLCM, currentLCM,	list[i][k].get_den_mpz_t());
		}//for k

		assert(mpz_class(currentLCM) > 0);

		//mult each element in this row by the lcm.
		if (mpz_class(currentLCM) != 1)
		{
			for (int k = 0; k < list[i].size(); ++k)
			{
				list[i][k] = list[i][k] * mpz_class(currentLCM);
				assert(list[i][k].get_den() == 1);
			}//for k
		}//divide by the gcd
	}//for i

	//now print it out!
	out.open(output);
	if ( !out.is_open())
	{
		cout << "cannot write to file " << output << endl;
		exit(1);
	}

	out << list.size() << " " << list[0].size() << endl;
	for(int i = 0; i < (int) list.size(); ++i)
	{
		for(int j = 0; j < (int) list[i].size(); ++j)
			out << list[i][j] << " ";
		out << endl;
	}


}//writeLatteFacets


void convertCDDineToLatte(const char * input, const char *output)
{
vector<vector<mpq_class> > facets;
facets = readCDDfacets(input);

writeLatteFacets(output, facets);

}

int main(int argc, char * argv[])
{
if (argc != 3)
{
	cout << "Usage: " << argv[0] << " input-file.[ext | ine] output-file"
			<< endl;
	return 0;
}

string fileName(argv[1]);

if (string::npos != fileName.rfind(".ext"))
{
	cout << "Sorry, this is not implemented" << endl;
}//this is a *.ext* file
else if (string::npos != fileName.rfind(".ine"))
{
	convertCDDineToLatte(argv[1], argv[2]);

}//this is a *.ine* file
else
{
	cout << "Sorry, " << argv[1] << " does not end in .ext or .ine" << endl;
}

return 0;
}//main
