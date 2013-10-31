#include <iostream>
#include <fstream>
#include <vector>
#include "TopKnapsack.h"
#include "latte_ntl.h"
#include "integration/burstTrie.h"
#include "integration/multiply.h"
#include <cstdlib>
#include <ctime>
#include <climits>

#include "cone.h"
#include "preprocess.h"
#include "barvinok/barvinok.h"
#include "barvinok/dec.h"
#include "integration/multiply.h"
#include "PeriodicFunction.h"
#include "print.h"
#include "dual.h"
#include "timing.h"
#include <sstream>
#include <getopt.h>
using namespace std;

int main(int argc, char *argv[]) {
	if (1) {

		cout << "Invocation: ";
		for(int i = 0; i < argc; ++i)
			cout << argv[i] << " ";
		cout << endl;

		string inFile, outFile;
		int k = -1;

		//process each option.
		while (1)
		{
			static struct option long_options[] =
			{
			{ "out",	  			  required_argument, 0, 'o' },
			{ "k",  				  required_argument, 0, 'k' },
			{ "help",				  no_argument,       0, 'h' },
			{ "file",				  required_argument, 0, 'f' },
			{ 0, 0, 0, 0 } };
			/* getopt_long stores the option index here. */

			int option_index = 0;
			int c;

			//single-character=short option-name. x: means x takes a required argument, x:: means x takes an optional argument.
			c = getopt_long(argc, argv, "o:k:hf:", long_options, &option_index);

			if (c == -1)
				break;

			switch (c)
			{
				case 0:
					// If this option set a flag, do nothing
					break;
				case 'o':
					outFile = optarg;
					break;
				case 'k':
					k = atoi(optarg);
					break;
				case 'h':
					cout << "TODO: print help menu" << endl;
					break;
				case 'f':
					inFile = optarg;
					break;
				default:
					cout << "main: Unknown case" << endl; //todo: throw exception in my fancy class.
					exit(1);
			}
		}//while.

		if (inFile.length() == 0 || k < 0)
		{
			cout << "Input error: run " << argv[0] << " -h for help" << endl;
			exit(1);
		}

		ifstream file;
		file.open(inFile.c_str());
		int n;
		file >> n;
		vec_ZZ alpha;
		alpha.SetLength(n);

		for (int i = 0; i < n; ++i)
			file >> alpha[i];

		Timer time("Total");
		time.start();
		TopKnapsack tk;
		tk.set(alpha);
		tk.coeff_NminusK(k);
		time.stop();

		if ( outFile.length())
		{
			ofstream f(outFile.c_str());
			tk.printAnswer(f);
			f << "#Total Time: " << time.get_seconds() << endl;
			f.close();
		}
		else
		{
			tk.printAnswer(cout);
			cout << "Total Time: " << time.get_seconds() << endl;
		}

	}
	else if (0)
	{
		TopKnapsack tk;

		GeneralMonomialSum<PeriodicFunction, int> a;
		vector<ZZ> ea, ee;
		ea.push_back(to_ZZ(1)); ee.push_back(to_ZZ(2));
		ea.push_back(to_ZZ(2)); ee.push_back(to_ZZ(0));
		ea.push_back(to_ZZ(3)); ee.push_back(to_ZZ(4));
		ea.push_back(to_ZZ(4)); ee.push_back(to_ZZ(5));

		cout << "got here" << endl;
		vec_ZZ junk;
		junk.SetLength(ea.size());


		tk.set(junk);
		tk.order = 4;
		ZZ bottom;
		tk.expandPeriodicPart(bottom, a, 1, ea, ee);

		cout << "anser is sign/bottom*" << a.printMonomials().c_str() << endl;

	}
}
//main()
