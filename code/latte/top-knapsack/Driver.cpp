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

#include "rational.h"
#include "cone.h"
#include "preprocess.h"
#include "barvinok/barvinok.h"
#include "barvinok/dec.h"
#include "integration/multiply.h"
#include "PeriodicFunction.h"
#include "print.h"
#include "dual.h"
#include "timing.h"
#include "config.h"
#include <sstream>
#include <getopt.h>
using namespace std;


void printHelpMenu()
{
	cout << 
	"Required parameters:\n"
	"  --file, -f FILENAME.{knap}               Input knapsack file.\n"
    "Options that control what to compute:\n"
    "  -k INT                                   Computes the INT-th coefficient: T^{N-INT+1}\n"
    "  --all-k INT                              Computes the top INT many coefficients: T^N, ..., T^{N-INT+1}\n"
    "Other options:\n"
    "  --random INT                             Seeds srand with INT or time(0) if INT = -1\n"
    "  --gcd-polynomial [0|1]                   If 1, uses a polynomial time algorithm in k for finding poles. Default is 1.\n"
    "  --help, -h                               Prints this help message\n"
    "  --out, -o FILENAME                       Saves the Ehrhart polynomial to a file.\n"   
    "\nExamples:\n"
    "  ./top-ehrhart-knapsack -o results.mpl -k 1 -f partition.knap;\n"
    "  (Will compute the top coefficient of the knapsack given in partition.knap and save the answer to a text file called results.mpl. \n"
    << endl;

}

int main(int argc, char *argv[]) {

	cout << "Invocation: ";
	for(int i = 0; i < argc; ++i)
		cout << argv[i] << " ";
	cout << endl;

	string inFile, outFile;
	int k = -1;
	int allk = -1;
	int seed = 654354;
	bool gcdPolynomial = true;

	//process each option.
	while (1)
	{
		static struct option long_options[] =
		{
		{ "out",	  			  required_argument, 0, 'o' },
		{ "k",  				  required_argument, 0, 'k' },
		{ "all-k",				  required_argument, 0, 0x100 }, //hex value larger than 1 byte/char
		{ "help",				  no_argument,       0, 'h' },
		{ "file",				  required_argument, 0, 'f' },
		{ "random",				  required_argument, 0, 'r' },
		{ "gcd-polynomial",           required_argument, 0, 0x101},
		{ 0, 0, 0, 0 } };
		/* getopt_long stores the option index here. */

		int option_index = 0;
		int c;

		//single-character=short option-name. x: means x takes a required argument, x:: means x takes an optional argument.
		c = getopt_long(argc, argv, "o:k:hf:r:", long_options, &option_index);

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
				k = atoi(optarg)-1;
				break;
			case 0x100:
				allk = atoi(optarg)-1;
				break;
			case 'h':
				printHelpMenu();
				exit(0);
				break;
			case 'f':
				inFile = optarg;
				break;
			case 'r':
				seed = atoi(optarg);
				break;
			case 0x101:
				gcdPolynomial = atoi(optarg);
				break;
			default:
				cout << "main: Unknown case" << endl; //todo: throw exception in my fancy class.
				exit(1);
		}
	}//while.

	if (inFile.length() == 0 || (k < 0 && allk < 0))
	{
		cout << "Input error: run " << argv[0] << " -h for help" << endl;
		THROW_LATTE( LattException::ue_BadCommandLineOption);
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
	tk.useSubsetsForGCD(gcdPolynomial);
	tk.seed(seed);
	if ( k > -1)
		tk.coeff_NminusK(k);
	else
		tk.coeff_topK(allk);
	time.stop();

	cout << "Printing answer..." << endl;
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
//main()
