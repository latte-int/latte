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
    "  --help, -h                               Prints this help message\n"
    "  --out, -o FILENAME                       Saves the Ehrhart polynomial to a file.\n"   
    "\nExamples:\n"
    "  ./top-ehrhart-knapsack -o results.mpl -k 1 -f partition.knap;\n"
    "  (Will compute the top coefficient of the knapsack given in partition.knap and save the answer to a text file called results.mpl. \n"
    << endl;

}

int main(int argc, char *argv[]) {

	//ZZ aa(NTL::INIT_VAL,   "533765551555708028075229384988162533893304992101320918011350006754330831749448685737984256335377984983889916831716778019758878193530263849346252117733839462400");
	//ZZ aa(NTL::INIT_VAL,   "4520293857002398570928708520349769486720967028972394570923857023895709283562385793465120893570294630160562306032857023896209857857203462304985702349860234950");
	//ZZ aa(NTL::INIT_VAL, "4609533765551555708028075229384988162533893304992101320918011350006754330831749448685737984256335377984983889916831716778019758878193530263849346252117733839462400");
	//ZZ bb(NTL::INIT_VAL, "-14587749093744857173029509001388066559542146001301398349519941981826976622712267583318429498029400");
	//cout << "aa:=" << aa << ";\nbb:=" << bb << ";" << endl;

	//ZZ g;
	//g = GCD(abs(aa),abs(bb));
	//cout << "g:=" << myGCD(aa,bb) << ";" << endl;
	//cout << "YES!!!\n" << g <<  endl;

	//cout << "hi" << endl;
	//mpz_class aaa("4609533765551555708028075229384988162533893304992101320918011350006754330831749448685737984256335377984983889916831716778019758878193530263849346252117733839462400",10);
	//mpz_class aaa = convert_ZZ_to_mpz(aa);
	//cout << "world" << endl;
	//mpz_class bbb("-14587749093744857173029509001388066559542146001301398349519941981826976622712267583318429498029400",10);
	//mpz_class bbb = convert_ZZ_to_mpz(bb);
	//cout << "aaa:=" << flush << aaa << flush << ";\nbbb:=" << bbb << ";" << endl;
	//mpz_class thegcd;
	//mpz_gcd(thegcd.get_mpz_t(), aaa.get_mpz_t(),bbb.get_mpz_t());
	//cout << "YES!" << endl;

	//mpz_t at, bt;
	//mpz_init(at);
	//mpz_init(bt);
	//cout << "hi" << endl;
	//mpz_set_str(at, "4609533765551555708028075229384988162533893304992101320918011350006754330831749448685737984256335377984983889916831716778019758878193530263849346252117733839462400", 10);
	//cout << "world" << endl;
	//mpz_set_str(bt, "-14587749093744857173029509001388066559542146001301398349519941981826976622712267583318429498029400", 10);
	//cout << "b = " << flush; mpz_out_str(stdout, 10, bt); cout << "\na = " << flush; mpz_out_str(stdout,10,at); cout << endl;
	if (1) {

		cout << "Invocation: ";
		for(int i = 0; i < argc; ++i)
			cout << argv[i] << " ";
		cout << endl;

		string inFile, outFile;
		int k = -1;
		int allk = -1;
		int seed = 654354;

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
			f << "\n\n\n#Total Time: " << time.get_seconds() << endl;
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
		ea.push_back(to_ZZ(2)); ee.push_back(to_ZZ(1));
		ea.push_back(to_ZZ(3)); ee.push_back(to_ZZ(4));
		ea.push_back(to_ZZ(4)); ee.push_back(to_ZZ(5));

		PeriodicFunction pf;
		pf.setToConstant(20);
		for(int j = 0; j < ea.size(); ++j)
		{
			pf.add(PeriodicFunction(RationalNTL(ea[j], ee[j]), false));
		}
		pf.pow(34);
		PeriodicFunction pf2(pf);
		pf2.times(pf);
		
		cout << pf2 << endl;
		
	}
}
//main()
