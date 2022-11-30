/**
 * To compile: g++ -o make-cube make-cube.cpp
 */
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>

using namespace std;


struct QuickFraction
{
	long numerator, denominator;
	QuickFraction(const string & number)
	{
		for (int i = 0; i < number.length(); ++i)
			if (number[i] == '/')
			{
				numerator = atol(number.substr(0, i).c_str());
				denominator = atol(
						number.substr(i + 1, number.length() - i - 1).c_str());
				return;
			}
		numerator = atol(number.c_str());
		denominator = 1;

	}
};


/**
 * Makes a cube, [0, scale]^n, as a latte facet file.
 * @param filename: output file name
 * @param dim: dim = n from above.
 * @param scale: we scale the unit cube.
 */
void makeCubeFile(const char *fileName, const long &dim, const QuickFraction& scale)
{
	ofstream out;
	long b, x;
	b = scale.numerator;
	x = scale.denominator; //coeff. of the x_ith term. 0 < x_i < b/x
	

	out.open(fileName);

	out << 2 * dim << " " << dim + 1 << endl; //first line.

	for(int i = 0; i < dim; ++i)
	{
		// 0 < b 0 0 0 0 -x
		out << b << " ";
		for(int k = 0; k < dim; ++k)
			if ( k == i)
				out << -1 * x << " ";
			else
				out << "0 ";
		out << endl;

		// 0 < 0 0000 1
		out << "0 ";
		for(int k = 0; k < dim; ++k)
			if ( k == i)
				out << "1 ";
			else
				out << "0 ";
		out << endl;
	}//for i.
	out.close();
}//makeCubeFile




int main(const int argc, const char **argv )
{
	if ( argc <= 1)
	{
		cout << "usage: " << argv[0] << " <output file name> <dim> [scale factor] " << endl;
	}
	if ( argc == 4 )
		makeCubeFile(argv[1], atol(argv[2]), QuickFraction(argv[3]));
	else
		makeCubeFile(argv[1], atol(argv[2]), QuickFraction("1"));

	return 0;
}



