/*
 * BuildRandomPolytope.cpp
 *
 *  Created on: June 15, 2010
 *      Author: bedutra
 */
#include "BuildRandomPolytope.h"
#include <iostream>

using namespace std;


BuildRandomPolytope::BuildRandomPolytope(int ambient_dim)
{
	ambientDim = ambient_dim;
	stringstream ss;
	ss << time(0);
	fileName = "BuildRandomPolytope." + ss.str() + ".temp";
	//fileName = "BuildRandomPolytope.temp"; //for now, keep the file names the same (for debugging)
	PolytopeComments = "Random edge polytope"; //default label. caller should call setComments for a better label.
	maxInteger = 50;
	probNegative = .5;
	numAffineHull = 0;

}//buildRandomPolytope


BuildRandomPolytope::~BuildRandomPolytope()
{
	stringstream command;
	command << "rm -f " << fileName.c_str();
	cout << "Deleting file " << command.str().c_str() << endl;
	system(command.str().c_str());
	command.str("");
	command.clear();
	command << "rm -f " << latteFile.c_str();
	cout << "Deleting file " << command.str().c_str() << endl;
	system(command.str().c_str());
}//deconstructor.

/**
 * Generates random points.
 * @parm int numPoints: number of random points to include.
 */
void BuildRandomPolytope::buildPolymakeFile(const int numPoints)
{
	ofstream file;
	
	file.open(fileName.c_str());
	file << "POINTS" << endl;
	
	srand(time(0)); //after testing, change to just srand();
	for(int k = 0; k < numPoints; ++k)
	{
		file << 1 << ' ';
		for(int vectorElm = 0; vectorElm < ambientDim; ++vectorElm)
		{

			int number = rand() % maxInteger;
			if ( rand() < RAND_MAX * probNegative)
				number = -1 * number;
			file << number << ' ';
		}//vectorElm
		file << endl;
	}//for k
	
	file.close();
	
	//ambientDim += 1;

}//BuildRandomPolytope


/**
 * Calls polymake to find the dim., ambient dim., and facets of the polytope.
 */
void BuildRandomPolytope::callPolymake()
{
	string command;
	cout << "callPolymake(): starting to call polymake" << endl;
	command = "polymake " + fileName + " DIM AMBIENT_DIM";
	system(command.c_str());
	command = "polymake " + fileName + " FACETS";
	system(command.c_str());
	cout << "callPolymake(): polymake finished, now reading from file" << endl;
	
	//next, read the polymake file.
	ifstream file;
	string line;
	file.open(fileName.c_str());
	for(getline(file, line, '\n'); line != "DIM"; getline(file, line, '\n'))
		;
	file >> dim;
	
	for(getline(file, line, '\n'); line != "AMBIENT_DIM"; getline(file, line, '\n'))
		;
	file >> ambientDim;

	for(getline(file, line, '\n'); line != "FACETS"; getline(file, line, '\n'))
		;
		
		
	facets.clear();
	
	string term;	
	file >> term;
	do {
		vector<mpq_class> oneFacet;
		for(int i = 0; i < ambientDim+1; ++i)
		{
			oneFacet.push_back(mpq_class(term));
			file >> term;
		}//for i
		//cout << "one facet =";
		//for(int i = 0; i < (int) oneFacet.size(); ++i)
		//	cout << oneFacet[i] << ", ";
		//cout << endl;
		facets.push_back(oneFacet);
	} while ( term != "AFFINE_HULL");
	
	file.ignore(10000, '\n');

	getline(file, line, '\n');
	while ( line != "")
	{
		stringstream ss(line);
		vector<mpq_class> oneFacet;

		for(int k = 0; k < ambientDim+1; ++k)
		{
			ss >> term;
			oneFacet.push_back(mpq_class(term));
		}//for k
		facets.push_back(oneFacet);
		numAffineHull += 1;
		getline(file, line, '\n');
	}//while

	file.close();
}//callPolymake()

/**
 * Given array of facets, will reduce the elements in the vectors from rationals to integers
 */
void BuildRandomPolytope::convertFacetEquations()
{


	for(int i = 0; i < (int) facets.size(); ++i)
	{
		mpz_class product(1);
		
		//find the prod. of all the denominators
		for(int k = 0; k < ambientDim+1; ++k)
		{
			product = product * facets[i][k].get_den();
		}//for k
		
		mpz_class currentGCD(facets[i][0]*product);
		mpz_t ans;
		mpz_init(ans);
		for(int k = 0; k < ambientDim+1; ++k)
		{
			facets[i][k] = product * facets[i][k];
			if ( facets[i][k] != mpq_class(0))
			{
				mpz_gcd(ans, currentGCD.get_mpz_t(), facets[i][k].get_num_mpz_t());
				currentGCD = mpz_class(ans);
			}//if
		}//for k	
		
		if (currentGCD != mpz_class(1))
		{
			for(int k = 0; k < ambientDim+1; ++k)
			{
				facets[i][k] = facets[i][k] / currentGCD;	
			}//for k	
		}//divide by the gcd
	}//for i
	
	//check everything looks ok.
	//for(int i = 0; i < (int) facets.size(); ++i)
	//{
	//	cout << "reduced facet " << i << " = ";
	//	for(int k = 0; k < ambientDim+1; ++k)
	//		cout << facets[i][k] << ", ";
	//	cout << endl;
	//}//for i

}//convertFacetEquations();


const string & BuildRandomPolytope::getLatteFile() const
{
	return latteFile;
}//getLatteFile();


/**
 * Assume: buildPolymakeFile() was called first.
 * Calls ehrhart2 to get and print the Ehrhart poly!
 */
void BuildRandomPolytope::findEhrhardPolynomial()
{
	findFacetEquations();
	
	stringstream d;
	d << dim;
	string command = "./ehrhart2 " + d.str() + " " + latteFile.c_str() + " -Rs " + "\"" + PolytopeComments.c_str()+ "\"";
	cout << "findEhrhardPolynomial(): ehrhart calling: " << command.c_str() << endl;
	system(command.c_str());
}//findEhrhardPolynomial

/**
 * Assume: buildPolymakeFile() was called first.
 * Runs polymake to find the facets and then extracts then for lattE
 */
void BuildRandomPolytope::findFacetEquations()
{
	//read the facets
	callPolymake();

	convertFacetEquations();
	
	printFacetEquationsForLattE();
	
}//findFacetEquations


/**
 * Prints the polytope's volume. We do not currently save or use this information.
 * This method allows us to check our valuation/volume computation in /valuation
 */
void BuildRandomPolytope::findVolumeWithPolymake()
{
	stringstream command; //used stringstream instead of sting because of a weird compile error.
	command << "polymake " << fileName.c_str() << " VOLUME ";
	system(command.str().c_str());

}//findVolumeWithPolymake()

/**
 * Prints a file containing the equations of the polytope for Latte.
 * Includes affine hull linearity information if needed.
 */
void BuildRandomPolytope::printFacetEquationsForLattE()
{
	ofstream file;
	
	latteFile = fileName + ".latte";
	file.open(latteFile.c_str());
	
	file << facets.size() << " " << ambientDim+1 << endl;
	for(int i = 0; i < (int) facets.size(); ++i)
	{
		for(int k = 0; k < ambientDim+1; ++k)
			file << facets[i][k] << " ";
		file << endl;
	}//for i

	if ( numAffineHull > 0)
	{
		file << "linearity " << numAffineHull << " " ;
		for(int i = facets.size() - numAffineHull; i < facets.size(); ++i)
		{
			file << i + 1 << " ";
		}//for the affine hull equations.
		file << endl;
	}//if there are affine hulls!
	file.close();
}//printFacetEuationsForLatte()

/**
 * This string is passed to ehrhart2 and is printed in the log files to make things easier to read.
 */
void BuildRandomPolytope::setComments(const string& newComments)
{
	PolytopeComments = newComments;
}//setComments

	
