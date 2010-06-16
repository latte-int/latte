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
	//fileName = "BuildRandomPolytope." + ss.str() + ".temp";
	fileName = "BuildRandomPolytope.temp"; //for now, keep the file names the same (for debugging)
	maxInteger = 200;
	probNegative = .5;
}//buildRandomPolytope

void BuildRandomPolytope::buildPolymakeFile(const int numPoints)
{
	ofstream file;
	
	file.open(fileName.c_str());
	file << "POINTS" << endl;
	
	srand(44); //after testing, change to just srand();
	for(int k = 0; k < numPoints; ++k)
	{
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
	
}//BuildRandomPolytope



void BuildRandomPolytope::callPolymake()
{
	string command;
	command = "polymake " + fileName + " DIM";
	system(command.c_str());
	command = "polymake " + fileName + " FACETS";
	system(command.c_str());
	
//cout << "callPolymake() returning early" << endl;
//return;
	//next, read the polymake file.
	ifstream file;
	string line;
	file.open(fileName.c_str());
	for(getline(file, line, '\n'); line != "DIM"; getline(file, line, '\n'))
		;
	file >> dim;
	
	for(getline(file, line, '\n'); line != "FACETS"; getline(file, line, '\n'))
		;
		
		
	facets.clear();
	
	string term;	
	file >> term;
	do {
		vector<mpq_class> oneFacet;
		for(int i = 0; i < dim + 1; ++i)
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
	

	file.close();
}//callPolymake()

void BuildRandomPolytope::convertFacetEquations()
{
	for(int i = 0; i < (int) facets.size(); ++i)
	{
		mpz_class product(1);
		
		//find the prod. of all the denominators
		for(int k = 0; k < ambientDim; ++k)
		{
			product = product * facets[i][k].get_den();
		}//for k
		
		mpz_class currentGCD(facets[i][0]*product);
		mpz_t ans;
		mpz_init(ans);
		for(int k = 0; k < ambientDim; ++k)
		{
			facets[i][k] = product * facets[i][k];
			mpz_gcd(ans, currentGCD.get_mpz_t(), facets[i][k].get_num_mpz_t());
			currentGCD = mpz_class(ans);	
		}//for k	
		
		if (currentGCD != mpz_class(1))
		{
			for(int k = 0; k < ambientDim; ++k)
			{
				facets[i][k] = facets[i][k] / currentGCD;	
			}//for k	
		}//divide by the gcd
	}//for i
	
	//check everything looks ok.
	for(int i = 0; i < (int) facets.size(); ++i)
	{
		cout << "reduced facet " << i << " = "; 
		for(int k = 0; k < ambientDim; ++k)
			cout << facets[i][k] << ", ";
		cout << endl;
	}//for i

}//convertFacetEquations();


void BuildRandomPolytope::findEhrhardPolynomial()
{
	findFacetEquations();
	
	stringstream d;
	d << dim;
	string command = "./ehrhart2 " + d.str() + " " + latteFile.c_str();
	system(command.c_str());
}//findEhrhardPolynomial


void BuildRandomPolytope::findFacetEquations()
{
	//read the facets
	callPolymake();

	convertFacetEquations();
	
	printFacetEquationsForLattE();
	
	
}//findFacetEquations

void BuildRandomPolytope::printFacetEquationsForLattE()
{
	ofstream file;
	
	latteFile = fileName + ".latte";
	file.open(latteFile.c_str());
	
	file << facets.size() << " " << dim << endl;
	for(int i = 0; i < (int) facets.size(); ++i)
	{
		for(int k = 0; k < ambientDim; ++k)
			file << facets[i][k] << " ";
		file << endl;
	}//for i
	file.close();
}//printFacetEuationsForLatte()
	
