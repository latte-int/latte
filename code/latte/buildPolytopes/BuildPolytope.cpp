/*
 * BuildPolytope.cpp
 *
 *  Created on: June 15, 2010
 *      Author: bedutra
 */
#include "BuildPolytope.h"
#include <cstring>
#include <iostream>
#include <cassert>
#include <ctime>

using namespace std;


BuildPolytope::BuildPolytope():
	ambientDim(0), 
	dim(0), 
	integerPoints(true), 
	createdPolymakeFile(false), createdLatteVRepFile(false), createdLatteHRepFile(false),
	numAffineHull(0)
{
	 time_t rawtime;
	 struct tm * timeinfo;
	 time ( &rawtime );
	 timeinfo = localtime ( &rawtime );

	 stringstream ss;
	 ss << "buildpolytope_";
	 ss << timeinfo->tm_min << "_" << timeinfo->tm_hour << "_" << timeinfo->tm_mday << "_" << timeinfo->tm_year +1990;
	fileBaseName = ss.str();
}

/**
 * Deletes the generated file if it has been created.
 */
void BuildPolytope::deletePolymakeFile() { if (createdPolymakeFile) system((string("rm -f ") + getPolymakeFile()).c_str());}
void BuildPolytope::deleteLatteVRepFile(){ if (createdLatteVRepFile) system((string("rm -f ") + getLatteVRepFile()).c_str());}
void BuildPolytope::deleteLatteHRepFile(){ if (createdLatteHRepFile) system((string("rm -f ") + getLatteHRepFile()).c_str());}

/**
 * Builds the polymake file containing the points information.
 * filename: fileBaseName.polymake.
 */
void BuildPolytope::buildPolymakeFile()
{
	ofstream file;
	
	if ( createdPolymakeFile == true)
		return; //already make polymake file!
		
	createdPolymakeFile = true;

	file.open(getPolymakeFile().c_str());
	file << "POINTS" << endl;

	
	for (int k = 0; k < (int)points.size(); ++k)
	{
		file << 1 << ' ';

		for (int vectorElm = 0; vectorElm < ambientDim; ++vectorElm)
		{
			file << points[k][vectorElm] << ' ';
		}//vectorElm
		file << endl;
	}//for k

	file.close();
}//buildPolymakeFile


/**
 * Prints a file containing the equations of the polytope for Latte.
 * Includes affine hull linearity information if needed.
 * filename: fileBaseName.hrep.latte.
 */
void BuildPolytope::buildLatteHRepFile()
{
	if ( createdLatteHRepFile == true)
		return; //already did this.
		
	createdLatteHRepFile = true;
	findFacets(); //facets information is saved in facets.
	
	//now print to a latte file.
	ofstream file;

	file.open(getLatteHRepFile().c_str());

	file << facets.size() << " " << ambientDim + 1 << endl;
	for (int i = 0; i < (int) facets.size(); ++i)
	{
		for (int k = 0; k < ambientDim + 1; ++k)
			file << facets[i][k] << " ";
		file << endl;
	}//for i

	if (numAffineHull > 0)
	{
		file << "linearity " << numAffineHull << " ";
		for (int i = (int)facets.size() - numAffineHull; i < (int)facets.size(); ++i)
		{
			file << i + 1 << " ";
		}//for the affine hull equations.
		file << endl;
	}//if there are affine hulls!
	file.close();}

/**
 * Prints a file containing the vertices of the polytope for Latte.
 * filename: fileBaseName.vrep.latte.
 */
void BuildPolytope::buildLatteVRepFile()
{
	if ( createdLatteVRepFile == true)
		return; //already did this.
		
	createdLatteVRepFile = true;
	findVertices(); //facets information is saved in facets.
	
	//now print to a latte file.
	ofstream file;

	file.open(getLatteVRepFile().c_str());

	file << points.size() << " " << ambientDim + 1 << endl;
	for (int i = 0; i < (int) points.size(); ++i)
	{
		for (int k = 0; k < ambientDim + 1; ++k)
			file << points[i][k] << " ";
		file << endl;
	}//for i

	file << endl;
	file.close();
}


/*
 * Assumes the polymake file exist and polymake has already ran.
 */
void BuildPolytope::findDimentions()
{
	if ( dim > 0)
		return; //already found it, don't waste time.

	ifstream file;
	string line;
	file.open(getPolymakeFile().c_str());
	for (getline(file, line, '\n'); line != "DIM"; getline(file, line, '\n'))
		;
	file >> dim;
	file.close();

	assert(0 < dim && dim <= ambientDim);
}

/**
 * Reads the polymake facet information into the facet vector.
 * Builds the polymake file if needed.
 */
void BuildPolytope::findFacets()
{
	buildPolymakeFile();
	
	//ask polymake for the 
	system((string("polymake ") + getPolymakeFile() + " DIM AMBIENT_DIM FACETS AFFINE_HULL").c_str());

	//next, read the polymake file.

	findDimentions();
	ifstream file;
	string line;
	file.open(getPolymakeFile().c_str());
	for (getline(file, line, '\n'); line != "FACETS"; getline(file, line, '\n'))
		;

	facets.clear();

	string term;
	file >> term;
	do
	{				
		vector<mpq_class> oneFacet;
		for (int i = 0; i < ambientDim + 1; ++i)
		{
			oneFacet.push_back(mpq_class(term));
			file >> term;
			//cout << term << ". " ;
		}//for i

		facets.push_back(oneFacet);//hyperplane has interger coeffs.
	} while (term != "AFFINE_HULL");


	file.ignore(10000, '\n');

	getline(file, line, '\n');
	numAffineHull = 0;
	while (line != "")
	{
		stringstream ss(line);
		vector<mpq_class> oneFacet;

		for (int k = 0; k < ambientDim + 1; ++k)
		{
			ss >> term;
			oneFacet.push_back(mpq_class(term));
		}//for k
		facets.push_back(oneFacet);
		numAffineHull += 1;
		getline(file, line, '\n');
	}//while

	file.close();	
	convertFacetEquations(); //make coeffs. integer.
}//findFacets

/**
 * Given array of facets, will reduce the elements in the vectors from rationals to integers
 */
void BuildPolytope::convertFacetEquations()
{

	for (int i = 0; i < (int) facets.size(); ++i)
	{
		mpz_class product(1);

		//find the prod. of all the denominators
		for (int k = 0; k < ambientDim + 1; ++k)
		{
			product = product * facets[i][k].get_den();
		}//for k

		mpz_class currentGCD(facets[i][0] * product);
		mpz_t ans;
		mpz_init(ans);
		for (int k = 0; k < ambientDim + 1; ++k)
		{
			facets[i][k] = product * facets[i][k];
			if (facets[i][k] != mpq_class(0))
			{
				mpz_gcd(ans, currentGCD.get_mpz_t(),
						facets[i][k].get_num_mpz_t());
				currentGCD = mpz_class(ans);
			}//if
		}//for k	

		if (currentGCD != mpz_class(1))
		{
			for (int k = 0; k < ambientDim + 1; ++k)
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


/**
 * Reads the polymake vertex information into the points vector.
 * Builds the polymake file if needed.
 */
void BuildPolytope::findVertices()
{
	buildPolymakeFile();
	
	//ask polymake for the 
	system((string("polymake ") + getPolymakeFile() + " DIM AMBIENT_DIM VERTICES ").c_str());


	//next, read the polymake file.
	findDimentions();
	ifstream file;
	string line;
	file.open(getPolymakeFile().c_str());

	for (getline(file, line, '\n'); line != "VERTICES"; getline(file, line, '\n'))
		;

	points.clear(); //remove old points information.

	string term;
	file >> term;
	do
	{				
		vector<mpq_class> onePoint;
		for (int i = 0; i < ambientDim + 1; ++i)
		{
			onePoint.push_back(mpq_class(term));
			file >> term;
		}//for i

		points.push_back(onePoint);//hyperplane has interger coeffs.
	} while ( strcmp(term.c_str(), "1") == 0);

	file.close();	

}
	
//int BuildPolytope::gcd(int a, int b) const;
int BuildPolytope::getAmbientDim() const {return ambientDim;}
int BuildPolytope::getDim() const {return dim;}
	
string BuildPolytope::getLatteVRepFile() const {return fileBaseName + ".vrep.latte";}
string BuildPolytope::getLatteHRepFile() const { return fileBaseName + ".hrep.latte";}
string BuildPolytope::getPolymakeFile() const { return fileBaseName + ".polymake";}

	
	
void BuildPolytope::setBaseFileName(const string & n) {fileBaseName = n;}
void BuildPolytope::setIntegerPoints(bool t) { integerPoints = t; }


void BuildPolytope::forDebugging()
{
	ambientDim = 3;
	
	for (int i = 0; i < ambientDim*4; ++i)
	{
		vector<mpq_class> onePoint;
		for(int j = 0; j < ambientDim; ++j)
			onePoint.push_back(mpq_class(rand()%10, rand()%25));
		points.push_back(onePoint);
		
	}
}





