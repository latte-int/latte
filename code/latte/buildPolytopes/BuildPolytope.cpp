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
#include <unistd.h>

using namespace std;


BuildPolytope::BuildPolytope():
	ambientDim(0), 
	dim(0), 
	integerPoints(true), 
	createdPolymakeFile(false), createdPolymakeDualFile(false),
	createdLatteVRepFile(false), createdLatteHRepFile(false),
	createdLatteHRepDualFile(false), createdLatteVRepDualFile(false),
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
 * Takes a point and adds it to our point vector.
 * @onePoint: one point to add. Assumes this does not have a leading 1.
 * We will add [1 onePoint] to the point vector
 */
void BuildPolytope::addPoint(vector<mpq_class> onePoint)
{
	onePoint.insert(onePoint.begin(), 1);
	points.push_back(onePoint);
}//addPoint


/**
 * Checks to see if this polytope is centered, and then
 * centers it if not.
 * The center command keeps important properties like vertices and facets, but some are
 * lost like SIMPLE and SIMPLICIAL.
 *
 * The origional un-centered polytope is lost. Also, we do not automatically re-read vertex and facet information.
 */
void BuildPolytope::centerPolytope()
{
	if ( isCentered() )
		return;
	//I do not know why, but in the terminal I can type center f.polymake f.polymake,
	//but here, it does not like it when the output and input files are the same.
	//Hence the ".temp" file that is made and renamed.
	system((string("center ") + getPolymakeFile() + ".temp " + getPolymakeFile()).c_str());
	rename((string(getPolymakeFile()+".temp ")).c_str(), getPolymakeFile().c_str());
	points.clear();
	facets.clear();
	dualFacets.clear();
	dualVertices.clear();
}

void BuildPolytope::clearPoints()
{
	points.clear();
}

/**
 * Deletes the generated file if it has been created.
 */
void BuildPolytope::deletePolymakeFile() { if (createdPolymakeFile) unlink(getPolymakeFile().c_str());}
void BuildPolytope::deletePolymakeDualFile() { if ( createdPolymakeDualFile) unlink(getPolymakeDualFile().c_str());}
void BuildPolytope::deleteLatteVRepFile(){ if (createdLatteVRepFile) unlink(getLatteVRepFile().c_str());}
void BuildPolytope::deleteLatteVRepDualFile(){ if (createdLatteVRepDualFile) unlink(getLatteVRepDualFile().c_str());}
void BuildPolytope::deleteLatteHRepFile(){ if (createdLatteHRepFile) unlink(getLatteHRepFile().c_str());}
void BuildPolytope::deleteLatteHRepDualFile(){ if (createdLatteHRepDualFile) unlink(getLatteHRepDualFile().c_str());}
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
		for (int vectorElm = 0; vectorElm < ambientDim+1; ++vectorElm)
		{
			file << points[k][vectorElm] << ' ';
		}//vectorElm
		file << endl;
	}//for k

	file.close();
}//buildPolymakeFile


/**
 * Finds the dual vertices and prints them to a new polymake file.
 * If this polytope is not centered, will center it.
 */
void BuildPolytope::buildPolymakeDualFile()
{

	if (createdPolymakeDualFile == true)
		return;

	findVerticesDual();//make sure we already know the dual vertices.
	fstream file;
	file.open(getPolymakeDualFile().c_str(), ios::out);
	file << "VERTICES" << endl;
	for(int i = 0; i < (int) dualVertices.size(); ++i)
	{
		for(size_t k = 0; k < dualVertices[i].size(); ++k)
			file << dualVertices[i][k] << " ";
		file << endl;
	}//for i.

	createdPolymakeDualFile = true;
	file.close();
}

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
	makeIntegerRows(facets);
	
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
 * Makes a latte dual h-rep file.
 */
void BuildPolytope::buildLatteHRepDualFile()
{
	centerPolytope();
	findFacetsDual();

	if ( createdLatteHRepDualFile == true)
		return; //already did this.

	createdLatteHRepDualFile = true;
	findFacetsDual(); //facets information is saved in facets.
	makeIntegerRows(dualFacets);

	//now print to a latte file.
	ofstream file;

	file.open(getLatteHRepDualFile().c_str());


	file << dualFacets.size() << " " << ambientDim + 1 << endl;
	for (int i = 0; i < (int) dualFacets.size(); ++i)
	{
		for (int k = 0; k < ambientDim + 1; ++k)
			file << dualFacets[i][k] << " ";
		file << endl;
	}//for i

	file.close();
}//buildLatteHRepDualFile

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
	makeIntegerList(points);
	
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

/**
 * Makes a latte dual v-rep from the facets information.
 */
void BuildPolytope::buildLatteVRepDualFile()
{
	if ( createdLatteVRepDualFile == true)
		return; //already did this.

	createdLatteVRepDualFile = true;
	findVerticesDual();
	makeIntegerList(dualVertices);//after this, dualVertices are in the form [a, vij], vij are integer

	//now print to a latte file.
	ofstream file;

	file.open(getLatteVRepDualFile().c_str());


	file << dualVertices.size() << " " << ambientDim + 1 << endl;
	for (int i = 0; i < (int) dualVertices.size(); ++i)
	{
		file << dualVertices[i][0] << " ";
		for (int k = 1; k < ambientDim + 1; ++k)
			file << dualVertices[i][k] << " ";
		file << endl;
	}//for i

	file << endl;
	file.close();
}//buildLatteVRepFualFIle()

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

	file.open(getPolymakeFile().c_str());
	for (getline(file, line, '\n'); line != "AMBIENT_DIM"; getline(file, line, '\n'))
		;
	file >> ambientDim;
	file.close();

	assert(0 < dim && dim <= ambientDim);
}

/**
 * Reads the polymake facet information into the facet vector.
 * Builds the polymake file if needed.
 * @rationalize: default is true. If true, each facet equation is mult. by a different positive number s.t. each slot is integer.
 */
void BuildPolytope::findFacets()
{
	buildPolymakeFile();
	
	//ask polymake for the 

	system((string("polymake ") + getPolymakeFile() + " DIM AMBIENT_DIM FACETS AFFINE_HULL > /dev/null ").c_str());

	//next, read the polymake file.

	findDimentions();
	ifstream file;
	string line;
	file.open(getPolymakeFile().c_str());
	for (getline(file, line, '\n'); line != "FACETS"; getline(file, line, '\n'))
		;

	facets.clear();

	getline(file, line, '\n');
	string term;
	while (line != "")
	{
		stringstream ss(line);
		vector<mpq_class> oneFacet;

		for (int k = 0; k < ambientDim + 1; ++k)
		{
			ss >> term;
			oneFacet.push_back(mpq_class(term));
		}//for k
		facets.push_back(oneFacet);//hyperplane has rational coeffs.
		getline(file, line, '\n');
	}//while
	/*
	 *
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

		facets.push_back(oneFacet);
	} while (! file.eof() && isStringNumber(term)); // != "AFFINE_HULL");

	*/
	file.close();
	findAffineHull(); //we need to close and reopen the file to find this because we do not know which one (facets or affineHull) comes first, and I'm not worried about the time cost of doing this.


}//findFacets


/**
 * Assumes facets and amb. dim has already been found and read in.
 * inserts the equations at the end of the facets vector.
 */

void BuildPolytope::findAffineHull()
{
	ifstream file;
	string line;
	string term;
	file.open(getPolymakeFile().c_str());
	for (getline(file, line, '\n'); line != "AFFINE_HULL"; getline(file, line, '\n'))
			;

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
}//findAffineHull

void BuildPolytope::findFacetsDual()
{
	centerPolytope();
	dualFacets = getVertices();
}//findFacetsDual()





/**
 * Reads the polymake vertex information into the points vector.
 * Builds the polymake file if needed.
 */
void BuildPolytope::findVertices()
{
	buildPolymakeFile();
	
	//ask polymake for the 
	system((string("polymake ") + getPolymakeFile() + " DIM AMBIENT_DIM VERTICES > /dev/null ").c_str());


	//next, read the polymake file.
	findDimentions();
	ifstream file;
	string line;
	file.open(getPolymakeFile().c_str());

	for (getline(file, line, '\n'); line != "VERTICES"; getline(file, line, '\n'))
		;

	points.clear(); //remove old points information.

	getline(file, line, '\n');
	string term;
	while (line != "")
	{
		stringstream ss(line);
		vector<mpq_class> oneRay;

		for (int k = 0; k < ambientDim + 1; ++k)
		{
			ss >> term;
			oneRay.push_back(mpq_class(term));
		}//for k
		points.push_back(oneRay);//points are rational.
		getline(file, line, '\n');
	}//while
	file.close();	
}

/**
 * finds the dual vertices.
 * After this call, dualVertices have the form [1, vi], vi are rational.
 */
void BuildPolytope::findVerticesDual()
{
	if ( dualVertices.size() > 0)
		return;

	centerPolytope();
	assert( getDim() == getAmbientDim() ); //duals are only defined for full dim.

	findFacets();

	dualVertices = getFacets();

	homogenizeDualVertices(); //after this, dualVertices are in the form [1, vij], vij are rational
}


//int BuildPolytope::gcd(int a, int b) const;
int BuildPolytope::getAmbientDim() const {return ambientDim;}
int BuildPolytope::getDim() const {return dim;}

vector<vector<mpq_class> > BuildPolytope::getFacets() const {return facets;}


string BuildPolytope::getLatteVRepFile() const {return fileBaseName + ".vrep.latte";}
string BuildPolytope::getLatteVRepDualFile() const {return fileBaseName + ".vrep.dual.latte";}
string BuildPolytope::getLatteHRepFile() const { return fileBaseName + ".hrep.latte";}
string BuildPolytope::getLatteHRepDualFile() const { return fileBaseName + ".hrep.dual.latte";}
string BuildPolytope::getPolymakeFile() const { return fileBaseName + ".polymake";}
string BuildPolytope::getPolymakeDualFile() const { return getDualFileBaseName() +".polymake";}
string BuildPolytope::getDualFileBaseName() const { return fileBaseName + ".dual";} //private function.


/**
 * Returns the homogenize vertices
 */
vector<vector<mpq_class> > BuildPolytope::getVertices()
{
	findVertices();
	return points;
}//getVertices

int BuildPolytope::getVertexCount()
{
	findVertices();
	return points.size();
}

int BuildPolytope::getVertexDualCount()
{
	findVerticesDual();
	return dualVertices.size();
}

/**
 * Right now the dualVertices are in the form [b v1 ... vn], after this call
 * they will be in form [1, v1/b, ..., vn/b].
 */
void BuildPolytope::homogenizeDualVertices()
{
	for(size_t i = 0; i < dualVertices.size(); ++i)
	{
		assert(dualVertices[i][0] > 0);

		for(size_t j = 1; j < dualVertices[i].size(); ++j)
			dualVertices[i][j] /= dualVertices[i][0];
		dualVertices[i][0] = 1;
	}//for i;
}//homogenizeDualVertices

bool BuildPolytope::isCentered()
{
	buildPolymakeFile();

	//ask polymake
	system((string("polymake ") + getPolymakeFile() + " CENTERED > /dev/null ").c_str());

	//next, read the polymake file.
	ifstream file;
	string line;
	file.open(getPolymakeFile().c_str());

	for (getline(file, line, '\n'); line != "CENTERED"; getline(file, line, '\n'))
		;

	char ans;
	ans = file.get();
	file.close();

	return (ans == '1');
}//isCentered


/**
 * true if s is a number in for form 4, 4.5, or 4/5

bool BuildPolytope::isStringNumber(const string & s) const
{
	bool signCount =0, decimalCount=0, rationalCount=0;

	if ( s.length() == 0)
		return false;

	for(int i = 0; i < s.length(); ++i)
	{
		if ( isdigit(s[i]))
			continue;

		if ( s[i] == '-')
		{
			if (signCount)
				return false; //too many signs
			else
				signCount++;
		}
		else if (s[i] == '.')
		{
			if (decimalCount)
				return false; //too many decimal points
			else
				decimalCount++;
		}
		else if ( s[i] == '/')
		{
			if (rationalCount)
				return false; //too many /
			else
				rationalCount++;
		}
		else
			return false;

	}//for i
	return true;
}//isStringNumber
 */

bool BuildPolytope::isSimplicial()
{
	buildPolymakeFile();

	//ask polymake
	system((string("polymake ") + getPolymakeFile() + " SIMPLICIAL > /dev/null ").c_str());


	//next, read the polymake file.
	ifstream file;
	string line;
	file.open(getPolymakeFile().c_str());

	for (getline(file, line, '\n'); line != "SIMPLICIAL"; getline(file, line, '\n'))
		;

	char ans;
	ans = file.get();
	file.close();
	
	return (ans == '1');
}

bool BuildPolytope::isSimple()
{
	buildPolymakeFile();

	//ask polymake
	system((string("polymake ") + getPolymakeFile() + " DIM > /dev/null ").c_str());
	system((string("polymake ") + getPolymakeFile() + " SIMPLE > /dev/null ").c_str());


	//next, read the polymake file.
	ifstream file;
	string line;
	file.open(getPolymakeFile().c_str());

	for (getline(file, line, '\n'); line != "SIMPLE"; getline(file, line, '\n'))
		;

	char ans;
	ans = file.get();
	file.close();

	return (ans == '1');
}


/**
 * If a polytope is simplicial, then its dual is simple.
 */
bool BuildPolytope::isDualSimplicial()
{
	return ! isSimple();
}//isDualSimplicial().

/**
 * If a polytope is simple, then its dual is simplicial.
 */
bool BuildPolytope::isDualSimple()
{
	return ! isSimplicial();
}//isDualSimplicial().


/**
 * @list: a vector of facets or vertices. Will mult each row by a different number to clear each denominators for latte
 *
 */
void BuildPolytope::makeIntegerRows(vector<vector<mpq_class> > &list)
{

	for (int i = 0; i < (int) list.size(); ++i)
	{
		mpz_class currentLCM(1);  //sorry, I was having a hard time calling with mpz_class
								  //it seems to only work with mpz_t, hence these two variables. (currentLCM and currentlcm)


		//find the lcm of al the denominator for this row only.
		for (int k = 0; k < ambientDim + 1; ++k)
		{
			mpz_t currentlcm; //set the output variable to 1.
			mpz_init_set_si(currentlcm, 1);

			if ( list[i][k] == mpz_class(0))
				continue;

			mpz_lcm(currentlcm, currentLCM.get_mpz_t(), list[i][k].get_den_mpz_t());
			currentLCM = mpz_class(currentlcm);//save the current lcm in the c++ object.
		}//for k

		assert(currentLCM > 0);

		//mult each element in this row by the lcm.
		if (currentLCM != mpz_class(1))
		{
			for (int k = 0; k < ambientDim + 1; ++k)
			{
				list[i][k] = list[i][k] * currentLCM;
				assert(list[i][k].get_den() == mpz_class(1));
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

}//makeIntegerRows();

/**
 * @list: a vector of facets or vertices. Will mult each row by the same number to clear each denominators for latte
 *
 */
void BuildPolytope::makeIntegerList(vector<vector<mpq_class> > &list)
{
	mpz_t currentLCM;
	mpz_init_set_si(currentLCM, 1);

	//currentLCM = 1
	//loop over everything and find the lcm of the denom.
	for (int i = 0; i < (int) list.size(); ++i)
	{
		for (int k = 0; k < (int) list[i].size(); ++k)
		{
			if (list[i][k] != mpq_class(0))
			{
				mpz_lcm(currentLCM, currentLCM,
						list[i][k].get_den_mpz_t());
			}//if not zero.
		}//for k
	}//for i

	//cout << "makeIntegerList::" << endl;
	//debugPrintList(list);
	//cout << "makeIntegerList:: the lcm is " << mpz_class(currentLCM) << endl;

	assert(currentLCM > 0);
	
	//loop over everything again and times by the lcm.
	for (size_t i = 0; i < list.size(); ++i)
	{
		for (size_t k = 0; k < list[i].size(); ++k)
		{
			list[i][k] = list[i][k]* mpz_class(currentLCM);
		}//for k
	}//for i
}//makeIntegerList


void BuildPolytope::setBaseFileName(const string & n) {fileBaseName = n;}

/**
 * Children of this class can use this method when building polytopes.
 */
void BuildPolytope::setIntegerPoints(bool t) { integerPoints = t; }
void BuildPolytope::setBuildPolymakeFile(bool t) { createdPolymakeFile = t;}
void BuildPolytope::setBuildLatteVRepDualFile(bool t) { createdLatteVRepDualFile = t;}

/**
 * This method makes a random polytope.
 */
void BuildPolytope::forDebugging()
{
	ambientDim = 3;
	srand(time(NULL));
	
	for (int i = 0; i < ambientDim+1+5; ++i)
	{
		vector<mpq_class> onePoint;
		onePoint.push_back(1);
		for(int j = 0; j < ambientDim; ++j)
			if ( integerPoints )
				onePoint.push_back(mpq_class(rand()%100,1));
			else
				onePoint.push_back(mpq_class(rand()%100, rand()%25));
		points.push_back(onePoint);
	}
}

/**
 * Prints a facet, vertex, dual facet, or dual vertex list.
 */
void BuildPolytope::debugPrintList(const vector<vector<mpq_class> > &list)
{
	for(int i = 0; i < (int) list.size(); ++i)
	{
		cout << "i " << i << "= ";
		for(int j = 0; j < (int) list[i].size(); ++j)
			cout << list[i][j] << " ";
		cout << endl;
	}//for i
}//debugPrintList




