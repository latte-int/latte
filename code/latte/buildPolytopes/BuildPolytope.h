/*
 *
 *
 *  Created on: June 15, 2010
 *
 *
 *  Base class for making polytopes.
 *
 *  h-reps: will scale the coeff. to integers
 *  v-reps: no scaling is performed.
 *
 *  h-reps dual: scales the vertices to integers.
 *  v-reps dual: scales the facets to integers (be finding the lcm of every denominator)
 */

#ifndef BUILD_POLYTOPE_H_
#define BUILD_POLYTOPE_H_

#include <iomanip>
#include <string>
#include <fstream>
#include <ctime>
#include <vector>
#include <sstream>
#include <cstdlib>
#include "gmp.h"
#include <gmpxx.h>



using namespace std;

class BuildPolytope
{
protected:
    int ambientDim; 		//dim of space the poly lives in.
    int dim;				//dim of the polytope found by polymake
    bool integerPoints;		//if true, makes random integer points; if false, makes random rational points.
    
    string fileBaseName; 		//base file name for polymake and latte.
    							//fileBaseName.polymake
    							//fileBaseName.vrep.latte
    							//fileBaseName.hrep.latte    	
    bool createdPolymakeFile;	//ture if these files have been created with this object.
    bool createdPolymakeDualFile;
    bool createdLatteVRepFile;
    bool createdLatteHRepFile;

    bool createdLatteHRepDualFile; //save as above, only for a dual polytope
    bool createdLatteVRepDualFile;


    vector< vector<mpq_class> > facets; //facets of the polytope found by polymake
    vector< vector<mpq_class> > dualVertices;//assumes points DO have a leading "1"..so the length is (dim+1)
    vector< vector<mpq_class> > dualFacets;//could be different from the vertices if the org. polytope is not centered.
    int numAffineHull;		//number of affine hull facets found by polymake (saved at the end of the facets vector)

    string getDualFileBaseName() const;
private:
    vector< vector<mpq_class> > points; //each point lives in ambientDim space. Assumes the points DO have a leading 1.
public:
    BuildPolytope();
	
	//A B C D E F G H I J K L M N O P Q R S T U V W X Y Z

    void addPoint(vector<mpq_class> onePoint); //adds [1 onePoint] ot the point vector
    void centerPolytope();
    void clearPoints();

    //Delete the generated files.
	void deletePolymakeFile();
	void deletePolymakeDualFile();
	void deleteLatteVRepFile();
	void deleteLatteVRepDualFile();
	void deleteLatteHRepFile();
	void deleteLatteHRepDualFile();


	void buildPolymakeFile(); //builds the polymake file.
	void buildPolymakeDualFile(); //build the dual polymake file.
	void buildLatteHRepFile();//builds the h-rep file in latte style
	void buildLatteHRepDualFile(); //builds h-rep of the dual.
	void buildLatteVRepFile();//build the v-rep file in latte style.
	void buildLatteVRepDualFile(); //build the v-rep for the dual (dilated) polytope


	void findDimentions(); //finds the dim of the polytope form polymake.
	void findFacets();     //saves the facets in our vector.
	void findFacetsDual();//find the dual facets.
	void findAffineHull();//find the rests of the facets.
	void findVertices();   //saves the vertices in points (overrides the org. points vector).
	void findVerticesDual();//finds the facet equations.
	

	int getAmbientDim() const; //returns abmient dim.
	int getDim() const; //returns the dim. Assumes polymake has been called.
	
	vector<vector<mpq_class> > getFacets() const;
	string getLatteVRepFile() const; //Return file names.
	string getLatteVRepDualFile() const;
	string getLatteHRepFile() const;
	string getLatteHRepDualFile() const;
	string getPolymakeFile() const;
	string getPolymakeDualFile() const;

	vector<vector<mpq_class> > getVertices();
	int getVertexCount();
	int getVertexDualCount();

	void homogenizeDualVertices(); //divides so the first slot is 1.

	bool isCentered();
	bool isSimplicial();
	bool isSimple();
	bool isDualSimplicial();
	bool isDualSimple();
	
	void makeIntegerRows(vector<vector<mpq_class> > &list);//mult. each equations by a (different) number to clear the denominators for latte.
	void makeIntegerList(vector<vector<mpq_class> > &list);//mult. each equations by a number to clear the denominators for latte.
	void setBaseFileName(const string & n); //sets the file name root.
	void setIntegerPoints(bool t); //should the polytope be interger?
	void setBuildPolymakeFile(bool t); //used for finding the dual vertices.
	
	void forDebugging();
	void debugPrintList(const vector<vector<mpq_class> > &list);
};//BuildRandomPolytope



#endif /* BUILD_POLYTOPE_H_ */

