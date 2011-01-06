/*
 *
 *
 *  Created on: June 15, 2010
 *      Author: bedutra
 *
 *  Base class for making polytopes.
 *
 *  h-reps: will scale the coeff. to integers
 *  v-reps: no scaling is performed.
 *
 *  h-reps dual: not implemented...I never needed this yet!
 *  v-reps dual: finds the scaled facet equations (so the dual is dilated to integer points).
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

    vector< vector<mpq_class> > points; //each point lives in ambientDim space. Assumes the points do not have a leading 1.
    vector< vector<mpq_class> > facets; //facets of the polytope found by polymake
    int numAffineHull;		//number of affine hull facets found by polymake (saved at the end of the facets vector)

    string getDualFileBaseName() const;
public:
    BuildPolytope();
	
	//A B C D E F G H I J K L M N O P Q R S T U V W X Y Z

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
	void buildLatteVRepFile();//build the v-rep file in latte style.
	void buildLatteVRepDualFile(); //build the v-rep for the dual (dilated) polytope


	void findDimentions(); //finds the dim of the polytope form polymake.
	void findFacets(bool dilateFacets = false);     //saves the facets in our vector.
	void findVertices();   //saves the vertices in points (overrides the org. points vector).
	void findVerticesDual();//finds the facet equations.
	

	int getAmbientDim() const; //returns abmient dim.
	int getDim() const; //returns the dim. Assumes polymake has been called.
	
	string getLatteVRepFile() const; //Return file names.
	string getLatteVRepDualFile() const;
	string getLatteHRepFile() const;
	string getLatteHRepDualFile() const;
	string getPolymakeFile() const;
	string getPolymakeDualFile() const;

	int getVertexCount();
	int getVertexDualCount();

	bool isSimplicial();
	bool isDualSimplicial(); //makes a temp file. I could not find a polymake command that would give this w/o making a temp file.
	
	void convertFacetEquations();//mult. each equations by a (different) number to clear the denominators for latte.
	void dilateFacetEquations(); //mult. each equations by the same number  to clear the denominators for latte. This is important for V-reps!
	void setBaseFileName(const string & n); //sets the file name root.
	void setIntegerPoints(bool t); //should the polytope be interger?
	
	void forDebugging();
};//BuildRandomPolytope



#endif /* BUILD_POLYTOPE_H_ */

