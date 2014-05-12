/*
 * nonlinearOptimizationDriver.cpp
 *
 *  Created on: Apr 23, 2014
 *      Author: bedutra
 */

#include <iostream>
#include <fstream>
#include <string>

#include "banner.h"
#include "barvinok/barvinok.h"
#include "ReadPolyhedron.h"
#include "LattException.h"
#include "nonlinearOptimization/WeightedExponentialSubs.h"
#include "nonlinearOptimization/BoxOptimization.h"
#include "print.h"
#include "dual.h"
#include "rational.h"
#include "integration/burstTrie.h"
#include "integration/PolyTrie.h"

using namespace std;

int main2(int argc, const char *argv[]) ;
int main1(int argc, const char *argv[]) ;

int main1(int argc, const char *argv[]) {
	if (argv[1][0] == '1')
		main1(argc - 1, argv + 1);
	else
		main2(argc -1, argv + 1);

	return 0;
}


int main(int argc, const char *argv[]) {

	string linFormFileName, boxFileName;



	latte_banner(cerr);

	cerr << "Invocation: ";
	for (int i = 0; i < argc; i++) {
		cerr << argv[i] << " ";
	}
	cerr << endl;

	for (int i = 1; i < argc; i++) {
		if (strncmp(argv[i], "--linFile=", 10) == 0){
			linFormFileName = string(argv[i] + 10);
		} else if (strncmp(argv[i], "--boxFile=", 10) == 0) {
			boxFileName = string(argv[i] + 10);
		} else {
			cerr << "Unknown command/option " << argv[i] << endl;
			THROW_LATTE_MSG(LattException::ue_BadCommandLineOption, argv[i]);
		}
	} //for i.

	if (linFormFileName.length() == 0 || boxFileName.length() == 0) {
		cerr << "Files are missing" << endl;
		THROW_LATTE(LattException::ue_FileNameMissing);
	}

	ifstream linFormFile(linFormFileName.c_str());
	string linFormStr;
	getline(linFormFile, linFormStr);
	cout << "lin form str: " << linFormStr.c_str() << endl;
	linFormSum originalLinearForm;
	loadLinForms(originalLinearForm, linFormStr.c_str());
	linFormFile.close();

	ifstream boxFile(boxFileName.c_str());
	int dim;
	boxFile >> dim;

	vec_ZZ lowerBound, upperBound;
	lowerBound.SetLength(dim);
	upperBound.SetLength(dim);
	for(int i = 0; i < dim; ++i)
		boxFile >> lowerBound[i] >> upperBound[i];

	cout << "lb: " << lowerBound << endl;
	cout << "ub: " << upperBound << endl;


	mpq_class weightedCount = computeWeightedCountingBox(lowerBound, upperBound, originalLinearForm);
	cout << "Final count: " << weightedCount << endl;



	return 0;
}


int main2(int argc, const char *argv[]) {
	ReadPolyhedronData read_polyhedron_data;
	string linFormFileName;

	struct BarvinokParameters *params = new BarvinokParameters;

	latte_banner(cerr);

	cerr << "Invocation: ";
	for (int i = 0; i < argc; i++) {
		cerr << argv[i] << " ";
	}
	cerr << endl;

	params->substitution = BarvinokParameters::PolynomialSubstitution;
	params->decomposition = BarvinokParameters::DualDecomposition;
	params->max_determinant = 1;
	for (int i = 1; i < argc; i++) {
		if (read_polyhedron_data.parse_option(argv[i])) {
		} 
		else if (strncmp(argv[i], "--linFile=", 10) == 0){
			linFormFileName = string(argv[i] + 10);
		} else {
			cerr << "Unknown command/option " << argv[i] << endl;
			THROW_LATTE_MSG(LattException::ue_BadCommandLineOption, argv[i]);
		}
	} //for i.

	if (read_polyhedron_data.expect_filename) {
		cerr << "Filename missing" << endl;
		THROW_LATTE(LattException::ue_FileNameMissing);
	}

	const char *fileName = read_polyhedron_data.filename.c_str();

	Polyhedron *Poly = read_polyhedron_data.read_polyhedron(params);

	cout << "fix me**********************************************" << endl;
	params->Number_of_Variables = Poly->numOfVars;
	
	if (Poly->cones != NULL && Poly->cones->rays == NULL) {

			// Only facets computed, for instance by using the 4ti2
			// method of computing vertex cones.  So dualize twice to
			// compute the rays.
			cerr << "(First computing their rays... ";
			cerr.flush();
			dualizeCones(Poly->cones, Poly->numOfVars, params);
			dualizeCones(Poly->cones, Poly->numOfVars, params); // just swaps
			cerr << "done; rays are now computed) \n";
			cerr.flush();

	}
	//printListCone(Poly->cones, Poly->numOfVars);
	//cout << "*****************" << endl;
			
	//cout << "numrays=" << lengthListVector(Poly->cones->rays) << endl;
	//cout << "num vars" << params->Number_of_Variables << endl;
	//cout << Poly->numOfVars;
	
	assert(Poly->cones->rays != NULL);

	//Weighted_Exponential_Single_Cone_Parameters *exp_param =
	//		new Weighted_Exponential_Single_Cone_Parameters(*params);
	//delete params;
	//params = exp_param;

	if ( linFormFileName.length() == 0) {
		cerr << "lin fomr file missing: " << endl;
		THROW_LATTE(LattException::ue_FileNameMissing);
	}
	
	ifstream linFormFile(linFormFileName.c_str());
	string linFormStr;
	getline(linFormFile, linFormStr);
	cout << "lin form str: " << linFormStr.c_str() << endl;
	linFormSum originalLinearForm;
	loadLinForms(originalLinearForm, linFormStr.c_str());
		

	cout << "what should happen next in the general setting?" << endl;
	mpq_class weighted_count = computeWeightedExponentialResidue(Poly->cones, params, originalLinearForm);
	cout << "Final count: " << weighted_count << endl;
	destroyLinForms(originalLinearForm);

	return 0;
}
