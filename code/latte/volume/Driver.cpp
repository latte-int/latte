/*
 * Driver.cpp
 *
 *  Created on: Jun 24, 2010
 *      Author: bedutra
 */



#include <cstdlib>
#include <iostream>
#include <string>

#include "barvinok/barvinok.h"
#include "ReadPolyhedron.h"
#include "triangulation/triangulate.h"
#include "convert.h"
#include "print.h"

#include "gnulib/progname.h"

using namespace std;

BarvinokParameters parameters;
ReadPolyhedronData read_polyhedron_data;
string output_filename;

int main(int argc, char *argv[])
{
  set_program_name(argv[0]);

  int i;
  for (i = 1; i<argc; i++) {
    if (read_polyhedron_data.parse_option(argv[i])) {}
    else if (strncmp(argv[i], "--output-cones=", 15) == 0) {
      output_filename = argv[i] + 15;
    }
    else {
      cerr << "Unknown argument: " << argv[i] << endl;
      exit(1);
    }
  }
  Polyhedron *poly = read_polyhedron_data.read_polyhedron(&parameters);

parameters.Number_of_Variables  = poly->numOfVars;

  printListConeToFile(output_filename.c_str(), poly->cones, poly->numOfVars);

listCone * triangulatedCones;
ZZ valutation = ZZ();
	for(listCone *currentCone = poly->cones; currentCone; currentCone = currentCone->rest)
	{
		triangulatedCones = triangulateCone(currentCone, poly->numOfVars, &parameters);

		for( listCone * simplicalCone = triangulatedCones; simplicalCone; simplicalCone = simplicalCone->rest)
		{

			ZZ det = determinant((createConeDecMatrix(simplicalCone, poly->numOfVars, poly->numOfVars)));
			det = abs(det);
			det *= simplicalCone->coefficient;
			cout << "this det=" << det << endl;
			valutation += det;

		}//for every simply cone.
	}//for every cone

	cout << "TOTAL: " << valutation << endl;

  printListConeToFile((output_filename + ".triangulated").c_str(), triangulatedCones, 4); //poly->numOfVars);

  return 0;
}




