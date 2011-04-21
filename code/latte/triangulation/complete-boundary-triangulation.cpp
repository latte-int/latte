/* complete-boundary-triangulation.cpp -- 
	       
   Copyright 2007 Matthias Koeppe
   Copyright 2011 Christof Soeger

   This file is part of LattE.
   
   LattE is free software; you can redistribute it and/or modify it
   under the terms of the version 2 of the GNU General Public License
   as published by the Free Software Foundation.

   LattE is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with LattE; if not, write to the Free Software Foundation,
   Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
*/

#include <cstdlib>
#include <iostream>
#include <string>

#include "barvinok/barvinok.h"
#include "print.h"
#include "ReadPolyhedron.h"
#include "ReadSubcones.h"
#include "triangulation/BoundaryTriangulation.h"

using namespace std;

BarvinokParameters parameters;
ReadPolyhedronData read_polyhedron_data;

string output_filename;
string output_completion_data_filename;
string vector_input_filename;

int main(int argc, char *argv[])
{
  int i;
  for (i = 1; i<argc; i++) {
    if (read_polyhedron_data.parse_option(argv[i])) {}
    else if (strncmp(argv[i], "--output-completion-data=", 25) == 0) {
      output_completion_data_filename = argv[i] + 25;
    }
    else if (strncmp(argv[i], "--output-cones=", 15) == 0) {
      output_filename = argv[i] + 15;
    }
    else if (strncmp(argv[i], "--input-completion-vector=", 26) == 0) {
      vector_input_filename = argv[i] + 26;
    }
    else {
      cerr << "Unknown argument: " << argv[i] << endl;
      exit(1);
    }
  }
  if (output_filename.size() == 0 && output_completion_data_filename.size() == 0) {
    cerr << "The option --output-cones=FILENAME or --output-completion-data=FILENAME must be given." << endl;
    exit(1);
  }
  Polyhedron *Poly = read_polyhedron_data.read_polyhedron(&parameters);
  listCone *cone = Poly->cones;
  cerr << lengthListCone(cone) << " cones." << endl;
  ConeConsumer *consumer = new PrintingConeConsumer(output_filename);

  parameters.Number_of_Variables
    = cone->rays->first.length();
  
  if (output_completion_data_filename.size() != 0) { // generate data which helps to find a vector to complete triangulation
    mat_ZZ alpha; //determinantes
    mat_ZZ basis; //LLL-Basis
    mat_ZZ F;     //avoid zeros
    generate_interior_vector_simplified_data
      (Poly->cones, parameters.Number_of_Variables, alpha, basis, F);
    //output
    ofstream out_alpha((output_completion_data_filename+".alpha").c_str());
	out_alpha << alpha;
	out_alpha.close();
    ofstream out_basis((output_completion_data_filename+".basis").c_str());
	out_basis << basis;
	out_basis.close();
    ofstream out_F(output_completion_data_filename.append(".F").c_str());
	out_F << F;
	out_F.close();
  }
  else if (vector_input_filename.size() == 0) { //try to find a vector to complete triangulation
    complete_boundary_triangulation_of_cone_with_subspace_avoiding_facets
      (Poly->cones, &parameters, *consumer);
  }
  else {    //use the vector from file
    ifstream cv_in(vector_input_filename.c_str());

//  listVector *cvl = readListVector(cvin);
//	 vec_ZZ completion_vector = cvl->first;
	 vec_ZZ completion_vector;
	 cv_in >> completion_vector;
	 if(!cv_in.good()) {
      cerr << "Problems reading completion_vector!"<<endl;
      exit(1);
    }
    complete_boundary_triangulation_of_cone_with_subspace_avoiding_facets
      (Poly->cones, &parameters, completion_vector, *consumer);
  }
  // this consumes Poly->cones, so:
  Poly->cones = NULL;
  
  delete Poly;
  delete consumer;
  if (output_filename.size() != 0) {
  cerr << "Wrote complete boundary triangulation to file `"
       << output_filename << "'" << endl;
  }
  return 0;
};


