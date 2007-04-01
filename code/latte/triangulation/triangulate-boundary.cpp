/* triangulate-boundary.cpp -- Triangulate the boundary of a cone
	       
   Copyright 2007 Matthias Koeppe

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

bool output_subcones;
string output_filename;

int main(int argc, char *argv[])
{
  int i;
  output_subcones = false;
  for (i = 1; i<argc; i++) {
    if (read_polyhedron_data.parse_option(argv[i])) {}
    else if (strncmp(argv[i], "--output-cones=", 15) == 0) {
      output_filename = argv[i] + 15;
    }
    else if (strncmp(argv[i], "--output-subcones=", 18) == 0) {
      output_filename = argv[i] + 18;
      output_subcones = true;
    }
    else {
      cerr << "Unknown argument: " << argv[i] << endl;
      exit(1);
    }
  }
  if (output_filename.size() == 0) {
    cerr << "One of the options --output-cones=FILENAME and --output-subcones=FILENAME must be given." << endl;
    exit(1);
  }
  Polyhedron *Poly = read_polyhedron_data.read_polyhedron(&parameters);
  listCone *cone = Poly->cones;
  cerr << lengthListCone(cone) << " cones." << endl;
  ConeConsumer *consumer;
  if (output_subcones)
    consumer = new SubconePrintingConeConsumer(copyCone(cone), output_filename);
  else
    consumer = new PrintingConeConsumer(output_filename);

  parameters.Number_of_Variables
    = cone->rays->first.length();
  
  for (cone = Poly->cones; cone!=NULL; cone=cone->rest)
    compute_triangulation_of_boundary(cone, &parameters, *consumer);
  
  delete Poly;
  delete consumer;
  cerr << "Wrote triangulation of boundary to file `"
       << output_filename << "'" << endl;
  return 0;
};


