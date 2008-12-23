/* polyhedron-to-cones.cpp -- Compute the vertex cones or the homogenization of a polyhedron
	       
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
#include "ReadPolyhedron.h"
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
  if (output_filename.size() == 0) {
    cerr << "The option --output-cones=FILENAME is mandatory." << endl;
    exit(1);
  }
  printListConeToFile(output_filename.c_str(), poly->cones, poly->numOfVars);
  cerr << "Cones (homogenized: " << boolalpha << poly->homogenized
       << ", dualized: " << poly->dualized
       << ") printed to file `" << output_filename << "'." << endl;
  return 0;
}
