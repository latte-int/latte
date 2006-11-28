/* test-triangulation.cpp -- Stand-alone program for triangulation
	       
   Copyright 2006 Matthias Koeppe

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

#include <iostream>
#include <string>
#include <cctype>
#include "print.h"
#include "triangulate.h"
#include "dual.h"

using namespace std;

// Read a cone in the format of `printListCone'.
// However, this is NOT a general function; we only read the extreme
// rays. 
listCone *
read_cone(istream &in)
{
  string s;
  while (in.good()) {
    in >> s;
    if (s == "rays:") break;
  }
  if (!in.good()) {
    return NULL;
  }
  listCone *cone = createListCone();
  while (in.good()) {
    vec_ZZ v;
    while (isspace(in.peek())) {
      char c;
      in.get(c);
    }
    if (in.peek() != '[') break;
    in >> v;
    if (in.good()) {
      cone->rays = appendVectorToListVector(v, cone->rays);
    }
  }
  return cone;
}

int main(int argc, char **argv)
{
  if (argc != 2) {
    cerr << "usage: triangulate LATTE-CONE-FILE" << endl;
    exit(1);
  }
  ifstream in(argv[1]);
  if (in.bad()) {
    cerr << "triangulate: Unable to open " << argv[1] << endl;
    exit(1);
  }
  listCone *cone = read_cone(in);
  if (!cone) {
    cerr << "triangulate: Parse error in file." << endl;
    exit(1);
  }

  if (cone->rays == NULL) {
    cerr << "triangulate: No rays." << endl;
    exit(2);
  }
  BarvinokParameters params;
  params.Number_of_Variables = cone->rays->first.length();
  cone->vertex = new Vertex(new rationalVector(params.Number_of_Variables));

  // Compute facets.
  cone = dualizeCones(cone, params.Number_of_Variables);
  cone = dualizeBackCones(cone, params.Number_of_Variables);
  
  cout << "*** Input cone:" << endl;
  printListCone(cone, params.Number_of_Variables);
  
  //   params.triangulation
  //     = BarvinokParameters::RegularTriangulationWithCdd;
  params.triangulation = BarvinokParameters::PlacingTriangulationWithTOPCOM;
  listCone *triang
    = triangulateCone(cone, params.Number_of_Variables, &params);
  listCone *t;
  for (t = triang; t!=NULL; t = t->rest) {
    computeDetAndFacetsOfSimplicialCone(t, params.Number_of_Variables);
  }
  cout << "*** Triangulation:" << endl;
  printListCone(triang, params.Number_of_Variables);
  return 0;
}
