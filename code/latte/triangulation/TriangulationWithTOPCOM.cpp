/* TriangulationWithTOPCOM.cpp -- Compute a triangulation using TOPCOM

   Copyright 2006 Matthias Koeppe
   Derived from the TOPCOM source files points2placingtriang.cc, ComputeTriangs.cc
     which are copyright 1999 Joerg Rambau  

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

#include <list>
#include <iostream>
#include "config.h"
#include "triangulation/TriangulationWithTOPCOM.h"
#include "latte_gmp.h"
#include "latte_system.h"
#ifdef HAVE_TOPCOM_LIB
#include <Matrix.hh>
#include <PointConfiguration.hh>
#include <PlacingTriang.hh>
#include <Symmetry.hh>
#include <SymmetricBFS.hh>
#endif

using namespace std;

static vector<listVector *>
ray_array(listCone *cone)
{
  int num_rays = lengthListVector(cone->rays);
  vector<listVector *> rays(num_rays);
  int j;
  listVector *ray;
  for (j = 0, ray = cone->rays; ray!=NULL; j++, ray = ray->rest)
    rays[j] = ray;
  return rays;
}


//
// Functions using the TOPCOM library.
//
#ifdef HAVE_TOPCOM_LIB
static listCone *
triangulation_to_cones(const SimplicialComplex &sc,
		       listCone *cone)
{
  // Copy rays into an array, so we can index them.
  vector<listVector *> rays = ray_array(cone);
  // Iterate through simplices and collect them.
  listCone *cones = NULL;
  for (SimplicialComplex::const_iterator i = sc.begin(); i!=sc.end(); ++i) {
    const Simplex &simplex = *i;
    listCone *c = createListCone();
    c->vertex = new Vertex(*cone->vertex);
    for (Simplex::const_iterator j = simplex.begin(); j!=simplex.end(); ++j) {
      size_t index = *j;
      c->rays = new listVector(rays[index]->first, c->rays);
    }
    c->rest = cones;
    cones = c;
  }
  return cones;
}

static PointConfiguration
cone_to_point_configuration(listCone *cone, int numOfVars)
{
  int num_rays = lengthListVector(cone->rays);
  Matrix matrix(numOfVars, num_rays);
  size_type i, j;
  listVector *ray;
  for (j = 0, ray = cone->rays; ray!=NULL; j++, ray = ray->rest) {
    for (i = 0; i<numOfVars; i++) {
      mpz_class mpz = convert_ZZ_to_mpz(ray->first[i]);
      matrix(i, j) = Rational(mpz.get_mpz_t());
    }
  }
  PointConfiguration points(matrix);
  return points;
}

listCone *
triangulate_cone_with_TOPCOM(listCone *cone, int numOfVars)
{
  PointConfiguration points = cone_to_point_configuration(cone, numOfVars);
  if (points.rank() > points.no()) {
    std::cerr << "rank must not be larger than no of points." << std::endl;
    exit(1);
  }
  Chirotope chiro(points, /*preprocess:*/ false);
  PlacingTriang triang(chiro);
  std::cout << "Triangulation: " << triang << std::endl;
  return triangulation_to_cones(triang, cone);
}
#endif


//
// Functions using the TOPCOM binaries
//

#ifdef HAVE_TOPCOM_BIN

static listCone *
read_TOPCOM_triangulation(istream &in,
			  listCone *cone)
{
  // Copy rays into an array, so we can index them.
  vector<listVector *> rays = ray_array(cone);
  // Read a triangulation in the TOPCOM format:
  // T[1]:={{0,1,2,3,4,5,6,8},{0,1,2,3,4,5,7,8}};
  char c;
  // Find and consume left brace of triangulation.
  do { in.get(c); } while (in.good() && c != '{');
  listCone *cones = NULL;
  do {
    // Consume left brace of simplex.
    in.get(c);
    if (!(in.good() && c == '{')) break;
    listCone *simp = createListCone();
    simp->vertex = new Vertex(*cone->vertex);
    while (in.good() && in.peek()!='}') {
      int index;
      in >> index;
      simp->rays = new listVector(rays[index]->first, simp->rays);
      // Skip the ',' of simplex.
      if (in.peek()==',') in.get(c);
    }
    simp->rest = cones;
    cones = simp;
    if (!in.good()) break;
    // Consume the '}' of simplex.
    in.get(c);
    // Expect ',' or '}' of triangulation.
    in.get(c);
    if (!in.good()) break;
    if (c == '}') return cones;
  } while (c == ',');
  // Error situation:
  freeListCone(cones);
  return NULL;
}

static void
write_TOPCOM_point_configuration(ostream &out,
				 listCone *cone, int numOfVars)
{
  out << "[";
  listVector *ray;
  for (ray = cone->rays; ray!=NULL; ray=ray->rest) {
    out << "[";
    int i;
    out << ray->first[0];
    for (i = 1; i<numOfVars; i++)
      out << ", " << ray->first[i];
    out << "]";
    if (ray->rest != NULL)
      out << "," << endl;
  }
  out << "]" << endl;
}

#ifndef HAVE_TOPCOM_LIB
// Alternative implementation if the TOPCOM library cannot be used.
listCone *
triangulate_cone_with_TOPCOM(listCone *cone, int numOfVars)
{
  string topcom_input_file_name = temporary_file_name("topcom_input");
  string topcom_output_file_name = temporary_file_name("topcom_output");
  {
    ofstream topcom_input(topcom_input_file_name.c_str());
    write_TOPCOM_point_configuration(topcom_input, cone, numOfVars);
  }
  system_with_error_check(TOPCOM_POINTS2PLACINGTRIANG " < " + topcom_input_file_name
			  + " > " + topcom_output_file_name);
  ifstream topcom_output(topcom_output_file_name.c_str());
  listCone *triangulation = read_TOPCOM_triangulation(topcom_output, cone);
  if (!triangulation) {
    cerr << "Failed to read triangulation from TOPCOM output" << endl;
    exit(1);
  }
  return triangulation;
}
#endif

list<listCone *>
all_triangulations_of_cone_with_TOPCOM(listCone *cone, int numOfVars)
{
  string topcom_input_file_name = temporary_file_name("topcom_input");
  string topcom_output_file_name = temporary_file_name("topcom_output");
  {
    ofstream topcom_input(topcom_input_file_name.c_str());
    write_TOPCOM_point_configuration(topcom_input, cone, numOfVars);
  }
  system_with_error_check(TOPCOM_POINTS2TRIANGS " < " + topcom_input_file_name
			  + " > " + topcom_output_file_name);
  ifstream topcom_output(topcom_output_file_name.c_str());
  list<listCone *> triangulations;
  while (topcom_output.good()) {
    listCone *triangulation = read_TOPCOM_triangulation(topcom_output, cone);
    if (!triangulation) break;
    triangulations.push_front(triangulation);
  }
  return triangulations;
}

#endif
