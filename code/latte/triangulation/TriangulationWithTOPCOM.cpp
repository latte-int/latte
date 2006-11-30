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

#include "triangulation/TriangulationWithTOPCOM.h"
#include "latte_gmp.h"
#include <Matrix.hh>
#include <PointConfiguration.hh>
#include <PlacingTriang.hh>
#include <Symmetry.hh>
#include <SymmetricBFS.hh>
#include <list>

static listCone *
triangulation_to_cones(const SimplicialComplex &sc,
		       listCone *cone)
{
  // Copy rays into an array, so we can index them.
  int num_rays = lengthListVector(cone->rays);
  vector<listVector *> rays(num_rays);
  size_type j;
  listVector *ray;
  for (j = 0, ray = cone->rays; ray!=NULL; j++, ray = ray->rest)
    rays[j] = ray;
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

list<listCone *>
all_triangulations_of_cone_with_TOPCOM(listCone *cone, int numOfVars)
{
  PointConfiguration points = cone_to_point_configuration(cone, numOfVars);
  if ((points.no() < 2) || (points.rank() < 2)) {
    std::cerr << "no of points and rank must be at least two." << std::endl;
    exit(1);
  }
  if (points.rank() > points.no()) {
    std::cerr << "rank must not be larger than no of points." << std::endl;
    exit(1);
  }
  Chirotope chiro(points, /*preprocess:*/ false);
  size_type no(chiro.no());
  size_type rank(chiro.rank());
  SymmetryGroup symmetries(no);
  SimplicialComplex seed = PlacingTriang(chiro);
  const SymmetryGroup seed_symmetries(symmetries, seed);
  
  SymmetricBFS bfs(no, rank, points, chiro, 
		   symmetries, seed, seed_symmetries, 
		   /*output_triangs:*/ true, /*fine_only:*/ false);
  
  return list<listCone *>();
}
