/* TriangulationWithTOPCOM.cpp -- Compute a triangulation using TOPCOM

   Copyright 2006 Matthias Koeppe
   Derived from the TOPCOM source file points2placingtriang.cpp,
     which is copyright 1999 Joerg Rambau  

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

listCone *
triangulation_to_cones(const SimplicialComplex &sc,
		       const Vertex &vertex,
		       const vector <listVector *> &rays)
{
  listCone *cones = NULL;
  for (SimplicialComplex::const_iterator i = sc.begin(); i!=sc.end(); ++i) {
    const Simplex &simplex = *i;
    listCone *cone = createListCone();
    cone->vertex = new Vertex(vertex);
    for (Simplex::const_iterator j = simplex.begin(); j!=simplex.end(); ++j) {
      size_t index = *j;
      cone->rays = new listVector(rays[index]->first, cone->rays);
    }
    cone->rest = cones;
    cones = cone;
  }
  return cones;
}

listCone *
triangulate_cone_with_TOPCOM(listCone *cone, int numOfVars)
{
  int num_rays = lengthListVector(cone->rays);

  Matrix matrix(numOfVars, num_rays); // for TOPCOM
  vector<listVector *> rays(num_rays); // for us
  size_type i, j;
  listVector *ray;
  for (j = 0, ray = cone->rays; ray!=NULL; j++, ray = ray->rest) {
    rays[j] = ray;
    for (i = 0; i<numOfVars; i++) {
      mpz_class mpz = convert_ZZ_to_mpz(ray->first[i]);
      matrix(i, j) = Rational(mpz.get_mpz_t());
    }
  }
  PointConfiguration points(matrix);
  //std::cout << matrix << endl;
  if (points.rank() > points.no()) {
    std::cerr << "rank must not be larger than no of points." << std::endl;
    exit(1);
  }
  bool preprocess = false;
  Chirotope chiro(points, preprocess);
  PlacingTriang triang(chiro);
  std::cout << "Triangulation: " << triang << std::endl;
//   for (i = 0; i<=numOfVars+1; i++) {
//     std::cout << "Triangulation[" << i << "]: " << triang[i] << std::endl;
//   }
  return triangulation_to_cones(triang, *cone->vertex, rays);
}
