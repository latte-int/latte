/* RegularTriangulation.cpp -- Support for regular triangulations
	       
   Copyright 2006, 2007 Matthias Koeppe

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

#include <cassert>
#include <vector>
#include "latte_gmp.h"
#include "latte_random.h"
#include "RegularTriangulationWithCddlib.h"

using namespace std;

vector<listVector *>
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

typedef void
height_function_type(mpq_t height, const vec_ZZ &ray, void *data);

void
random_height(mpq_t height, const vec_ZZ &ray, void *data)
{
  int max_height = * (int *) data;
  int h = uniform_random_number(1, max_height);
  mpq_set_si(height, h, 1);
}

void
biased_random_height(mpq_t height, const vec_ZZ &ray, void *data)
{
  int percentage = * (int *) data;
  int x = uniform_random_number(0, 99);
  mpq_set_si(height, (x < percentage) ? 2 : 1, 1);
}

void
delone_height(mpq_t height, const vec_ZZ &ray, void *data)
{
  ZZ h;
  int i;
  for (i = 0; i<ray.length(); i++) {
    h += ray[i] * ray[i];
  }
  mpq_class hq = convert_ZZ_to_mpq(h);
  mpq_set(height, hq.get_mpq_t());
}

void
prescribed_height(mpq_t height, const vec_ZZ &ray_vector, void *data)
{
  prescribed_height_data *height_data = (prescribed_height_data *) data;
  ZZ h;
  h = 0;
  listVector *ray;
  int index;
#if 0
  cerr << "ray_vector: " << ray_vector << endl;
#endif
  for (ray = height_data->special_rays, index = 0; ray != NULL; ray = ray->rest, index++) {
#if 0
    cerr << "special ray: " << ray->first << endl;
#endif
    if (ray->first == ray_vector) {
      h = (*height_data->special_heights)[index];
      break;
    }
  }
  mpq_class hq = convert_ZZ_to_mpq(h);
  mpq_set(height, hq.get_mpq_t());
}
