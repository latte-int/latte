/* ReadSubcones.cpp -- Read/write a simple file format describing a collection of subcones
	       
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

#include <cassert>
#include <vector>
#include <iomanip>
#include "ReadSubcones.h"

using namespace std;

static void check_stream(const istream &f, const char *fileName, const char *proc)
{
  if (!f.good()) {
    cerr << "Read error on input file " << fileName << " in " << proc << "." << endl;
    exit(1);
  }
};

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

static listCone *cone_from_ray_indicator(const vector<listVector *> &ray_array,
					 const vector<bool> &ray_indicator,
					 listCone *master_cone)
{
  listCone *result = createListCone();
  assert(ray_array.size() == ray_indicator.size());
  int i;
  for (i = 0; i<ray_array.size(); i++) {
    if (ray_indicator[i]) {
      result->rays = new listVector(ray_array[i]->first, result->rays);
    }
  }
  result->vertex = new Vertex(*master_cone->vertex);
  return result;
}

void
ReadSubcones(listCone *master_cone, int numOfVars,
	     ifstream &f, const char *fileName,
	     ConeConsumer &consumer)
{
  int numOfSubcones, numOfRays;
  f >> numOfSubcones >> numOfRays;
  check_stream(f, fileName, "ReadSubcones");
  if (numOfRays != lengthListVector(master_cone->rays)) {
    cerr << "Wrong subcones file dimensions:"
	 << "Master cone has " << lengthListVector(master_cone->rays) << " rays, "
	 << "subcones file specified " << numOfRays << " rays."
	 << endl;
    exit(1);
  }
  int i;
  vector<listVector *> rays = ray_array(master_cone);
  vector<bool> ray_indicator(numOfRays);
  for (i = 0; i<numOfSubcones; i++) {
    int j;
    for (j = 0; j<numOfRays; j++) {
      int x;
      f >> x;
      if (x != 0 && x != 1) {
	cerr << "Subcone file contains bad numbers, only 0 and 1 are defined." << endl;
	exit(1);
      }
      ray_indicator[j] = x;
    }
    check_stream(f, fileName, "ReadSubcones");
    consumer.ConsumeCone(cone_from_ray_indicator(rays, ray_indicator,
						 master_cone));
  }
}

void
ReadSubcones(listCone *master_cone, int numOfVars,
	     const string &fileName,
	     ConeConsumer &consumer)
{
  ifstream file(fileName.c_str());
  ReadSubcones(master_cone, numOfVars, file, fileName.c_str(), consumer);
}


SubconeReadingConeProducer::SubconeReadingConeProducer
(listCone *a_master_cone, const string &a_filename, int a_size_estimate)
  : master_cone(a_master_cone),
    filename(a_filename),
    size_estimate(a_size_estimate)
{
}

void SubconeReadingConeProducer::Produce(ConeConsumer &consumer)
{
  if (size_estimate)
    consumer.SetNumCones(size_estimate);
  ReadSubcones(master_cone, master_cone->rays->first.length(),
	       filename, consumer);
}


SubconePrintingConeConsumer::SubconePrintingConeConsumer(const listCone *master_cone, 
							 const std::string & filename)
  : cone_count(0), stream(filename.c_str())
{
  listVector *ray;
  int index;
  for (ray = master_cone->rays, index = 0; ray!=NULL; ray=ray->rest, index++) {
    std::map<vector<mpz_class>, int>::value_type p(convert_vec_ZZ_to_mpz_vector(ray->first),						  
						   p.second = index);
    index_map.insert(p);
  }
  // We fill in the correct number of lines later.
  stream << setw(16) << left << -1 << setw(0) << right
	 << " " << index << endl;
}

SubconePrintingConeConsumer::~SubconePrintingConeConsumer()
{
  // Go to the beginning and fill in the correct number of lines.
  stream.seekp(0, ios::beg);
  stream << setw(16) << left << cone_count;
  
}

int SubconePrintingConeConsumer::ConsumeCone(listCone *cone)
{
  cone_count++;
  int num_master_rays = index_map.size();
  vector<bool> ray_indicator(num_master_rays);
  listVector *ray;
  for (ray = cone->rays; ray!=NULL; ray=ray->rest) {
    vector<mpz_class> v = convert_vec_ZZ_to_mpz_vector(ray->first);
    std::map<vector<mpz_class>, int>::iterator i
      = index_map.find(v);
    if (i == index_map.end()) {
      cerr << "Cone has a ray that does not belong to the master cone; cannot print it as a subcone."
	   << endl;
      exit(1);
    }
    int index = (*i).second;
    ray_indicator[index] = true;
  }
  int index;
  for (index = 0; index<num_master_rays; index++) {
    stream << ray_indicator[index] << " ";
  }
  stream << endl;
  freeCone(cone);
  return 1;
}
