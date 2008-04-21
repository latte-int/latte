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
  unsigned int i;
  for (i = 0; i<ray_array.size(); i++) {
    if (ray_indicator[i]) {
      result->rays = new listVector(ray_array[i]->first, result->rays,
				    /* index hint: */ i);
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
  consumer.SetNumCones(numOfSubcones);
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
    listCone *cone = cone_from_ray_indicator(rays, ray_indicator,
					     master_cone);
    cone->index_hint = i;
    consumer.ConsumeCone(cone);
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


IncrementalVectorFileWriter::IncrementalVectorFileWriter(const std::string &filename, int a_dimension)
  : num_vectors(0), stream(filename.c_str()), dimension(a_dimension)
{
  // We fill in the correct number of lines later.
  if (!stream.good()) {
    cerr << "Cannot write to file " << filename << endl;
    exit(1);
  }
  stream << setw(16) << left << -1 << setw(0) << right
	 << " " << dimension << endl;
}

IncrementalVectorFileWriter::~IncrementalVectorFileWriter()
{
  UpdateNumVectors();
}

void
IncrementalVectorFileWriter::WriteVector(const vec_ZZ &v)
{
  int index;
  assert(dimension == v.length());
  for (index = 0; index<dimension; index++) {
    stream << v[index] << " ";
  }
  stream << endl;
  num_vectors++;
  if (!stream.good()) {
    cerr << "Error writing to vector file" << endl;
    exit(1);
  }
}

void
IncrementalVectorFileWriter::WriteVector(const std::vector<bool> &v)
{
  int index;
  assert(dimension == v.size());
  for (index = 0; index<dimension; index++) {
    stream << v[index] << " ";
  }
  stream << endl;
  num_vectors++;
  if (!stream.good()) {
    cerr << "Error writing to vector file" << endl;
    exit(1);
  }
}

void
IncrementalVectorFileWriter::WriteVector(const std::vector<int> &v)
{
  int index;
  assert(dimension == v.size());
  for (index = 0; index<dimension; index++) {
    stream << v[index] << " ";
  }
  stream << endl;
  num_vectors++;
  if (!stream.good()) {
    cerr << "Error writing to vector file" << endl;
    exit(1);
  }
}

void
IncrementalVectorFileWriter::UpdateNumVectors()
{
  stream.seekp(0, ios::beg);
  stream << setw(16) << left << num_vectors;
  stream.seekp(0, ios::end);
  stream.flush();
  if (!stream.good()) {
    cerr << "Error writing to vector file" << endl;
    exit(1);
  }
}


SubconePrintingConeConsumer::SubconePrintingConeConsumer(const listCone *master_cone, 
							 const std::string & filename)
  : cone_count(0), master_rays(lengthListVector(master_cone->rays))
{
  listVector *ray;
  int index;
  for (ray = master_cone->rays, index = 0; ray!=NULL; ray=ray->rest, index++) {
    std::map<vector<mpz_class>, int>::value_type p(convert_vec_ZZ_to_mpz_vector(ray->first), 
						   index);
    index_map.insert(p);
    master_rays[index] = ray->first;
  }
  file_writer = new IncrementalVectorFileWriter(filename, index);
}

SubconePrintingConeConsumer::~SubconePrintingConeConsumer()
{
  delete file_writer;
}

int SubconePrintingConeConsumer::ConsumeCone(listCone *cone)
{
  cone_count++;
  int num_master_rays = index_map.size();
  vector<bool> ray_indicator(num_master_rays);
  listVector *ray;
  for (ray = cone->rays; ray!=NULL; ray=ray->rest) {
    if (ray->index_hint >= 0
	&& ray->index_hint < master_rays.size()
	&& ray->first == master_rays[ray->index_hint])
      ray_indicator[ray->index_hint] = true;
    else {
      // index_hint not available or wrong
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
  }
  file_writer->WriteVector(ray_indicator);
  freeCone(cone);
  return 1;
}
