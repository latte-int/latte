// This is a -*- C++ -*- header file.

/* ReadSubcones.h -- Read/write a simple file format describing a collection of subcones
	       
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

#ifndef READSUBCONES_H
#define READSUBCONES_H

#include <fstream>
#include <map>
#include "latte_gmp.h"
#include "cone.h"
#include "cone_consumer.h"

void
ReadSubcones(listCone *master_cone, int numOfVars,
	     ifstream &f, const char *fileName,
	     ConeConsumer &consumer);

void
ReadSubcones(listCone *master_cone, int numOfVars,
	     const string &fileName,
	     ConeConsumer &consumer);

class SubconeReadingConeProducer : public ConeProducer {
  listCone *master_cone;
  string filename;
  int size_estimate;
public:
  SubconeReadingConeProducer(listCone *a_master_cone, const string &a_filename,
			     int a_size_estimate = 0);
  void Produce(ConeConsumer &consumer);
};

class IncrementalVectorFileWriter {
public:
  long int num_vectors;
  IncrementalVectorFileWriter(const std::string &filename, int a_dimension);
  ~IncrementalVectorFileWriter();
  void WriteVector(const vec_ZZ &v);
  void WriteVector(const std::vector<bool> &v);
  void WriteVector(const std::vector<int> &v);
  void UpdateNumVectors();
private:
  std::ofstream stream;
  int dimension;
};

class SubconePrintingConeConsumer : public ConeConsumer {
public:
  int cone_count;
  SubconePrintingConeConsumer(const listCone *master_cone, const std::string &filename);
  ~SubconePrintingConeConsumer();
  int ConsumeCone(listCone *cone);
private:
  IncrementalVectorFileWriter *file_writer;
  std::map<vector<mpz_class>, int> index_map;
  vector< vec_ZZ > master_rays;
};

#endif

