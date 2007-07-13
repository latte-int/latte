/* cone_consumer.cpp -- Cone producer/consumer pattern
	       
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
#include "cone_consumer.h"
#include "print.h"

ConeConsumer::~ConeConsumer()
{
}

void
ConeConsumer::SetNumCones(size_t num_cones)
{
  // Do nothing.
}

CollectingConeConsumer::CollectingConeConsumer()
  : Collected_Cones(NULL)
{
}

int CollectingConeConsumer::ConsumeCone(listCone *cone)
{
  assert(cone->rest == NULL);
  cone->rest = Collected_Cones;
  Collected_Cones = cone;
  return 1; // means "success, please continue"
}

PrintingConeConsumer::PrintingConeConsumer(string filename)
  : stream(filename.c_str()), cone_count(0)
{}

int
PrintingConeConsumer::ConsumeCone(listCone *cone)
{
  assert(cone->rest == NULL);
  int numOfVars = cone->rays->first.length();
  cone_count++;
  printConeToFile(stream, cone, numOfVars);
  freeCone(cone);
  return 1; // means "success, please continue"
}



ConeProducer::~ConeProducer()
{}

SingletonConeProducer::SingletonConeProducer(listCone *a_cone)
  : cone(a_cone)
{}

void SingletonConeProducer::Produce(ConeConsumer &consumer)
{
  consumer.ConsumeCone(cone);
}

ListConeReadingConeProducer::ListConeReadingConeProducer
(const string &a_filename, int a_size_estimate)
  : filename(a_filename), size_estimate(a_size_estimate)
{
}

void ListConeReadingConeProducer::Produce(ConeConsumer &consumer)
{
  if (size_estimate)
    consumer.SetNumCones(size_estimate);
  ifstream file(filename.c_str());
  if (!file.good()) {
    cerr << "Error opening file `" << filename << "'" << endl;
    exit(1);
  }
  readListConeFromFile(file, consumer);
}



ConeTransducer::ConeTransducer()
  : consumer(0)
{}
  
void
ConeTransducer::SetConsumer(ConeConsumer *a_consumer)
{
  consumer = a_consumer;
}

CompositeConeProducer::CompositeConeProducer(ConeProducer *a_producer, 
					     ConeTransducer *a_transducer)
  : producer(a_producer), transducer(a_transducer)
{
}

void
CompositeConeProducer::Produce(ConeConsumer &consumer)
{
  transducer->SetConsumer(&consumer);
  producer->Produce(*transducer);  
}

CompositeConeConsumer::CompositeConeConsumer(ConeTransducer *a_transducer,
					     ConeConsumer *a_consumer)
  : transducer(a_transducer), consumer(a_consumer)
{
  transducer->SetConsumer(consumer);
}

int CompositeConeConsumer::ConsumeCone(listCone *cone)
{
  return transducer->ConsumeCone(cone);
}

ConeProducer *
compose(ConeProducer *a_producer, ConeTransducer *a_transducer)
{
  return new CompositeConeProducer(a_producer, a_transducer);
}

ConeConsumer *
compose(ConeTransducer *a_transducer, ConeConsumer *a_consumer)
{
  return new CompositeConeConsumer(a_transducer, a_consumer);
}

