// This is a -*- C++ -*- header file.

/* cone_consumer.h -- Cone producer/consumer pattern
	       
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

#ifndef CONE_CONSUMER_H
#define CONE_CONSUMER_H

#include "cone.h"

// ConeConsumers consume cones.

class ConeConsumer {
public:
  // Take CONE and consume it.
  virtual int ConsumeCone(listCone *cone) = 0;
  // Take notice of the expected number of cones to be produced.
  virtual void SetNumCones(size_t num_cones);
  virtual ~ConeConsumer();
};

class CollectingConeConsumer : public ConeConsumer {
public:
  CollectingConeConsumer();
  listCone *Collected_Cones;
  int ConsumeCone(listCone *cone);
};

// A consumer that just prints to a file.
class PrintingConeConsumer : public ConeConsumer {
  ofstream stream;
public:
  int cone_count;
  PrintingConeConsumer(string filename);
  int ConsumeCone(listCone *cone);
};

// ConeProducers produce cones.

class ConeProducer {
public:
  virtual void Produce(ConeConsumer &consumer) = 0;
  virtual ~ConeProducer();
};

class SingletonConeProducer : public ConeProducer {
public:
  SingletonConeProducer(listCone *a_cone);
  void Produce(ConeConsumer &consumer);
private:
  listCone *cone;
};

class ListConeReadingConeProducer : public ConeProducer {
public:
  ListConeReadingConeProducer(const string &a_filename, int a_size_estimate = 0);
  void Produce(ConeConsumer &consumer);
private:
  string filename;
  int size_estimate;
};

// ConeTransducer consume cones and produce other cones;
// which they send to another consumer.

class ConeTransducer : public ConeConsumer {
public:
  ConeTransducer();
  int ConsumeCone(listCone *cone) = 0;
  void SetConsumer(ConeConsumer *a_consumer);
protected:
  ConeConsumer *consumer;
};

// The three types of operations can be composed.

class CompositeConeProducer : public ConeProducer {
public:
  CompositeConeProducer(ConeProducer *a_producer, 
			ConeTransducer *a_transducer);
  void Produce(ConeConsumer &consumer);
private:
  ConeProducer *producer;
  ConeTransducer *transducer;
};

class CompositeConeTransducer : public ConeTransducer {
public:
  CompositeConeTransducer(ConeTransducer *a_first,
			  ConeTransducer *a_second);
  int ConsumeCone(listCone *cone);
private:
  ConeTransducer *first;
  ConeTransducer *second;
};

class CompositeConeConsumer : public ConeConsumer {
public:
  CompositeConeConsumer(ConeTransducer *a_transducer,
			ConeConsumer *a_consumer);
  int ConsumeCone(listCone *cone);
private:
  ConeTransducer *transducer;
  ConeConsumer *consumer;
};

ConeProducer *
compose(ConeProducer *a_producer, ConeTransducer *a_transducer);

ConeTransducer *
compose(ConeTransducer *a_transducer_1, ConeTransducer *a_transducer_1);

ConeConsumer *
compose(ConeTransducer *a_transducer, ConeConsumer *a_consumer);

#endif
