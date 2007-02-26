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
  PrintingConeConsumer(string filename);
  int ConsumeCone(listCone *cone);
};

#endif
