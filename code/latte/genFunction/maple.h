// This is a -*- C++ -*- header file.

/* maple.h -- Create Maple input

   Copyright 2002-2004 Jesus A. De Loera, David Haws, Raymond
      Hemmecke, Peter Huggins, Jeremy Tauzer, Ruriko Yoshida
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

#ifndef GENFUNCTION_MAPLE_H
#define GENFUNCTION_MAPLE_H

#include <string>
#include <fstream>

#include "cone.h"
#include "cone_consumer.h"
#include "barvinok/barvinok.h"

void writeTermToFile(ofstream & out, const vec_ZZ &, int);
void writeTermOfGeneratingFunctionToFile(ofstream & out, listCone*, int);
void createGeneratingFunctionAsMapleInput(const char*, listCone*, int);
void createGeneratingFunctionAsMapleInputGrob(listCone* cones, 
					      int numOfVars, ofstream & out);

class GeneratingFunctionWritingConeConsumer : public ConeConsumer {
public:
  GeneratingFunctionWritingConeConsumer(const std::string &genfun_filename);
  virtual int ConsumeCone(listCone *cone);
private:
  ofstream genfun_stream;
  bool first_term;
};
  
#endif
