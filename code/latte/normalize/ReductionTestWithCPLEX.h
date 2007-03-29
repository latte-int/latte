// This is a -*- C++ -*- header file.

/* ReductionTestWithCPLEX.h -- Check reducibility of vectors using CPLEX
	       
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

#ifndef REDUCTIONTESTWITHCPLEX_H
#define REDUCTIONTESTWITHCPLEX_H

#include <cplex.h>

#include "ReductionTest.h"

class ReductionTestWithCPLEX : public ReductionTest {
public:
  ReductionTestWithCPLEX(const ReductionTestFactory &data);
  virtual ~ReductionTestWithCPLEX();
  virtual bool IsReducible(const std::vector<int> &v);
  virtual int GetDimension();
private:
  bool verbose;
  CPXENVptr env;
  CPXLPptr lp;
  int dim;
  int num_generators;
  double *multipliers;
  int *indices;
  char *ctype;
};

#endif
