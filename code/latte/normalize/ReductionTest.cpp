/* ReductionTest.cpp -- Check reducibility of vectors
	       
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

#include <cstring>
#include <cstdlib>
#include <iostream>

#include "config.h"
#include "ReductionTest.h"

#ifdef HAVE_CPLEX
#include "ReductionTestWithCPLEX.h"
#endif

using namespace std;

ReductionTest::~ReductionTest()
{
}

int
ReductionTest::GetDimension()
{
  return -1;
}


NoReductionTest::NoReductionTest(const ReductionTestFactory &)
{
}

bool
NoReductionTest::IsReducible(const vector<int> &v)
{
  return false;
}


ReductionTestFactory::ReductionTestFactory()
  : type(NoReduction), verbose(false)
{
}

void
ReductionTestFactory::show_options(ostream &stream)
{
  stream << "  --reduction={none,cplex,facets}          Use a reduction method." << endl
	 << "  --reduction-verbose                      Talk about the reduction." << endl
	 << "  --reduction-rays-file=FILE               Reduce using generators (for --reduction=cplex)." << endl
	 << "  --reduction-facets-file=FILE             Reduce using facets (for --reduction=facets)." << endl;
}

bool
ReductionTestFactory::parse_option(const char *arg)
{
  if (strncmp(arg, "--reduction=", 12) == 0) {
    if (strcmp(arg + 12, "none") == 0)
      type = NoReduction;
    else if (strcmp(arg + 12, "cplex") == 0
	     || strcmp(arg + 12, "CPLEX") == 0)
      type = ReductionWithCPLEX;
    else if (strcmp(arg + 12, "facets") == 0)
      type = ReductionWithFacets;
    else {
      cerr << "Unknown reduction type: " << arg + 12 << endl;
      exit(1);
    }
  }
  else if (strcmp(arg, "--reduction-verbose") == 0) 
    verbose = true;
  else if (strncmp(arg, "--reduction-rays-file=", 22) == 0) {
    rays_filename = arg + 22;
  }
  else if (strncmp(arg, "--reduction-facets-file=", 24) == 0) {
    facets_filename = arg + 24;
  }
  else
    return false;
  return true;
}

ReductionTest *
ReductionTestFactory::CreateReductionTest()
{
  switch (type) {
  case NoReduction:
    return new NoReductionTest(*this);
  case ReductionWithCPLEX:
#ifdef HAVE_CPLEX
    return new ReductionTestWithCPLEX(*this);
#else
    cerr << "ReductionWithCPLEX not compiled in, sorry." << endl;
    exit(1);
#endif
  case ReductionWithFacets:
    cerr << "ReductionWithFacets not written yet, sorry." << endl;
    exit(1);
  }
}
