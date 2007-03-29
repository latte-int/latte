// This is a -*- C++ -*- header file.

/* ReductionTest.h -- Check reducibility of vectors
	       
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

#ifndef REDUCTIONTEST_H
#define REDUCTIONTEST_H

#include <vector>
#include <string>
#include <ostream>

class ReductionTest {
public:
  virtual ~ReductionTest();
  virtual bool IsReducible(const std::vector<int> &v) = 0;
};

class ReductionTestFactory {
public:
  enum ReductionType {
    NoReduction, ReductionWithCPLEX, ReductionWithFacets
  };
  ReductionType type;
  std::string rays_filename;
  std::string facets_filename;
  bool verbose;
public:
  ReductionTestFactory();
  void show_options(std::ostream &stream);
  bool parse_option(const char *arg);
  ReductionTest *CreateReductionTest();
};

class NoReductionTest : public ReductionTest {
public:
  NoReductionTest(const ReductionTestFactory &);
  virtual bool IsReducible(const std::vector<int> &v);
};

#endif
