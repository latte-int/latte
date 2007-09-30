/* special-simplex.cpp -- Main program for testing special simplex construction
	       
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

#include "SpecialSimplex.h"
#include "print.h"

int main(int argc, char **argv)
{
  if (argc < 2) {
    cerr << "usage: special-simplex LATTE-CONE-FILE " << endl;
    exit(1);
  }
  listCone *cone;
  /* Input in LattE cone format. */
  ifstream in(argv[argc-1]);
  if (in.bad()) {
    cerr << "special-simplex: Unable to open " << argv[argc-1] << endl;
    exit(1);
  }
  cone = readConeFromFile(in);
  if (!cone) {
    cerr << "special-simplex: Parse error in file." << endl;
    exit(1);
  }

  if (cone->rays == NULL) {
    cerr << "special-simplex: No rays." << endl;
    exit(2);
  }
  int numOfVars = cone->rays->first.length();

  printListCone(FindSpecialSimplex(cone, numOfVars), numOfVars);
}
