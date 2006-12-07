/* latte2ine.cpp -- Convert a LattE file to an ine file
	       
   Copyright 2006 Matthias Koeppe

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

#include <cstdlib>
#include <iostream>
#include "ReadLatteStyle.h"

using namespace std;

int main(int argc, char **argv)
{
  if (argc != 1) {
    cerr << "usage: " << argv[0] << " < INPUTFILE.latte > OUTPUTFILE.ine" << endl;
    exit(1);
  }
  dd_MatrixPtr m = ReadLatteStyleMatrix(cin, /* vrep: */ false,
					"standard input");
  dd_WriteMatrix(stdout, m);
  dd_FreeMatrix(m);
  return 0;
}
