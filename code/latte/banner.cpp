/* banner.cpp -- Print the LattE banner

   Copyright 2006, 2011 Matthias Koeppe

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

#include "config.h"
#include "banner.h"

using namespace std;

void
latte_banner(ostream &s)
{
  s << "This is LattE integrale " << VERSION << endl;
  s << "Available from http://www.math.ucdavis.edu/~latte/" << endl << endl;
}

