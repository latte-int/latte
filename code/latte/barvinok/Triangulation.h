// This is a -*- C++ -*- header file.

/* Triangulation.cpp -- Triangulations of cones.

   Copyright 2002 Ruriko Yoshida

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

#ifndef TRIANGULATION__H
#define TRIANGULATION__H

#include <list>
#include "latte_ntl.h"

using namespace std;

int Triangulation_Load_Save(const mat_ZZ &, const int &, const int &, char*, list< int >&, char *File_Name, int Cone_Index, unsigned int Flags);


#endif
