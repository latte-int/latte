/* latte_random.cpp -- Interface to random number generators

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

#include <stdlib.h>
#include "latte_random.h"

void
seed_random_generator(unsigned int seed)
{
  srand(seed);
  // Also seed NTL's pseudo-random generator.
  ZZ z_seed;
  z_seed = seed;
  SetSeed(z_seed);
}

int
uniform_random_number(int from, int to)
{
  int range = to - from + 1;
  return from + (int) ((double)range * rand()/(RAND_MAX+1.0));
}

ZZ
uniform_random_number(ZZ from, ZZ to)
{
  ZZ range = to - from + 1;
  return from + RandomBnd(range);
}
