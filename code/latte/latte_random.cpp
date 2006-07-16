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
