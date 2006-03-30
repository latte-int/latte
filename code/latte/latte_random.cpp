#include <stdlib.h>
#include "latte_random.h"

int
uniform_random_number(int from, int to)
{
  int range = to - from + 1;
  return from + (int) ((double)range * rand()/(RAND_MAX+1.0));
}
