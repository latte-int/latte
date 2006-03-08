#include <iostream>
using namespace std;

#include "Irrational.h"

void
irrationalizeCone(listCone *cone, int numOfVars)
{
  cerr << "Irrationalizing not implemented, you should see wrong results " << endl
       << "due to lower-dimensional cones." << endl;
}

void
irrationalizeCones(listCone *cones, int numOfVars)
{
  listCone *cone;
  for (cone = cones; cone != NULL; cone=cone->rest)
    irrationalizeCone(cone, numOfVars);
}

