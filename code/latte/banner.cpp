#include "config.h"
#include "banner.h"

using namespace std;

void
latte_banner(ostream &s)
{
  cout << "This is LattE " << VERSION << endl;
  cout << "Derived from the official LattE release 1.2 (August 18, 2005)" << endl;
  cout << "as available from http://www.math.ucdavis.edu/~latte/" << endl << endl;
}

