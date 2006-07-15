#include "config.h"
#include "banner.h"

using namespace std;

void
latte_banner(ostream &s)
{
  cout << "This is LattE macchiato " << VERSION << endl;
  cout << "Available from http://www.math.uni-magdeburg.de/~mkoeppe/latte/" << endl;
  cout << "Derived from the official LattE release 1.2 (August 18, 2005)" << endl;
  cout << "as available from http://www.math.ucdavis.edu/~latte/" << endl << endl;
}

