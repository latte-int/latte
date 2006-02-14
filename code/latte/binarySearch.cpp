#include <iostream>
#include <vector>
using namespace std;

  template <class Comp>
  int binarySearch( const vector<Comp> & a, const Comp & x)
  {
      int low = 0, high = a.size() - 1;
      while(low <= high)
      {
      	int mid = (low + high) / 2;

         if(a[ mid ] < x)
            low = mid + 1;

         else if(x < a[ mid ])
             high = mid - 1;
         els return mid; //found

         }
      return NOT_FOUND;

  }