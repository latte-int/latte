#include <fstream.h>
#include <cstdlib>
#include <cstring>
#include <math.h>
#include <algorithm>
#include <time.h>
#include <iostream>

using namespace std;

int main(){


  int times = 0, seed, sum = 0;

  cin >> times >> seed;
  int A[times];
  srand(seed);
  for(int i = 0; i < times; i++)
   { A[i] = rand()/100;
   sum += A[i];
   cout << A[i] << " ";
   } cout << endl;
  cout << sum / 2 << endl;
  return 0;
}
