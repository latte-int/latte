/*******************************************************************
   Auther: Ruriko Yoshida
   Date: August 19th, 2003
	This is for computing random vertices.

*********************************************************************/
#include <list>

#include <fstream.h>
#include <cstdlib>
#include <cstring>
#include <math.h>
#include <algorithm>
#include <time.h>

#include <NTL/mat_ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_RR.h>
#include <NTL/vec_RR.h>
#include <NTL/config.h>
#include <NTL/LLL.h>
#include <NTL/HNF.h>
#include <NTL/ZZ.h>

int main(int argc, char *argv[]) {
  int dim, numOfVertices;
  dim = atoi(argv[1]);
  numOfVertices = atoi(argv[2]);
  ZZ seed;
  int seed2 = atoi(argv[3]);
  seed = atoi(argv[3]);
  SetSeed(seed);
  ZZ Num[numOfVertices];
  ZZ Dem[numOfVertices];

  for(int i = 0; i < numOfVertices; i++){
   	Num[i] = RandomBits_ZZ(20);
      Dem[i] = RandomBits_ZZ(20);
  }
  srand(seed2);
  ofstream out(argv[4]);
  out << "V-representation " << endl;
  out << "begin" << endl;
  out << numOfVertices << " " << dim + 1 << "rational" << endl;
  for(int i = 0; i < numOfVertices ; i++){
   	out << 1 << " ";
      for(int j = 0; j < dim - 1; j++){
       	if(rand()%2 == 1)
           	out << -1<< "/" << 2 << " ";
         else out << 1 << "/" << 2 << " ";
      }
      out << Num[i] << "/" << Dem[i] << endl;
  }
   out << "end" << endl;
   out << "hull" << endl;
 return 0;
}
