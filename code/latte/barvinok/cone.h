/***********************************************************************
  Author: Ruriko Yoshida
  July 24th, 2002
  Update: October 4th, 2002
  This is a program for Barvinok's decomposition of cones.
  This is a class file.

************************************************************************/
#ifndef CONE__H
#define CONE__H

#include <fstream.h>
#include <cstdlib>
#include <cstring>
#include <time.h>
#include <NTL/mat_ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_RR.h>
#include <NTL/vec_RR.h>
#include <NTL/config.h>
#include <NTL/LLL.h>
#include <NTL/HNF.h>
#include <NTL/ZZ.h>


 class Cone {
 public:
  mat_ZZ generator;
  bool sign;
  int number;
  Cone(){
  sign = 0;
  }
  Cone( const Cone& cone ){
     mat_ZZ one = cone.generator;
     generator = one;
/*       int onesign; */
     sign = cone.sign;
     }
  ~Cone(){
     generator.kill();
     sign = 0;}
  friend Cone Equal(const Cone & cone){
  Cone tmp;
  tmp.generator = cone.generator;
  tmp.sign = cone.sign;
 return tmp;
}
  friend ostream& operator << (ostream& out, const Cone& rec){
    long m = 0, n = 0;
    m = rec.generator.NumCols();
    n = rec.generator.NumRows();
    out << "Sign is: " << endl;
    out << rec.sign << endl; 
    out << "Cone is: " << endl;
    for(long i = 0; i < m; i++){
       for(long j = 0; j < n; j++)
         out << rec.generator[i][j] << " ";
       out << endl;
       }
    return out;
    }
  };


vec_ZZ ComputeOmega( const mat_ZZ & B, long m, int x, int y);
vec_ZZ ComputeOmega_2(mat_ZZ &B, long m);



#endif

ZZ lcm(const ZZ&, const ZZ&);


