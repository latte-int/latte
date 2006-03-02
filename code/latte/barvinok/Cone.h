/*********************************************************** -*- c++ -*-
  Author: Ruriko Yoshida
  July 24th, 2002
  Update: October 4th, 2002
  This is a program for Barvinok's decomposition of cones.
  This is a class file.

************************************************************************/
#ifndef CONE__H
#define CONE__H

#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <time.h>
#include "latte_ntl.h"

using namespace std;

class Cone {
public:
  bool sign;
  mat_ZZ generator;
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



