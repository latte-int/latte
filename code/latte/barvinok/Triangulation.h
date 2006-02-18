/***********************************************************************************
   Author: Ruriko Yoshida
   Date: August 18th, 2002
   Update: October 9th, 2002
   Project for LattE
   This program computes a triangulation of a cone
   Input: the number of vectors, the dimension of the vectors and
          the cone generators.
   Output: the cone triangulation in R^n.
************************************************************************************/
#ifndef TRIANGULATION__H
#define TRIANGULATION__H
#include <list>
#include <NTL/mat_ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_RR.h>
#include <NTL/vec_RR.h>
#include <NTL/config.h>
#include <NTL/LLL.h>
#include <NTL/HNF.h>
#include <NTL/ZZ.h>
using namespace std;
using namespace NTL_NAMESPACE;

int Triangulation(const mat_ZZ &, const int &, const int &, char*, list< int >&);


int Triangulation_Load_Save(const mat_ZZ &, const int &, const int &, char*, list< int >&, char *File_Name, int Cone_Index, unsigned int Flags);


#endif
