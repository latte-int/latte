
#include <NTL/vec_ZZ.h>
#include <NTL/ZZ.h>

void addToFraction(const int &s[50],int n, const ZZ &t, const int &index[50], const int &counter[50], ZZ &a, ZZ &b);
void enumerate(const int &s[50], int i, int k, int n, const ZZ &p, const ZZ &t, const int &index[50], const int &counter[50], ZZ &a, ZZ &b);
void computeResidue(int d, ZZ M, const vec_ZZ &innerProDiff, const ZZ &p, ZZ &a, ZZ &b);//compute residue given the inner products
//d is dimension, M is power raised to, innerProDiff is the difference of inner products w.r.t. one vertex
//p is the inner product value at this vertex
