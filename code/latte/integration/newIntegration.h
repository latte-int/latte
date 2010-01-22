#ifndef NEWINTEGRATION_H
#define NEWINTEGRATION_H
#include <NTL/vec_vec_ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/ZZ.h>

NTL_CLIENT

struct simplexZZ
{
	int d;
	vec_vec_ZZ s;
	ZZ v;
};

void delSpace(string &line);
void convertToSimplex(simplexZZ &mySimplex, string line);
void integrateList(ZZ &a, ZZ &b, string line, simplexZZ mySimplex);
void integrateFlatVector(ZZ& numerator, ZZ& denominator, const linFormSum &forms , simplexZZ mySimplex);
#endif
