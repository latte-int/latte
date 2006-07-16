// This is a -*- C++ -*- header file.

#ifndef LATTE_RANDOM_H
#define LATTE_RANDOM_H

#include "latte_ntl.h"

void
seed_random_generator(unsigned int seed);

int
uniform_random_number(int from, int to);

ZZ
uniform_random_number(ZZ from, ZZ to);

#endif
