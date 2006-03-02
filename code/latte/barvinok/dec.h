/*********************************************************** -*- C++ -*- */
#ifndef BARVINOK_DEC_H
#define BARVINOK_DEC_H

#include "../flags.h"

listCone* readListCone(rationalVector*, int);

listCone*
decomposeCones(listCone*, int, unsigned int Flags,
	       char *File_Name, int max_determinant = 1);

// Added by Peter/David
// 
void decomposeCones_Single (listCone *, int, int degree, unsigned int flags, char *File_Name);

#endif
