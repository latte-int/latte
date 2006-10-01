// This is a -*- C++ -*- header file.

/* latte_ntl.cpp -- Interface to NTL

   Copyright 2006 Matthias Koeppe

   This file is part of LattE.
   
   LattE is free software; you can redistribute it and/or modify it
   under the terms of the version 2 of the GNU General Public License
   as published by the Free Software Foundation.

   LattE is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with LattE; if not, write to the Free Software Foundation,
   Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
*/

// Include all required header files from NTL,
// and if NTL is configured to use C++ namespaces,
// use the NTL namespace.

#ifndef LATTE_NTL_H
#define LATTE_NTL_H

#include <NTL/config.h>
#include <NTL/ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include <NTL/RR.h>
#include <NTL/mat_RR.h>
#include <NTL/vec_RR.h>
#include <NTL/HNF.h>
#include <NTL/LLL.h>

#if defined(NTL_STD_CXX) || defined(NTL_PSTD_NNS)
using namespace NTL_NAMESPACE;
#endif

// Additional functions.
void
InnerProductModulo(ZZ &result, const vec_ZZ &a, const vec_ZZ &b, const ZZ &module);

#endif
