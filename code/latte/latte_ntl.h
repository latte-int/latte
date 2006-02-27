// This is a -*- C++ -*- header file.
//
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

#endif
