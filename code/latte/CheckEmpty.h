/*                                                                    *
 * Author: Ruriko Yoshida                                             *
 * Date: Febrary 10th, 2004                                           *
 * Update: Febrary 10th, 2004                                         *
 * This is for checking whether the input polytope is empty or not.   *
 *
 */
#ifndef CHECKEMPTY__H
#define CHECKEMPTY__H

#include <string.h>
#include <stdio.h>
#include <NTL/mat_ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_RR.h>
#include <NTL/vec_RR.h>
#include <NTL/config.h>
#include <NTL/LLL.h>
#include <NTL/HNF.h>
#include <NTL/ZZ.h>

#include "myheader.h"
#include "barvinok/dec.h"
#include "barvinok/barvinok.h"
#include "barvinok/Cone.h"
#include "barvinok/ConeDecom.h"
#include "barvinok/Triangulation.h"
#include "vertices/cdd.h"
#include "genFunction/maple.h"
#include "genFunction/piped.h"
#include "cone.h"
#include "dual.h"
#include "RudyResNTL.h"
#include "Grobner.h"

#include "preprocess.h"
#include "print.h"
#include "ramon.h"
#include "rational.h"
#include "timing.h"
#include "flags.h"

#include "IntegralHull.h"
#include "ReadingFile.h"
#include "binarySearchIP.h"

void CheckEmpty(char * filename);

#endif
