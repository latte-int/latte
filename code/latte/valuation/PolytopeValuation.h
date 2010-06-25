/*
 * PolytopeValuation.h
 *
 *  Created on: Jun 25, 2010
 *      Author: bedutra
 */

#ifndef POLYTOPEVALUATION_H_
#define POLYTOPEVALUATION_H_


#include <string.h>
#include <stdio.h>

#include "config.h"
#include "barvinok/dec.h"
#include "barvinok/barvinok.h"
#include "barvinok/Triangulation.h"
#include "vertices/cdd.h"
#include "genFunction/maple.h"
#include "genFunction/piped.h"
#include "cone.h"
#include "dual.h"
#include "RudyResNTL.h"
#include "Residue.h"
#include "Grobner.h"
//  #include "jesus.h"
#include "preprocess.h"
#include "print.h"
#include "ramon.h"
#include "rational.h"
#include "timing.h"
#include "flags.h"
//#include "testing.h"
#include "IntegralHull.h"
#include "ReadingFile.h"
#include "binarySearchIP.h"
#include "CheckEmpty.h"
#include "ProjectUp.h"

#include "banner.h"
#include "convert.h"
#include "latte_system.h"

#include "gnulib/progname.h"






class PolytopeValuation
{
public:
	PolytopeValuation();
	virtual ~PolytopeValuation();


	// A B C D E F G H I J K L M N O P Q R S T U V W X Y Z


	void findVolume(listCone *cones, int numOfVars, unsigned int Flags,
		       char *File_Name, int max_determinant,
		       bool dualize,
		       BarvinokParameters::DecompositionType decomposition);

	void findVolumeSingle(listCone *cones, int numOfVars,
			unsigned int flags, char *File_Name, int max_determinant, bool dualize,
			BarvinokParameters::DecompositionType decomposition);

	listCone * decomposeAndReturnCones(listCone *cones, bool dualize,
				   Standard_Single_Cone_Parameters &param);

};

#endif /* POLYTOPEVALUATION_H_ */
