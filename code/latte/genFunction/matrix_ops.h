/* matrix_ops.cpp -- Smith normal form

   Copyright 2006 Susan Margulies

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

#ifndef MATRIX_OPS_H 
#define MATRIX_OPS_H 

#include "cone.h"
#include "barvinok/barvinok.h"

BarvinokParameters::SmithFormType
smith_form_type_from_name(const char *name);

void
show_standard_smith_option(ostream &stream);

bool
parse_standard_smith_option(const char *arg,
			    BarvinokParameters *params);





mat_ZZ
SmithNormalFormLidia(const mat_ZZ &, mat_ZZ &, mat_ZZ &);

/** mat_ZZ **/
/** SmithNormalForm(listVector *, mat_ZZ &, mat_ZZ &); **/

mat_ZZ
SmithNormalFormIlio(const mat_ZZ &, mat_ZZ &, mat_ZZ &);

mat_ZZ
SmithNormalForm(const mat_ZZ &, mat_ZZ &, mat_ZZ &, BarvinokParameters *params);


#endif

