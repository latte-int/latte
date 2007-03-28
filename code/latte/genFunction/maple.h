/* maple.h -- Create Maple input

   Copyright 2002-2004 Jesus A. De Loera, David Haws, Raymond
      Hemmecke, Peter Huggins, Jeremy Tauzer, Ruriko Yoshida

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

rationalVector* addRationalVectorsWithUpperBoundOne(rationalVector*, 
						    rationalVector*, int);
rationalVector* subRationalVector(rationalVector*, rationalVector*, int);
listVector* readListVector(char*);
vec_ZZ movePoint(vec_ZZ, rationalVector*, rationalVector*, vec_ZZ*, int, int);
listVector* pointsInParallelepiped(rationalVector*, listVector*, int);
void writeTermToFile(FILE*, vec_ZZ, int);
void writeTermOfGeneratingFunctionToFile(FILE*, listCone*, int);
void createGeneratingFunctionAsMapleInput(const char*, listCone*, int);
void createGeneratingFunctionAsMapleInputGrob(listCone* cones, 
					      int numOfVars, ofstream & out);
