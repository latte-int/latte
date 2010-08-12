/* dual.cpp -- Dualize polyhedral cones

   Copyright 2002 Raymond Hemmecke, Ruriko Yoshida
   Copyright 2006, 2007 Matthias Koeppe

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

#include <cassert>
#include "config.h"
#include "cone.h"
#include "print.h"
#include "ramon.h"
#include "rational.h"
#include "vertices/cdd.h"
#include "barvinok/dec.h"
#include "convert.h"
#include <list>
#include "latte_system.h"
#include "latte_relocatable.h"

#ifdef HAVE_FORTYTWO_LIB
#include "DualizationWith4ti2.h"
#endif

using namespace std;

BarvinokParameters::DualizationType
dualization_type_from_name(const char *name)
{
  if (strcmp(name, "cdd") == 0) return BarvinokParameters::DualizationWithCdd;
  else if (strcmp(name, "4ti2") == 0) return BarvinokParameters::DualizationWith4ti2;
  else {
    cerr << "Unknown dualization type name: " << name << endl;
    exit(1);
  }
}

bool
parse_standard_dualization_option(const char *arg,
				  BarvinokParameters *params)
{
  if (strncmp(arg, "--dualization=", 14) == 0) {
    params->dualization = dualization_type_from_name(arg + 14);
  }
  else return false;
  return true;
}

void
show_standard_dualization_option(ostream &stream)
{
  stream << "  --dualization={cdd,4ti2}" << endl;
}

static void dualizeCone_with_cdd(listCone *tmp, int numOfVars)
{
  char cddInFileName[127];
  string tmpString;
  int i,j,tmpInt,len,numOfVertices;
  ZZ x,y;
  rationalVector *w;
  listVector *rays, *rays2, *facets, *endFacets;

  assert(tmp->subspace_generators == NULL);
  
  rays=tmp->rays;
  rays2=rays;
  len=lengthListVector(rays);

  strcpy(cddInFileName,"latte_cdd.ine");

  ofstream out(cddInFileName);
  if(!out){
    cerr << "Cannot open output file in dualizeCones." << endl;
    exit(1);
  }

  out << "H-representation\n";
  out << "begin\n";
  out << lengthListVector(rays) << " " << numOfVars+1 << "integer" << endl;

  while (rays) {
    out << "0 ";
    for (i=0; i<(numOfVars); i++) out << -(rays->first)[i] << " ";
    out << endl;
    rays=rays->rest;
  }
  out << "end\n";
  out.close();

  /*      printf("Computing facets with cdd..."); */
  system_with_error_check(relocated_pathname(CDD_PATH) + " latte_cdd.ine > latte_cdd.out");
  /*      printf("done.\n"); */

  strcpy(cddInFileName,"latte_cdd.ext");

  ifstream in(cddInFileName);
  if(!in){
    cerr << "Cannot open input file in dualizeCones." << endl;
    exit(1);
  }

  while (tmpString!="begin") getline(in,tmpString);

  in >> numOfVertices >> tmpInt >> tmpString;

  facets=createListVector(createVector(numOfVars));
  endFacets=facets;
      
  for (i=0; i<numOfVertices; i++) {
    w=createRationalVector(numOfVars);
    for (j=0; j<numOfVars+1; j++) {
      x=0;
      y=0;
      ReadCDD(in,x,y);
      if (j>0) {
	w->set_entry(j-1, x, y, true /* avoid recomputation of
					integer scale */);
      }
    }
    w=normalizeRationalVector(w,numOfVars);
    endFacets->rest=createListVector(w->numerators());
    delete w;
    endFacets=endFacets->rest;
  }
  in.close();
#if 0
  system_with_error_check("rm -f latte_cdd.*");
#endif
  tmp->facets=tmp->rays;    
  tmp->rays=facets->rest;
  delete facets; // only delete the dummy head
}

void computeDetAndFacetsOfSimplicialCone(listCone *cone, int numOfVars);

void dualizeCone(listCone *tmp, int numOfVars, BarvinokParameters *params)
{
  if ((tmp->rays != NULL && tmp->facets != NULL)
      /* Both rays and facets are already computed,
	 so just swap. */
      || (tmp->rays == NULL && tmp->facets != NULL)
      /* Facets are already computed, so just swap, so we have a ray
	 representation. */) {
    swap(tmp->determinant, tmp->dual_determinant);
    swap(tmp->rays, tmp->facets);
    swap(tmp->subspace_generators, tmp->equalities);
  }
  else if (lengthListVector(tmp->rays) == params->Number_of_Variables) {
    /* We assume full-dimensional cones, so this must be a simplicial
       cone. */
    computeDetAndFacetsOfSimplicialCone(tmp, params->Number_of_Variables);
    swap(tmp->determinant, tmp->dual_determinant);
    swap(tmp->rays, tmp->facets);
    swap(tmp->subspace_generators, tmp->equalities);
  }
  else {
    switch (params->dualization) {
    case BarvinokParameters::DualizationWithCdd:
      dualizeCone_with_cdd(tmp, params->Number_of_Variables);
      break;
    case BarvinokParameters::DualizationWith4ti2:
#ifdef HAVE_FORTYTWO_LIB
      dualizeCone_with_4ti2(tmp, params->Number_of_Variables);
#else
    cerr << "DualizationWith4ti2 not compiled in, sorry."
	 << endl;
    exit(1);
#endif
      break;
    default:
      cerr << "Unknown DualizationType" << endl;
      exit(1);
    }
  }
}

void
dualizeCones(listCone *cones, int numOfVars, BarvinokParameters *params)
{
  params->dualize_time.start();
  int numOfConesDualized,numOfAllCones;
  listCone *tmp;

  cerr << "Dualizing all cones...";
  cerr.flush();

  numOfConesDualized=0;
  numOfAllCones=lengthListCone(cones);

  tmp=cones;
  while (tmp) {
    dualizeCone(tmp, numOfVars, params);
    tmp=tmp->rest;
    numOfConesDualized++;
    if (numOfConesDualized==50*(numOfConesDualized/50)) {
    	cerr << numOfConesDualized << " / " << numOfAllCones << " done.\n";
    }
  }

  cerr << "All cones are now dualized." << endl;
  //removeListVector(cones->facets);
  params->dualize_time.stop();
  cerr << params->dualize_time;
}

/* ----------------------------------------------------------------- */

void computeDetAndFacetsOfSimplicialCone(listCone *cone, int numOfVars)
{
  int i;
  listVector *rays;
  mat_ZZ Mat, Inverse;

  assert(cone->facets == NULL);
  assert(lengthListVector(cone->rays) == numOfVars);
  
  Mat.SetDims(numOfVars, numOfVars);

  rays=cone->rays;
  for(i = 0; i < numOfVars; i++) {
    Mat[i] = rays->first;
    rays = rays -> rest;
  }

  ZZ det;
  // computes d = det A and solves X*A = I*d if d != 0.
  // thus X = det A * A^{-1}, det X = (det A)^{n-1}.
  inv(det, Inverse, Mat);
  assert(det != 0);
  cone->determinant = det;
  if (det < 0) {
    // Fix the sign.
    Inverse = -Inverse;
  }
  Inverse = - transpose(Inverse);
  cone->dual_determinant
    = determinant(Inverse); // FIXME: Easier to compute
  cone->facet_divisors.SetLength(numOfVars);
  cone->facets
    = transformArrayBigVectorToListVector(Inverse,
					  numOfVars, numOfVars);
  listVector *facet;
  for (i = 0, facet = cone->facets; i<numOfVars;
       i++, facet = facet->rest) {
    /* Cancel GCD: */
    ZZ gcd;
    int j;
    for (j = 0; j<numOfVars; j++)
      gcd = GCD(gcd, facet->first[j]);
    if (gcd != 0 && gcd != 1) {
      for (j = 0; j<numOfVars; j++)
	facet->first[j] /= gcd;
      cone->dual_determinant /= gcd;
    }
    ZZ remainder;
    DivRem(cone->facet_divisors[i], remainder, abs(cone->determinant), gcd);
    assert(IsZero(remainder));
  }
}

/* ----------------------------------------------------------------- */

void computeTightInequalitiesOfCones(listCone *cones,
				     listVector *inequalities,
				     int numOfVars)
{
  listCone *cone;
  for (cone = cones; cone; cone = cone->rest) {
    if (cone->facets == NULL) {
      listVector *inequality;
      listVector *tight_inequalities = NULL;
      ZZ vertex_scale_factor;
      vec_ZZ scaled_vertex
	= scaleRationalVectorToInteger(cone->vertex->vertex, numOfVars,
				       vertex_scale_factor);
      for (inequality = inequalities; inequality; inequality = inequality->rest) {
	int i;
	ZZ sp;
	vec_ZZ &ineq = inequality->first;
	sp = vertex_scale_factor * ineq[0];
	for (i = 0; i<numOfVars; i++)
	  sp += scaled_vertex[i] * ineq[i + 1];
	if (IsZero(sp)) {
	  vec_ZZ vec;
	  vec.SetLength(numOfVars);
	  for (i = 0; i<numOfVars; i++)
	    vec[i] = -ineq[i+1];
	  tight_inequalities = new listVector(vec,
					      tight_inequalities);
	}
      }
      cone->facets = tight_inequalities;
    }
  }
}
