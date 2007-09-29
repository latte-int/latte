/* SpecialSimplex.cpp -- Check for a special simplex using CPLEX
	       
   Copyright 2007 Matthias Koeppe

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

#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cassert>

#include <cplex.h>

#include "latte_gmp.h"
#include "latte_random.h"
#include "SpecialSimplex.h"
#include "triangulation/RegularTriangulationWith4ti2.h"

#include "dual.h"
#include "print.h"

using namespace std;

static double
convert_ZZ_to_double(const ZZ &zz)
{
  mpz_class mpz = convert_ZZ_to_mpz(zz);
  return mpz.get_d();
}

listCone *
FindSpecialSimplex(listCone *cone, int numOfVars)
{
  CPXENVptr env;
  int status;
  env = CPXopenCPLEX(&status);
  if (status != 0) {
    cerr << "Failed to obtain CPLEX environent." << endl;
    abort();
  }

  int num_rays = lengthListVector(cone->rays);

  CPXLPptr lp = CPXcreateprob(env, &status, "repr");
  if (status != 0) abort();
  
  status = CPXchgprobtype(env, lp, CPXPROB_MILP);
  if (status != 0) abort();

  /* Fill equations that express that e_n is a linear combination of
     the rays, using variables x_i as multipliers. */
  status = CPXnewrows(env, lp, numOfVars, /*rhs:*/ NULL, /*sense:*/ NULL,
		      /*rngval:*/ NULL, /*rownames:*/ NULL);
  if (status != 0) abort();

  {
    double *obj = new double[num_rays];
    int i;
    for (i = 0; i<num_rays; i++)
      obj[i] = 1.0;
    status = CPXnewcols(env, lp, num_rays, obj,
			/*lb:*/ NULL, /*ub:*/ NULL,
			/*ctype:*/ NULL, /*colname:*/NULL);
    delete[] obj;
    if (status != 0) abort();
  }

  listVector *ray;
  int j;
  for (ray = cone->rays, j = 0; ray!=NULL; ray = ray->rest, j++) {
    int i;
    for (i = 0; i<numOfVars; i++) {
      status = CPXchgcoef(env, lp, i, j, convert_ZZ_to_double(ray->first[i]));
      if (status != 0) abort();
    }
    /* Add constraints x_i <= y_i <= M x_i */
    double M = 1000;
    int beg[2];
    int ind[4];
    double val[4];
    char sense[2];
    beg[0] = 0;  ind[0] = j; val[0] = +1; ind[1] = num_rays + j; val[1] = -M; sense[0] = 'L';
    beg[1] = 2;  ind[2] = j; val[2] = -1; ind[3] = num_rays + j; val[3] = +1; sense[1] = 'L';
    status = CPXaddrows(env, lp, 1, 2, 4, /* rhs: */ NULL, sense, beg, ind, val,
			/* colname: */ NULL, /* rowname: */ NULL);
    if (status != 0) abort();
    int index = num_rays + j;
    char ctype = 'B';
    status = CPXchgctype(env, lp, 1, &index, &ctype);
    if (status != 0) abort();
  }
  { /* Enter -1000 * e_n as a right-hand side.
       FIXME: Actually need to try both positive and negative. */
    int index = numOfVars - 1;
    double value = -1000.0;
    status = CPXchgrhs(env, lp, 1, &index, &value);
    if (status != 0) abort();
  }
  /* Add a cardinality constraint. */
  {
    int beg = 0;
    int *ind = new int[num_rays];
    double *val = new double[num_rays];
    double rhs = numOfVars;
    char sense = 'E';
    int i;
    for (i = 0; i<num_rays; i++) {
      ind[i] = num_rays + i;
      val[i] = 1.0;
    }
    status = CPXaddrows(env, lp, 0, 1, num_rays, &rhs, &sense, &beg, ind, val,
			/* colname: */ NULL, /* rowname: */ NULL);
    delete[] ind;
    delete[] val;
    if (status != 0) abort();
  }
  status = CPXwriteprob(env, lp, "special-simplex.lp", "LP");
  if (status != 0) abort();

  status = CPXmipopt(env, lp);
  if (status != 0) abort();
  
  int stat = CPXgetstat(env, lp);
  if (stat != CPXMIP_OPTIMAL) {
    cerr << "Did not find special simplex." << endl;
    exit(1);
  }

  /* Inspect which rays form the special simplex. */

  listCone *special = createListCone();
  special->vertex = new Vertex(*cone->vertex);
  {
    double *x = new double[num_rays];
    status = CPXgetmipx(env, lp, x, num_rays, 2 * num_rays - 1);
    if (status != 0) abort();
    int i;
    listVector *ray;
    for (ray = cone->rays, i = 0; i<num_rays; ray = ray->rest, i++) {
      if (fabs(x[i] - 1.0) < 0.1) {
	special->rays = new listVector(ray->first, special->rays);
      }
    }
    delete[] x;
  }
  assert(lengthListVector(special->rays) == numOfVars);
  
  return special;
}

void
special_height(mpq_t height, const vec_ZZ &ray, void *data)
{
  listCone *special_cone = (listCone *) data;
  
  int max_height = 10000;
  listVector *r;
  int h = 1;
  for (r = special_cone->rays; r!=NULL; r = r->rest, h++) {
    if (r->first == ray) {
      /* Have a special ray, put it on low height. */
      mpq_set_si(height, h, 1);
      return;
    }
  }
  /* Choose a random height. */
  h = uniform_random_number(10000, 10000 + max_height);
  mpq_set_si(height, h, 1);
}

static bool
facets_ok(listCone *cone, int numOfVars)
{
  listVector *facet;
  for (facet = cone->facets; facet!=NULL; facet=facet->rest) {
    if (facet->first[numOfVars - 1] == 0)
      return false;
  }
  return true;
}

static void
check_facets(listCone *cone, int numOfVars)
{
  if (!facets_ok(cone, numOfVars)) {
    cerr << "The following cone has bad facets." << endl;
    printCone(cone, numOfVars);
    //abort();
  }
}

class FacetCheckingConeTransducer : public ConeTransducer {
public:
  FacetCheckingConeTransducer() {}
  int ConsumeCone(listCone *cone);
};

int FacetCheckingConeTransducer::ConsumeCone(listCone *cone)
{
  int numOfVars = cone->rays->first.length();
  if (cone->facets == NULL)
    computeDetAndFacetsOfSimplicialCone(cone, numOfVars);
  check_facets(cone, numOfVars);
  return consumer->ConsumeCone(cone);
}

void
special_triangulation_with_subspace_avoiding_facets
(listCone *cone, BarvinokParameters *Parameters, ConeConsumer &consumer)
{
  listCone *special_cone = FindSpecialSimplex(cone, Parameters->Number_of_Variables);
  cerr << "Found special cone: " << endl;
  printListCone(special_cone, Parameters->Number_of_Variables);
  FacetCheckingConeTransducer checking_transducer;
  checking_transducer.SetConsumer(&consumer);
  triangulate_cone_with_4ti2(cone, Parameters,
			     special_height, special_cone,
			     Parameters->Number_of_Variables, checking_transducer);
}
