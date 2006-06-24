/*******************************************************************
   Author: Ruriko Yoshida
   Date: July 25th, 2002
   Update: Febrary 3rd, 2003
   This program computes Barvinok's decomposition of a cone.
   This program is for the project "LattE."

*********************************************************************/
#include <list>

#include <fstream>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <math.h>
#include <time.h>

#include "Cone.h"
#include "barvinok.h"
#include "../myheader.h"
#include "../ramon.h"
#include "../RudyResNTL.h"
#include "rational.h"
#include "convert.h"
#include "dual.h"
#include "config.h"
#ifdef HAVE_EXPERIMENTS
#include "Irrational.h"
#endif

 /* Note:  We are dealing with the "Row space" of the
    input matrix due to NTL. */

/**********************************************************************/
vec_ZZ CheckOmega( const mat_ZZ & U, vec_ZZ & Z){

  int m;
  m = U.NumCols();
  mat_ZZ Dummy;
  Dummy.SetDims(m + 1, m);
  Dummy[0] = Z;
  ZZ d;

  for(int i = 0; i < m; i++)
    Dummy[i + 1] = U[i];

  mat_ZZ dummy;
  image(d, Dummy, dummy);

  int flag = 1, number = 0; 

  for(int i = 0; i <= m; i++)
      if(dummy[0][i] >= 0) number++;
  if(number == (m + 1))  flag = 0;

  if(flag != 0){
  number = 0;
  for(int i = 0; i <= m; i++)
     if(dummy[0][i] <= 0) number++;
  if(number == (m + 1)) flag = 0;
  }     
  if(flag == 0){
    Z = - Z; 
  }
  Dummy.kill();
  dummy.kill();
  return Z;

}
/**********************************************************************/
 
void MatrixGCD(mat_ZZ & B, long & m){
  ZZ gcds[m];
  for(int i = 1; i <= m; i++)
    for(int j = 1; j <= m; j++)
      if(B(i, j) != 0)
	gcds[i-1] = GCD(gcds[i-1], B(i, j));
  for(int i = 1; i <= m; i++)
    for(int j = 1; j <= m; j++)
      if(B(i, j) != 0)
	B(i, j) = B(i, j) / gcds[i-1];

}
/**********************************************************************/

/* Likewise barvinok_Single, but the cone is given as a listCone;
   the function consumes the cone. */
int
barvinok_DFS(listCone *cone, Single_Cone_Parameters *Parameters);

int
barvinok_Single(mat_ZZ B, Single_Cone_Parameters *Parameters,
		rationalVector *vertex)
{
	//cout << "barvinok_Single Called." << endl;;
	
	long m, n;
  	m = B.NumRows();
  	n = B.NumCols();

   	if( m != n)
   	{	
       		cerr << "Input must be square. " << endl;
       		exit(2);
   	}

   	ZZ D = determinant(B);

         if( D == 0)
   	{
       		cerr << "Input must be linearly independent. " << endl;
       		exit(3);
   	}

	 Parameters->Total_Simplicial_Cones++;
	 
   	/* The following routine is to get the minimal
      	integral generators for the cone.  */

   	MatrixGCD(B, m);

   	listCone *dummy = createListCone();
   	dummy->coefficient = 1;
	dummy->determinant = D;
	dummy->vertex = copyRationalVector(vertex);
	dummy->rays = transformArrayBigVectorToListVector(B, m, n);

	switch (Parameters->decomposition) {
	case BarvinokParameters::DualDecomposition:
	  // Keep the dual cones during Barvinok decomposition
	  computeDetAndFacetsOfSimplicialCone(dummy, n);
	  break;
	case BarvinokParameters::IrrationalPrimalDecomposition:
	  // Do Barvinok decomposition on the primal cones.
	  dualizeBackCones(dummy, Parameters->Number_of_Variables);
#ifdef HAVE_EXPERIMENTS
	  irrationalizeCone(dummy, Parameters->Number_of_Variables);
#else
	  cerr << "Irrationalization not configured in this build."
	       << endl;
	  exit(1);
#endif
	  break;
	default:
	  cerr << "Unknown BarvinokParameters::decomposition";
	  abort();
	}
	
	int result;
	result = barvinok_DFS(dummy, Parameters);

	return result;
}
	
/* Let GENERATOR and MAT be the same matrix, with determinant DET.
   Copy the vector Z into each row of the matrix (we are dealing with
   the row space) and compute the determinant of the resulting matrix;
   store the determinants in DETS[].  When any determinant is
   larger-or-equal than DET in absolute value, stop the computation
   and return false.  Otherwise return true.
*/
static bool
computeAndCheckDeterminants(const mat_ZZ &generator, const ZZ &Det,
			    const vec_ZZ &Z, int m, 
			    mat_ZZ &mat, ZZ Dets[])
{
  ZZ absDet = abs(Det);
  for (int i = 1; i <= m; i++) {
    /* Copy in the row */
    for(int j = 1; j <= m; j++)
      mat(i, j) = Z(j);
    /* Compute and store the determinant. */
    determinant(Dets[i - 1], mat);
    /* Restore the original row */
    for(int j = 1; j <= m; j++)
      mat(i, j) = generator(i, j);
    if (abs(Dets[i - 1]) >= absDet) 
      return false;
  }
  return true;
}

/* Decompose the cone spanned by GENERATOR (which has determinant DET)
   according to Barvinok's theory into M (the dimension) many CONES
   and store their determinants in DETS.

   Entries with Det[i] == 0 have Cones[i] == NULL (we don't generate
   lower-dimensional cones).
*/
static void
barvinokStep(const listCone *Cone, 
	     listCone *Cones[], ZZ Dets[],
	     int m)
{
  mat_ZZ generator = createConeDecMatrix(Cone, m, m);
  mat_ZZ dual = createFacetMatrix(Cone, m, m);
  /* ComputeOmega(const mat_ZZ &, long& ) computes
     an integral vector in the parallelogram. */
  vec_ZZ Z = ComputeOmega(generator, dual, m, 0, 0);
  Z = CheckOmega(generator, Z);
     
  mat_ZZ mat = generator;
  bool success
    = computeAndCheckDeterminants(generator, Cone->determinant, Z,
				  m, mat, Dets);
  if (!success) {
    cout << "Second loop... " << endl;
    Z = ComputeOmega(generator, dual, m, 2, 2);
    Z = CheckOmega(generator, Z);
    success = computeAndCheckDeterminants(generator, Cone->determinant, Z,
					  m, mat, Dets);
    assert(success);
  }

  for(int i = 0; i < m; i++) {
    if (Dets[i] == 0)
      Cones[i] = NULL;
    else {
      Cones[i] = createListCone();
      {
	/* Create the rays: */
	/* Copy in the row */
	for(int j = 1; j <= m; j++)
	  mat(i+1, j) = Z(j);
	Cones[i]->rays
	  = transformArrayBigVectorToListVector(mat, m, m);
	/* Restore the original row */
	for(int j = 1; j <= m; j++)
	  mat(i+1, j) = generator(i+1, j);
      }
      Cones[i]->determinant = Dets[i];
      {
	/* Compute the sign: */
	long signDet = sign(Cone->determinant);
	long signDeti = sign(Dets[i]);
	Cones[i]->coefficient = Cone->coefficient * signDet * signDeti;
      }
      Cones[i]->vertex = copyRationalVector(Cone->vertex);
      computeDetAndFacetsOfSimplicialCone(Cones[i], m);
    }
  }
}

int barvinok_DFS(listCone *C, Single_Cone_Parameters *Parameters)
{
  ZZ absDet;
  switch (Parameters->decomposition) {
  case BarvinokParameters::DualDecomposition:
    absDet = abs(C->dual_determinant);
    break;
  case BarvinokParameters::IrrationalPrimalDecomposition:
    absDet = abs(C->determinant);
    break;
  default:
    cerr << "Unknown BarvinokParameters::decomposition";
    abort();
  }
       	
  if (absDet == 0) {
    //cout << "barvinok_DFS: Det = 0." << endl;
    return 1;	
  }		     
  else if (Parameters->max_determinant == 0
	   || absDet <= Parameters->max_determinant) {
    Parameters->Total_Uni_Cones += 1;
    if ( Parameters->Total_Uni_Cones % 1000 == 0)
      cout << Parameters->Total_Uni_Cones
	   << (Parameters->max_determinant == 0
	       ? " simplicial cones done."
	       : (Parameters->max_determinant == 1
		  ? " unimodular cones done."
		  : " low-index cones done."))
	   << endl;
    switch (Parameters->decomposition) {
    case BarvinokParameters::DualDecomposition:
      C = dualizeBackCones(C, Parameters->Number_of_Variables);
      return Parameters->ConsumeCone(C);
    case BarvinokParameters::IrrationalPrimalDecomposition:
      return Parameters->ConsumeCone(C);
    default:
      cerr << "Unknown BarvinokParameters::decomposition";
      abort();
    }
  }	     
  
  //cout << "barvinok_DFS: non-uni cone." << endl;
     
  int result = 1;
  long m = Parameters->Number_of_Variables;

  ZZ Dets[m];	     
  listCone *cones1 [m];

  barvinokStep(C, cones1, Dets, m);
  
  ZZ max;
  max = -1;

#ifdef SHOWDETS
  cout << "Index " << absDet << " -> ";
#endif
  for(int i = 0; i < m; i++)
    {
      Dets[i] = abs(Dets[i]);
#ifdef SHOWDETS
      cout << Dets[i] << ", ";
#endif
      if(Dets[i] > max)
	max = Dets[i];
      
      if (Dets[i] > 0) {
	Parameters->Current_Simplicial_Cones_Total ++;

#ifdef HAVE_EXPERIMENTS
	if (Parameters->decomposition == BarvinokParameters::IrrationalPrimalDecomposition)
	  checkConeIrrational(cones1[i], Parameters->Number_of_Variables);
#endif
      }
    }
#ifdef SHOWDETS
  cout << endl;
#endif

  int current;
  ZZ min;

  if (Parameters->Current_Simplicial_Cones_Total > Parameters->Max_Simplicial_Cones_Total)
    Parameters->Max_Simplicial_Cones_Total = Parameters->Current_Simplicial_Cones_Total;
	
  do {
    min = max + 1;
    current = -1;
    for(int j = 0; j < m; j++) {
      if(Dets[j] < min && Dets[j] != 0)
	{
	  current = j;
	  min = Dets[j];
	}
    }
    if (current >= 0) {
      Dets[current] = 0; // mark done
      if(barvinok_DFS(cones1[current], Parameters) == -1)
	result = -1;
      Parameters->Current_Simplicial_Cones_Total--;
    }
  } while (current >= 0 && result == 1);
  freeListCone(C);
  return result;
}

