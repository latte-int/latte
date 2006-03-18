/* ----------------------------------------------------------------- */
/*                                                                   */
/* LattE (Lattice Point Enumeration)                                 */
/*                                                                   */
/* Computing all lattice points in fundamental parallelepiped        */
/*                                                                   */
/* Author     : Raymond Hemmecke, Ruriko Yoshida                     */
/*                                                                   */
/* Created    : 07-JUN-02                                            */
/* Last Update: 20-DEC-02 by Rudy                                    */
/*                                                                   */
/* ----------------------------------------------------------------- */
#include "config.h"
#include "../myheader.h"
#include "../cone.h"
#include "../print.h"
#include "../ramon.h"
#include "../rational.h"
#include "../vertices/cdd.h"
#include <string>

#include "IntCombEnum.h"
#include "NTL_to_LiDIA.h" 
#include "matrix_ops.h" 
#include <NTL/mat_ZZ.h>
#include <assert.h>

using namespace std;

/* Forward declare functions internal to this module */
vec_ZZ get_integer_comb(const mat_ZZ &, int *);
void get_multipliers_from_snf(const mat_ZZ &, int *);
mat_ZZ get_U_star(listCone *);
vec_ZZ translate_lattice_point(const vec_ZZ &, const mat_ZZ &, const mat_ZZ &,
   listCone *);  

/* ----------------------------------------------------------------- */
vec_ZZ movePoint(vec_ZZ x, rationalVector *coeffsX, 
		 rationalVector *coeffsVertex, vec_ZZ *matrix, 
		 int numOfRays, int numOfVars) {
  int i,j;
  vec_ZZ z, movement;
  rationalVector *difference;

  difference=subRationalVector(coeffsX,coeffsVertex,numOfVars);

  movement=createVector(numOfRays);
  for (i=0; i<numOfRays; i++) {
    if (difference->denominator[i]==1) 
      movement[i]=-difference->enumerator[i];
    else {
      movement[i]=-(difference->enumerator[i]/difference->denominator[i]);
      if (difference->enumerator[i]<0) movement[i]=movement[i]+1;
    }
  }

  z=copyVector(x,numOfVars);

  for (i=0; i<numOfRays; i++) {
    if (movement[i]!=0) {
      for (j=0; j<numOfVars; j++) 
/* Note that each COLUMN of matrix constitutes an extreme ray, not
   each row! */
	z[j]=z[j]+movement[i]*matrix[j][i];
    }
  }

  return (z);
}

/*
 * Description:
 * pointsInParallelepiped: 1) Find Smith Normal form of cone basis. In other
 * words, find a basis v and w such that v_i = n_i*w_i. 2) Ennumerate all
 * lattice points in fund||piped by taking all (bounded: 0 <= k <= n_{i-1})
 * integer combinations of w_i, and translating them accordingly
 *
 * Parameters:
 * rationalVector *vertex: vertex of cone
 * listVector *rays: cone is defined by non-negative linear combinations of
 *    these rays.
 * listVector *facets: not applicable
 * int numOfVars: dimension of cone
 *
 * Return: listVector* of lattice points
 */
listVector* pointsInParallelepiped(listCone *cone, int numOfVars) {
   mat_ZZ snf_U, B, B_inv, C, U_star;
   mat_ZZ U = convert_listVector_to_mat_ZZ(cone->rays);
   rationalVector *vertex = cone->vertex;
   vec_ZZ lat_pt, trans_lat_pt;
   int rows;
   int *n, *next;
   IntCombEnum *iter_comb;
   listVector *lat_points = NULL;

   cout << "Computing Smith-Normal form...\n";
   /* get Smith Normal form of matrix, Smith(A) = BAC */
   snf_U = SmithNormalForm(U, B, C);
   /*snf_U = SmithNormalForm(cone->rays, B, C);*/
   rows = snf_U.NumRows();

   /* extract n_i such that v_i = n_i*w_i from Smith Normal form */
   n = new int[rows]; 
   get_multipliers_from_snf(snf_U, n);

   /*
    * Smith(A) = BAC, and the diagonal entries of Smith(A) are the multipliers
    * n_i such that w_i*n_i= v_i. Therefore, AC = V, and B^-1 = W. Since we
    * must take the integer combinations of W, we must calculate B^-1.
    */
   B_inv = inv(B);

   /* U* = (U^-1)^T. Thus, we can calculate U* from the facets */
   U_star = get_U_star(cone);
   
   /*
    * enumerate lattice points by taking all integer combinations
    * 0 <= k <= (n_i - 1) of each vector.
    */
   cout << "Enumerating lattice points...\n";
   iter_comb = new IntCombEnum(n, rows);
   iter_comb->decrementUpperBound();
   while((next = iter_comb->getNext())) {
   cout << "Calculating integer combination ...\n";
      lat_pt = get_integer_comb(B_inv, next);
   cout << "Finished calculating, translating lattice point...\n";
      trans_lat_pt = translate_lattice_point(lat_pt, U, U_star, cone);
   cout << "Finished translating, adding lattice point...\n";
      lat_points = appendVectorToListVector(trans_lat_pt, lat_points);
   cout << "Sucessfully added lattice points...\n";
   }

   cout << "Cleanup...\n";
   /* cleanup */
   delete iter_comb;
   cout << "deleted iter_comb...\n";
   delete [] n;
   cout << "deleted n[]...\n";

   return (lat_points);
}

/*
 * calculates a integer combination of the specified lattice basis
 * and stores the result in the vector v. Because of incompatible datatypes,
 * it is easiest to simply do this directly.
 */
vec_ZZ get_integer_comb(const mat_ZZ & lat_basis, int *scalar) {
   int i,j;
   long conv_l;
   int row = lat_basis.NumRows();
   int col = lat_basis.NumCols();
   vec_ZZ v = createVector(col);

   for (i = 0; i < row; i++) {
      v[i] = 0;
      for (j = 0; j < col; j++) {
         conv_l = to_long(lat_basis[i][j]);
         v[i] += conv_l*scalar[j]; 
      }
   }
   return (v);
}

mat_ZZ get_U_star(listCone *cone) {
   /* Note that < RAY_i, FACET_j > = -FACET_DIVISOR_i * DELTA_{i,j}. */
   mat_ZZ U_star = convert_listVector_to_mat_ZZ(cone->facets);

   return (U_star);
}

/*
 * Translate point m using the following formula:
 * m' = v + sum {<m - v, u_i*>}u_i, where {} means the factional part 
 */
vec_ZZ translate_lattice_point(const vec_ZZ& m, const mat_ZZ & U,
   const mat_ZZ & U_star, listCone * cone) {
   int cols = U.NumCols();
   int len = m.length();
   mat_ZZ U_trans = transpose(U);
   mat_ZZ U_star_trans = transpose(U_star);
   rationalVector *v = cone->vertex;
   rationalVector w;

   /* intialized to one */
   ZZ dot_prod_enum = to_ZZ(0);
   ZZ dot_prod_denom = to_ZZ(1);

   /* initialized to zero */
   ZZ frac_enum, frac_denom, q, r, neg;
   vec_ZZ scaled_u_i, m_trans;
   vec_ZZ offset_enum, sum_enum, m_trans_enum;
   vec_ZZ offset_denom, sum_denom, m_trans_denom;

   /* allocate memory for intermediate calculation vectors */
   offset_enum.SetLength(len);
   offset_denom.SetLength(len);
   m_trans_enum.SetLength(len);
   m_trans_denom.SetLength(len);
   sum_enum.SetLength(len);
   sum_denom.SetLength(len);
   m_trans.SetLength(len);

cout << "translate_lattice_point: U\n";
print_debug_matrix(U);
cout << "translate_lattice_point: U_star\n";
print_debug_matrix(U_star);

   /* initialize sum */
   sum_enum = v->enumerator;
   sum_denom = v->denominator;

   /* w = m - v */
cout << "translate_lattice_point: m - v\n";
   for (int j = 0; j < len; j++) {
      /* m - x/y = (my - x)/y */
cout << "translate_lattice_point: m[" << j << "] = " << m[j] << "\n";
cout << "translate_lattice_point: v->denom[" << j << "] = " << v->denominator[j] << "\n";
cout << "translate_lattice_point: v->enum[" << j << "] = " << v->enumerator[j] << "\n";
      offset_enum[j] = (m[j] * v->denominator[j]) - v->enumerator[j];
      offset_denom[j] = v->denominator[j];
cout << "translate_lattice_point: offset_enum[" << j << "] = " << offset_enum[j] << "\n";
cout << "translate_lattice_point: offset_denom[" << j << "] = " << offset_denom[j] << "\n";
   }

   for (int i = 0; i < cols; i++) {
cout << "translate_lattice_point: < m - v, u_" << i << "* >\n";
      /* reset/initialize variables */
      dot_prod_enum = to_ZZ(0);
      dot_prod_denom = to_ZZ(1);

      /* < m - v, u_i* > */
      for (int j = 0; j < len; j++) {
         /* neg = - U_star[j][i] */
         //NTL::negate(neg, U_star[j][i]);
         /* x/y + r/s*z = x/y + rz/s = xs + rzy/ys */
         dot_prod_enum = (dot_prod_enum * offset_denom[j]) +
            (dot_prod_denom * offset_enum[j] * U_star[j][i]);
         dot_prod_denom *= (offset_denom[j]);
cout << "translate_lattice_point: facet_divisor[" << j << "] = " << cone->facet_divisors[j] << "\n";
      }
cout << "translate_lattice_point: before dot_prod_denom = " << dot_prod_denom << "\n";
      dot_prod_denom *= cone->facet_divisors[i];
cout << "translate_lattice_point: dot_prod_enum = " << dot_prod_enum << "\n";
cout << "translate_lattice_point: dot_prod_denom = " << dot_prod_denom << "\n";
   
      /*
       * {<m - v, u_i*>} : take the fractional part of this quantity
       * q = floor(dot_prod_nume/dot_prod_denom)
       * r = dot_prod_enum - dot_prod_denom*q
       */
cout << "translate_lattice_point: {< m - v, u_" << i << "* >}\n";
      DivRem(q, r, dot_prod_enum, dot_prod_denom);
      frac_enum = r;
      frac_denom = dot_prod_denom;
cout << "translate_lattice_point: frac_enum = " << frac_enum << "\n";
cout << "translate_lattice_point: frac_denom = " << frac_denom << "\n";

      /* {<m - v, u_i*>}u_i */
      scaled_u_i = U_trans[i] * frac_enum;

      /* sum {<m - v, u_i*>}u_i */
      for (int j = 0; j < len; j++) {
         /* w/z + x/y = (wy + xz)/zy */
         sum_enum[j] = (sum_enum[j] * frac_denom) +
            (sum_denom[j] * scaled_u_i[j]);
         sum_denom[j] *= frac_denom;
cout << "translate_lattice_point: sum_enum[" << j << "] = " << sum_enum[j] << "\n";
cout << "translate_lattice_point: sum_denom[" << j << "] = " << sum_denom[j] << "\n";
      }
   }

   for (int j = 0; j < len; j++) {
      /*
       * This MUST be an integer point!! This MUST divide evenly!
       * m_trans = floor(sum_enum/sum_denom)
       * r = sum_enum - sum_denom*q
       */
cout << "translate_lattice_point: sum_enum[" << j << "] = " << sum_enum[j] << "\n";
cout << "translate_lattice_point: sum_denom[" << j << "] = " << sum_denom[j] << "\n";
      DivRem(m_trans[j], r, sum_enum[j], sum_denom[j]);
cout << "translate_lattice_point: m_trans[" << j << "] = " << m_trans[j] << "\n";
      assert(IsZero(r));
   }
cout << "translate_lattice_point: EXITING!!!!!";
   return (m_trans);
}
  
/*
 * If S is the Smith Normal Form of A, then the multipliers n_i such that
 * v_i = n_iw_i are on the diagonal
 */
void get_multipliers_from_snf(const mat_ZZ & snf, int *n) {
   int j;
   int row = snf.NumRows();

   for (j = 0; j < row; j++) {
      //cout << "get_multipliers_from_snf:: snf[" << j << "," << j << "] = " << snf[j][j] << "\n"; 
      n[j] = to_int(snf[j][j]);
   }
   return;
}


/* ----------------------------------------------------------------- */
listVector* pointsInParallelepipedOfUnimodularCone(rationalVector *vertex, 
						  listVector *rays, 
						  int numOfVars) {
  int i,j,k,numOfRays;
  ZZ a,b;
  vec_ZZ lambda,w,z;
  vec_ZZ *matrix, *originalMatrix;
  listVector *points, *tmp;
  rationalVector *coeffs;

/*  printf("Computing point in parallelepiped of unimodular cone.\n"); */

/*  printRationalVector(vertex,numOfVars); */

  //points=createListVector(createVector(numOfVars));

  numOfRays=lengthListVector(rays);
  if (numOfRays!=numOfVars) {
    printf("Cone is NOT simplicial!\n");
    exit(0);
  }

/*  printf("numOfRays = %d, numOfVars = %d\n", numOfRays, numOfVars); */

  matrix=createArrayVector(numOfVars);
  originalMatrix=createArrayVector(numOfVars);

  for (i=0; i<numOfVars; i++) matrix[i]=createVector(numOfRays);

  k=0;
  tmp=rays;
  while (tmp) {
    for (i=0; i<numOfVars; i++) matrix[i][k]=(tmp->first)[i];
    k++;
    tmp=tmp->rest;
  }

 for (i=0; i<numOfRays; i++) 
    originalMatrix[i]=copyVector(matrix[i],numOfRays);

  for (i=0; i<numOfVars; i++) {
    for (j=0; j<numOfRays; j++) {
      matrix[i][j]=matrix[i][j]*(vertex->denominator[i]);
    }
  }
  w=copyVector(vertex->enumerator,numOfVars);

  /* Now we have to solve: matrix * x = w */

  coeffs=solveLinearSystem(matrix,w,numOfVars,numOfRays);

  lambda=createVector(numOfRays);

/* Note that we make heavy use of numOfVars=numOfRays! That is, we
   assume the cone to be simplicial! */

  for (i=0; i<numOfVars; i++) {
    a=coeffs->denominator[i];
    b=coeffs->enumerator[i];

    if (((b/a)*a)==b) {
      lambda[i]=b/a;
    } else {
      if (((b<0) && (a>0)) || ((b>0) && (a<0))) {
//          lambda[i]=-abs(b)/abs(a); 
         lambda[i]=abs(b)/abs(a); 
         lambda[i]=-lambda[i]; 
     } else 
	lambda[i]=b/a;
      if (b>0) lambda[i]=lambda[i]+1;
    }
  }

//  cout << "lambda = ";
//  printVector(lambda,numOfVars);

  z=createVector(numOfVars);
  for (i=0; i<numOfVars; i++) z[i]=0;

  for (i=0; i<numOfRays; i++) {
    for (j=0; j<numOfVars; j++) {
      z[j]=z[j]+lambda[i]*originalMatrix[j][i];
    }
  }

//  cout << "point = ";
//  printVector(z,numOfVars);

  //points->rest=createListVector(z);
  
  coeffs->denominator.kill ();
  coeffs->enumerator.kill ();

  delete coeffs;
  
  points=createListVector(z);
  delete [] matrix;
  delete [] originalMatrix;
  w.kill();
  z.kill();
  lambda.kill();

  return(points);
}
/* ----------------------------------------------------------------- */

void computePointsInParallelepiped(listCone *cone, int numOfVars)
{
/*
   int p1, p2;
   listVector *tmp;
   tmp = pointsInParallelepiped(cone,numOfVars);
   p1 = lengthListVector(tmp);
   printListVector(tmp, numOfVars);
   cout << "p1 = " << p1 << "\n";
   cone->latticePoints = pointsInParallelepipedOfUnimodularCone(cone->vertex,cone->rays,numOfVars);
   p2 = lengthListVector(cone->latticePoints);
   printListVector(cone->latticePoints, numOfVars);
   cout << "p2 = " << p2 << "\n";
   assert(p1 == p2);
   assert(isEqual(tmp, cone->latticePoints));
*/
   if (abs(cone->determinant) != 1) {
      cout << "Processing cone with determinant " << cone->determinant << endl;
      cone->latticePoints = pointsInParallelepiped(cone, numOfVars);
      cout << "Finished processing cone..." ;
   } else {
      cone->latticePoints = pointsInParallelepipedOfUnimodularCone(
         cone->vertex, cone->rays, numOfVars);
   }
}

void computePointsInParallelepipeds(listCone *cones, int numOfVars)
{
  listCone *tmp = cones;
  int Cones_Processed_Count = 0;
  while (tmp) {
    computePointsInParallelepiped(tmp, numOfVars);
    tmp=tmp->rest;
    Cones_Processed_Count++;
    if ((Cones_Processed_Count % 1000) == 0 )
      cout << Cones_Processed_Count << " cones processed." << endl;
  }
}
