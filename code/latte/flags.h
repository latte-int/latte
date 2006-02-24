/* This is a -*- C++ -*- header file. */
#ifndef FLAGS_H
#define FLAGS_H 1

#define PRINT 0x1
#define OUTPUT 0x6
#define OUTPUT0 0x2
#define OUTPUT1 0x4
#define DUAL_APPROACH 0x8
#define DECOMPOSE 0x10
#define LOAD	0x20
#define SAVE    0x40

struct BarvinokParameters {
  // Whether we use the
  //   - traditional LattE monomial substitution z_i |-> (1 + s)^(lambda_i) 
  //   - or the exponential substitution         z_i |-> exp(t lambda_i)
  enum { PolynomialSubstitution, ExponentialSubstitution } substitution;
  // The maximum determinant of cones that we do not subdivide
  // further.  Set to 1 to subdivide until we reach unimodular cones
  // only.
  int max_determinant;
};

#endif
