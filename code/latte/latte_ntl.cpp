#include <cassert>
#include "latte_ntl.h"

void
InnerProductModulo(ZZ &result, const vec_ZZ &a, const vec_ZZ &b, const ZZ &module)
{
#if 1
  InnerProduct(result, a, b);
  rem(result, result, module);
#else
  /* very slow! */
  result = 0;
  assert(a.length() == b.length());
  int dimension = a.length();
  int i;
  ZZ ai, bi, p;
  for (i = 0; i<dimension; i++) {
    rem(ai, a[i], module);
    rem(bi, b[i], module);
    MulMod(p, ai, bi, module);
    AddMod(result, result, p, module);
  }
#endif
}
