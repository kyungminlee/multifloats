#include "../src/multifloats.hh"
#include "../src/multifloats_c.h"
#include <cstdio>
#include <cmath>
int main() {
  complex64x2_t z;
  z.re.hi = 1.0; z.re.lo = 0.0;
  z.im.hi = 0.0; z.im.lo = 0.0;
  complex64x2_t r = catanhdd(z);
  printf("catanh(1+0i): re=(%g,%g) im=(%g,%g)\n",
         r.re.hi, r.re.lo, r.im.hi, r.im.lo);
  printf("  isinf(re.hi)=%d, re.hi>0=%d\n",
         std::isinf(r.re.hi), r.re.hi > 0);

  // also (-1, 0)
  z.re.hi = -1.0;
  r = catanhdd(z);
  printf("catanh(-1+0i): re=(%g,%g) im=(%g,%g)\n",
         r.re.hi, r.re.lo, r.im.hi, r.im.lo);
  return 0;
}
