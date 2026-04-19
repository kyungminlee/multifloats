#include "../src/multifloats.hh"
#include "../src/multifloats_c.h"
#include <cstdio>
#include <cmath>
namespace mf = multifloats;

static void trace(char const *tag, mf::float64x2 x) {
  printf("  %s: hi=%g lo=%g signbit(hi)=%d signbit(lo)=%d\n",
         tag, x._limbs[0], x._limbs[1],
         std::signbit(x._limbs[0]), std::signbit(x._limbs[1]));
}

int main() {
  double nz = std::copysign(0.0, -1.0);
  mf::float64x2 a(2.0), b;
  b._limbs[0] = nz; b._limbs[1] = 0.0;
  trace("a", a); trace("b", b);

  mf::float64x2 ab = a * b;
  trace("a*b", ab);

  mf::float64x2 sum = ab + ab;
  trace("a*b + a*b", sum);

  mf::float64x2 neg_sum;
  neg_sum._limbs[0] = -sum._limbs[0];
  neg_sum._limbs[1] = -sum._limbs[1];
  trace("-(a*b+a*b)", neg_sum);

  // Also via mf operator: DD unary -
  mf::float64x2 neg_via_op = mf::float64x2(0.0) - sum;
  trace("0 - sum", neg_via_op);

  // Now: 1 - a*a + b*b
  mf::float64x2 re_1mz2 = mf::float64x2(1.0) - a * a + b * b;
  trace("Re(1-z^2)", re_1mz2);

  // Full casin for (2, -0).
  complex64x2_t z;
  z.re.hi = 2.0; z.re.lo = 0.0;
  z.im.hi = nz;  z.im.lo = 0.0;
  complex64x2_t r = casindd(z);
  printf("\ncasin(2,-0): re=(%g,%g) im=(%g,%g) signbit(im.hi)=%d\n",
         r.re.hi, r.re.lo, r.im.hi, r.im.lo, std::signbit(r.im.hi));
  return 0;
}
