#include "../src/multifloats.hh"
#include "../src/multifloats_c.h"
#include <cstdio>
#include <cmath>
int main() {
  auto prt = [](char const *tag, complex64x2_t r) {
    printf("%s: re=(%g,%g) im=(%g,%g)  signbit(im.hi)=%d\n",
           tag, r.re.hi, r.re.lo, r.im.hi, r.im.lo,
           std::signbit(r.im.hi));
  };
  double pz = +0.0, nz = std::copysign(0.0, -1.0);

  for (double x : {2.0, -2.0}) {
    for (double im : {pz, nz}) {
      complex64x2_t z;
      z.re.hi = x;  z.re.lo = 0.0;
      z.im.hi = im; z.im.lo = 0.0;
      char tag[64];
      snprintf(tag, sizeof tag, "casin(%g,%c0)", x, std::signbit(im)?'-':'+');
      prt(tag, casindd(z));
    }
  }
  // Also casin in principal range: (0.5, 0.5)
  complex64x2_t z = {{0.5, 0.0}, {0.5, 0.0}};
  prt("casin(0.5,0.5)", casindd(z));
  return 0;
}
