// Compute π to full double-double precision as 4·atan(1) and print it.
//
// Build against an installed multifloats package (see README.md):
//   cmake -S example -B build/example -DCMAKE_PREFIX_PATH=<install-prefix>
//   cmake --build build/example
#include <multifloats/float64x2.h>  // links with multifloats::multifloats
#include <iostream>

int main() {
  using namespace multifloats;
  float64x2 x(1.0);
  float64x2 pi = float64x2(4.0) * atan(x);  // π to ~106 bits
  std::cout << "pi = " << to_string(pi, 32) << '\n';
  return 0;
}
