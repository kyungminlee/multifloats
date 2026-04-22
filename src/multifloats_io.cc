// Decimal I/O for float64x2 — scientific-notation formatting and operator<<
// for std::ostream. Moved out of multifloats.hh because this path is never
// on a hot loop; inlining is worth nothing here and the header is cleaner
// without a dependency on <string> / <ostream>.

#include "multifloats.hh"

#include <cmath>
#include <cstddef>
#include <cstring>
#include <ostream>
#include <string>

namespace multifloats {

namespace {

// Renormalize (hi, lo) into canonical DD form.
inline void io_renorm(double &hi, double &lo) {
  double s = hi + lo;
  double err = lo - (s - hi);
  hi = s;
  lo = err;
}

// (hi, lo) *= d, with d a double; result renormalized.
inline void io_dd_mul_d(double &hi, double &lo, double d) {
  double p = hi * d;
  double e = std::fma(hi, d, -p);
  lo = lo * d + e;
  hi = p;
  io_renorm(hi, lo);
}

} // anonymous namespace

std::string to_string(float64x2 const &x, int precision) {
  double hi = x._limbs[0];
  double lo = x._limbs[1];

  if (std::isnan(hi) || std::isnan(lo)) return "nan";
  if (std::isinf(hi)) return hi > 0 ? "inf" : "-inf";
  if (hi == 0.0) return std::signbit(hi) ? "-0e+00" : "0e+00";

  if (precision < 1) precision = 1;
  if (precision > 34) precision = 34;

  const bool neg = hi < 0.0;
  if (neg) { hi = -hi; lo = -lo; }

  // Decimal exponent estimate. log10 may be off by 1 due to rounding; a
  // post-scale fixup corrects.
  int e10 = static_cast<int>(std::floor(std::log10(hi)));

  // Scale to [1, 10). Multiplying by 10 is exact in DD; multiplying by 0.1
  // is not, but each iteration retains ~106 bits of precision, so 308
  // iterations (max for dp range) lose at most a handful of ulps — well
  // below the visible 34 digits.
  int shift = -e10;
  while (shift > 0) { io_dd_mul_d(hi, lo, 10.0); --shift; }
  while (shift < 0) { io_dd_mul_d(hi, lo, 0.1);  ++shift; }
  // Drift correction: at most ±1 from the log10 estimate.
  if (hi >= 10.0)     { io_dd_mul_d(hi, lo, 0.1);  ++e10; }
  else if (hi < 1.0)  { io_dd_mul_d(hi, lo, 10.0); --e10; }

  // Extract precision+2 digits; the last two guard round-half-to-even.
  const int ndigits = precision + 2;
  const std::size_t prec = static_cast<std::size_t>(precision);
  char digits[36] = {};
  for (int i = 0; i < ndigits; ++i) {
    int d = static_cast<int>(hi);
    if (d < 0) d = 0;
    else if (d > 9) d = 9;
    digits[i] = static_cast<char>('0' + d);
    hi -= d;
    io_renorm(hi, lo);
    io_dd_mul_d(hi, lo, 10.0);
  }

  int guard = (digits[precision] - '0') * 10 + (digits[precision + 1] - '0');
  bool round_up = guard > 50 ||
                  (guard == 50 && ((digits[precision - 1] - '0') & 1));
  if (round_up) {
    for (int i = precision - 1; i >= 0; --i) {
      if (digits[i] < '9') { ++digits[i]; break; }
      digits[i] = '0';
      if (i == 0) {
        std::memmove(digits + 1, digits, prec - 1);
        digits[0] = '1';
        ++e10;
      }
    }
  }
  digits[precision] = '\0';

  std::string out;
  out.reserve(prec + 8);
  if (neg) out.push_back('-');
  out.push_back(digits[0]);
  if (precision > 1) {
    out.push_back('.');
    out.append(digits + 1, prec - 1);
  }
  out.push_back('e');
  int abs_e = e10 < 0 ? -e10 : e10;
  out.push_back(e10 < 0 ? '-' : '+');
  if (abs_e < 10) out.push_back('0');
  out += std::to_string(abs_e);
  return out;
}

std::ostream &operator<<(std::ostream &os, float64x2 const &x) {
  int p = static_cast<int>(os.precision());
  if (p <= 17) p = 32;
  os << to_string(x, p);
  return os;
}

} // namespace multifloats
