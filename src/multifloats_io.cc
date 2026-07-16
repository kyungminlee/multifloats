// Decimal I/O for float64x2 — scientific-notation formatting.
//
// Layering. The core formatter `format_scientific_chars` is a file-local
// (anon-namespace) helper, called directly by three siblings in this TU:
//   * `multifloats::to_chars` (primary 5-arg overload) — wraps it with the
//     `std::to_chars_result` return convention and the chars_format check.
//   * `multifloats::operator<<` — writes the formatted bytes to an ostream.
//   * extern "C" `to_charsdd` — the shim for C / Fortran consumers
//     (declared in the C-ABI section of multifloats.h).
//
// Header-inline `to_string` (in multifloats.h) sits above `to_chars` and
// constructs the returned `std::string` entirely in the consumer's TU,
// keeping the library's link surface free of `std::string` and the
// libstdc++ `_GLIBCXX_USE_CXX11_ABI` dual-ABI mangling distinction.
// `std::to_chars_result` is ABI-stable (plain `{char*, errc}`), so the
// 5-arg `to_chars` is safe to define out-of-line here.

#include "multifloats.h"

#include <charconv>
#include <cmath>
#include <cstddef>
#include <cstdio>   // std::snprintf for exponent formatting
#include <cstring>
#include <ostream>
#include <system_error>

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

// Scientific-notation decimal formatter for float64x2. Writes into
// [first, last) without a NUL terminator, returns one-past-last-written
// on success or nullptr when the output wouldn't fit (in which case the
// buffer is left unchanged). `precision` is clamped to [1, 34]. Special
// values emit "nan", "inf", "-inf", "0e+00", "-0e+00".
//
// Strategy: produce the full output into a fixed 48-byte tmp[] buffer
// (sized to hold any legal output), then size-check once and memcpy into
// the caller's range. Any legal output fits in at most ~42 bytes.
char *format_scientific_chars(multifloats::float64x2 const &value,
                              int precision,
                              char *first, char *last) noexcept {
  double hi = value.limbs[0];
  double lo = value.limbs[1];

  char tmp[48];
  std::size_t w = 0;

  auto put_char = [&](char c) {
    if (w < sizeof(tmp)) tmp[w++] = c;
  };
  auto put_str = [&](const char *s) {
    while (*s && w < sizeof(tmp)) tmp[w++] = *s++;
  };

  if (std::isnan(hi) || std::isnan(lo)) {
    put_str("nan");
  } else if (std::isinf(hi)) {
    put_str(hi > 0 ? "inf" : "-inf");
  } else if (hi == 0.0) {
    put_str(std::signbit(hi) ? "-0e+00" : "0e+00");
  } else {
    if (precision < 1) precision = 1;
    if (precision > 34) precision = 34;

    const bool neg = hi < 0.0;
    if (neg) { hi = -hi; lo = -lo; }

    // Decimal exponent estimate. log10 may be off by 1 due to rounding; a
    // post-scale fixup corrects.
    int e10 = static_cast<int>(std::floor(std::log10(hi)));

    // Scale to [1, 10) by a power of ten held to full DD precision. For each
    // set bit k of |e10| we divide (e10 > 0) or multiply (e10 < 0) by the DD
    // constant 10^(2^k) — O(log|e10|) DD ops at ~1 ulp each. The previous code
    // iterated ×0.1, injecting ~1 ulp PER decimal place, which destroyed most
    // of the 32 significant digits once |x| ≥ 10. Scaling by one table entry at
    // a time (not forming 10^|e10| first) keeps the intermediate from
    // overflowing for tiny subnormal inputs. Constants are the correctly-
    // rounded DD of 10^(2^k), k=0..8 (10^1 … 10^256; |e10| ≤ ~323 < 512).
    static const double kP10Hi[9] = {
        0x1.4p+3, 0x1.9p+6, 0x1.388p+13, 0x1.7d784p+26, 0x1.1c37937e08p+53,
        0x1.3b8b5b5056e17p+106, 0x1.84f03e93ff9f5p+212, 0x1.27748f9301d32p+425,
        0x1.54fdd7f73bf3cp+850};
    static const double kP10Lo[9] = {
        0.0, 0.0, 0.0, 0.0, 0.0, -0x1.3107fp+52, -0x1.2ac340948e389p+157,
        -0x1.901cc86649e4ap+371, -0x1.7222446fe467p+795};
    {
      multifloats::float64x2 v(hi, lo);
      for (int k = 0, n = (e10 >= 0 ? e10 : -e10); n && k < 9; ++k, n >>= 1) {
        if (n & 1) {
          multifloats::float64x2 p(kP10Hi[k], kP10Lo[k]);
          v = (e10 >= 0) ? v / p : v * p;
        }
      }
      // Drift correction: the log10 estimate can be ±1.
      multifloats::float64x2 ten(kP10Hi[0], 0.0);
      if (v.limbs[0] >= 10.0)     { v = v / ten; ++e10; }
      else if (v.limbs[0] < 1.0)  { v = v * ten; --e10; }
      hi = v.limbs[0];
      lo = v.limbs[1];
    }

    // Extract precision+2 digits; the last two guard round-half-to-even.
    const int ndigits = precision + 2;
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
          std::memmove(digits + 1, digits,
                       static_cast<std::size_t>(precision) - 1);
          digits[0] = '1';
          ++e10;
        }
      }
    }
    digits[precision] = '\0';

    if (neg) put_char('-');
    put_char(digits[0]);
    if (precision > 1) {
      put_char('.');
      for (int i = 1; i < precision; ++i) put_char(digits[i]);
    }
    put_char('e');
    int abs_e = e10 < 0 ? -e10 : e10;
    put_char(e10 < 0 ? '-' : '+');
    // Two-digit minimum exponent (matches libc `%e` formatting).
    char ebuf[8];
    int ew = std::snprintf(ebuf, sizeof(ebuf), "%d", abs_e);
    if (ew == 1) put_char('0');
    for (int i = 0; i < ew; ++i) put_char(ebuf[i]);
  }

  std::size_t avail = (last > first) ? static_cast<std::size_t>(last - first) : 0;
  if (w > avail) return nullptr;
  std::memcpy(first, tmp, w);
  return first + w;
}

} // anonymous namespace

namespace multifloats {

std::to_chars_result
to_chars(char *first, char *last, float64x2 const &value,
         std::chars_format fmt, int precision) noexcept {
  if (fmt != std::chars_format::scientific)
    return {first, std::errc::invalid_argument};
  char *end = format_scientific_chars(value, precision, first, last);
  if (!end) return {last, std::errc::value_too_large};
  return {end, std::errc{}};
}

std::ostream &operator<<(std::ostream &os, float64x2 const &x) {
  int p = static_cast<int>(os.precision());
  if (p <= 17) p = 32;
  char buf[MULTIFLOATS_DD_CHARS_BUFSIZE];
  char *end = format_scientific_chars(x, p, buf, buf + sizeof(buf));
  // `buf` is sized to hold any legal output, so `end` is non-null here.
  if (end) os.write(buf, end - buf);
  return os;
}

// The C-ABI `to_charsdd` shim lives inside `namespace multifloats` so that
// its C++ qualified name matches the declaration in the public header.
// `extern "C"` keeps the linker symbol unmangled — C and Fortran callers
// see plain `to_charsdd` just the same.
extern "C" {
char *to_charsdd(float64x2 x, int precision, char *first, char *last) {
  return format_scientific_chars(x, precision, first, last);
}
// C23-style forwarder — portable across ELF/Mach-O. See note in
// float64x2/f64x2_aliases.inc on why we use wrappers instead of
// __attribute__((alias)).
char *to_charsf64x2(float64x2 x, int precision, char *first, char *last) {
  return to_charsdd(x, precision, first, last);
}
} // extern "C"

} // namespace multifloats
