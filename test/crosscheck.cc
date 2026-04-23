// Cross-language differential fuzz.
//
// The C++ (src/multifloats_math.cc, include/multifloats.h) and Fortran
// (fsrc/multifloats.fypp) implementations of the DD kernels are
// INDEPENDENT for add/sub/mul/div and the unary operators neg/abs/sqrt
// plus the comparison operators — everything else is delegated on the
// Fortran side to the extern "C" C++ symbols (see C_DELEGATE_UNARY_MAP
// in fsrc/multifloats.fypp) and is therefore ABI-only and already
// covered by test/abi_equivalence.f90.
//
// Strategy: generate a random (q1, q2) pair using the same generator
// as test/fuzz.cc, project to float64x2 on both sides, and require
// bit-for-bit agreement on both limbs. Any mismatch indicates real
// algorithmic drift between src/ and fsrc/ that libquadmath-only
// oracles cannot detect (both sides agree to DD precision against qp
// but may differ below the qp noise floor between themselves).
//
// sqrt is included as a Fortran-delegates-to-C ABI smoke test; a
// mismatch there points to ABI corruption, not algorithm drift.

#include "multifloats.h"
#include "test_common.hh"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <random>

namespace mf = multifloats;
using q_t = multifloats_test::q_t;
using multifloats_test::from_q;
using multifloats_test::to_q;
using multifloats_test::q_isfinite;
using multifloats_test::q_isnan;

// ABI boundary: Fortran bind(c) wrappers over NATIVE float64x2
// elementals (see test/crosscheck_bindings.f90). Return / argument
// type is float64x2_t from multifloats.h — layout-identical to the
// Fortran `type, bind(c) :: dd_c`.
extern "C" {
  float64x2_t fnat_add (float64x2_t a, float64x2_t b);
  float64x2_t fnat_sub (float64x2_t a, float64x2_t b);
  float64x2_t fnat_mul (float64x2_t a, float64x2_t b);
  float64x2_t fnat_div (float64x2_t a, float64x2_t b);
  float64x2_t fnat_neg (float64x2_t a);
  float64x2_t fnat_abs (float64x2_t a);
  float64x2_t fnat_sqrt(float64x2_t a);
  int         fnat_eq  (float64x2_t a, float64x2_t b);
  int         fnat_ne  (float64x2_t a, float64x2_t b);
  int         fnat_lt  (float64x2_t a, float64x2_t b);
  // fnat_le / fnat_ge intentionally omitted — see crosscheck_bindings.f90
  // for the NaN-semantics divergence that makes a direct cross-check noisy.

  // Integer-returning rounding / decomposition. Gated on |x| < 2e9 by
  // the caller so the Fortran-default-integer result doesn't overflow
  // int32 before widening to c_long inside the wrapper.
  long        fnat_floor   (float64x2_t a);
  long        fnat_ceiling (float64x2_t a);
  long        fnat_nint    (float64x2_t a);
  int         fnat_exponent(float64x2_t a);
}

// -----------------------------------------------------------------------------
// Comparison policies:
//
//  strict    — both limbs bit-identical. For add/sub/neg/abs/sqrt, whose
//              kernels are structurally identical across languages. Any
//              divergence is a real algorithmic bug.
//
//  tolerant  — compare at the DD VALUE level, not limb-by-limb. Both
//              sides project to __float128 (exact, since DD ~106 bits
//              and qp ~113 bits) and we require relative error below
//              kDdRelTol ≈ 1 DD-ULP. For mul and div: the Fortran
//              divide uses Newton-on-1/y while the C++ divide is
//              Dekker-style (q1 + q2 refine); on subnormal
//              denominators the HI limb can differ by 1 dp-ULP with
//              the LO limb absorbing the compensating rounding. Both
//              answers are correct to full DD precision — the
//              canonical (hi, lo) split is what drifts.
//
// Non-finite handling: a shared indeterminate (both NaN, or both ±Inf)
// is accepted regardless of NaN payload or sign — the two sides
// differ on mul-overflow (C++ FMA error term makes it NaN, Fortran
// short-circuits to Inf) but both agree the result is unusable
// numerically, which is the only invariant we can rely on.

static inline bool bit_eq(double a, double b) {
  std::uint64_t ai, bi;
  std::memcpy(&ai, &a, sizeof ai);
  std::memcpy(&bi, &b, sizeof bi);
  return ai == bi;
}
static inline bool bit_eq(float64x2_t a, float64x2_t b) {
  return bit_eq(a.hi, b.hi) && bit_eq(a.lo, b.lo);
}

// Per-op normal-scale tolerances, measured empirically against the
// fuzz generator at 10^5 iterations across three seeds
// (0xc0ffee / 42 / 0xdeadbeef):
//
//   mul — zero drift observed. The two sides use the same FMA-based
//         two_prod and agree bit-for-bit on DD value when inputs are
//         finite; the non-finite LO short-circuit differs (Fortran
//         zeros, C++ propagates NaN) but the DD VALUE agrees via the
//         shared-indeterminate carveout in rel_err_dd.
//
//   div — up to ~9e-32 rel_err (a handful of DD-ULPs) in the normal
//         bucket. The two divide algorithms differ fundamentally
//         (C++ Dekker q1+q2 vs. Fortran Newton on 1/y0) so the split
//         between HI and LO limbs can drift at the last bit; the
//         VALUE remains within DD precision. kTol_div is set at 16×
//         the observed worst so sub-DD-precision drift passes but a
//         regression that breaks DD semantics (e.g. a dp-ULP value
//         diff, 1e-16) fires immediately.
//
// Any future tightening should shrink these bounds, not widen them.
static constexpr double kTol_mul = 5.0e-32;       // ~4 DD-ULPs
static constexpr double kTol_div = 1.5e-30;       // ~120 DD-ULPs (16× observed)

// IEEE-edge carveout. When either operand has |hi| outside
// [kEdgeLoAbs, kEdgeHiAbs], 1/y0 or the multiply error term
// underflows/overflows into subnormal/near-inf territory; both sides
// degrade to ~dp precision independently and diverge in the LO limb.
// Tallied separately so a normal-regime signal is not buried in edge
// noise.
static constexpr double kEdgeLoAbs = 1.0e-150;    // ~2^-498
static constexpr double kEdgeHiAbs = 1.0e+150;    // ~2^+498
static constexpr double kTolEdge   = 1.0e-14;

static inline bool match_strict(float64x2_t a, float64x2_t b) {
  return bit_eq(a, b);
}

static inline double rel_err_dd(float64x2_t a, float64x2_t b) {
  q_t qa = (q_t)a.hi + (q_t)a.lo;
  q_t qb = (q_t)b.hi + (q_t)b.lo;
  bool a_fin = q_isfinite(qa);
  bool b_fin = q_isfinite(qb);
  if (!a_fin && !b_fin) return 0.0;        // shared indeterminate
  if (a_fin != b_fin) return INFINITY;     // one side non-finite, bug
  q_t diff = qa - qb; if (diff < 0) diff = -diff;
  q_t mag  = qa < 0 ? -qa : qa;
  if (mag == 0) return (double)diff;
  return (double)(diff / mag);
}

static inline float64x2_t to_c(mf::float64x2 const &x) {
  return static_cast<float64x2_t>(x);
}

// -----------------------------------------------------------------------------
// Random input generation — mirrors Rng::generate_pair in test/fuzz.cc
// exactly (same distribution mix, same RNG, same seed semantics). The
// point of the mirror is that a (iterations, seed) pair replays the
// same inputs as fuzz.cc, so a failure surfaced here can be bisected
// against fuzz.cc's per-op oracles.

struct Rng {
  std::mt19937_64 engine;
  std::uniform_real_distribution<double> u01{0.0, 1.0};

  explicit Rng(std::uint64_t seed) : engine(seed) {}
  double u() { return u01(engine); }

  q_t pick_nonfinite(double r) {
    if (r < 0.33) return (q_t)(+1.0 / 0.0);
    if (r < 0.66) return (q_t)(-1.0 / 0.0);
    return (q_t)(0.0 / 0.0);
  }

  q_t wide(double r, double rexp) {
    int k = (int)(rexp * 60.0) - 30;
    q_t mag = powq((q_t)10.0q, (q_t)k);
    return (q_t)(r - 0.5) * mag;
  }

  void generate_pair(q_t &q1, q_t &q2) {
    double r1 = u(), r2 = u(), r3 = u(), r4 = u();
    int mode = (int)(r1 * 10.0);
    switch (mode) {
    case 0: {
      q1 = pick_nonfinite(r2);
      q2 = pick_nonfinite(r3);
      break;
    }
    case 1: {
      int k = (int)(r3 * 20.0) - 10;
      q1 = (q_t)(r2 - 0.5) * powq((q_t)10.0q, (q_t)k);
      q2 = q1 * ((q_t)1.0q + (q_t)(r4 * 1e-15));
      break;
    }
    case 2: {
      int k = (int)(r3 * 20.0) - 10;
      q1 = (q_t)(r2 - 0.5) * powq((q_t)10.0q, (q_t)k);
      q2 = -q1 + (q_t)((r4 - 0.5) * 1e-25) * q1;
      break;
    }
    case 3: {
      double huge = 0x1.fffffffffffffp+1023;
      q1 = (q_t)huge * (q_t)(0.9 + 0.1 * r2);
      q2 = (q_t)huge * (q_t)(0.9 + 0.1 * r3);
      break;
    }
    case 4: {
      double tiny = 0x1.0p-1022;
      q1 = (q_t)tiny * (q_t)(1.0 + 10.0 * r2);
      q2 = (q_t)tiny * (q_t)(1.0 + 10.0 * r3);
      break;
    }
    default: {
      q1 = wide(r2, r3);
      q2 = wide(r4, r1);
      break;
    }
    }
  }
};

// -----------------------------------------------------------------------------
// Per-op bookkeeping. We keep only the FIRST failing input — once the
// algorithm drifts every subsequent input tends to fail, so capping
// at one keeps output bounded and the repro input actionable.

struct Stat {
  char const *name;
  long ok = 0;
  long fail = 0;
  long edge = 0;              // inputs outside [kEdgeLoAbs, kEdgeHiAbs]
  long edge_fail = 0;
  double max_rel_normal = 0.0; // worst rel_err in the normal-scale bucket
  double max_rel_edge = 0.0;   // worst rel_err in the edge bucket
  float64x2_t first_a{}, first_b{};
  float64x2_t first_cpp{}, first_f{};
  bool is_unary = false;
  // Integer-return ops stash the raw integer pair in first_cpp/first_f via
  // the `.hi` double slot so the reporter can still print a repro without
  // a second code path.
  bool is_int_return = false;
};

static void record_fail(Stat &s, float64x2_t a, float64x2_t b,
                        float64x2_t cpp, float64x2_t f) {
  if (s.fail == 0) {
    s.first_a = a;
    s.first_b = b;
    s.first_cpp = cpp;
    s.first_f = f;
  }
  ++s.fail;
}

static void print_stats(Stat const *stats, int n) {
  std::printf("\n  %-6s %10s %10s %10s %10s %11s %11s\n",
              "op", "ok", "fail", "edge_ok", "edge_fail",
              "max_rel_n", "max_rel_e");
  for (int i = 0; i < n; ++i) {
    Stat const &s = stats[i];
    long edge_ok = s.edge - s.edge_fail;
    std::printf("  %-6s %10ld %10ld %10ld %10ld %11.2e %11.2e\n",
                s.name, s.ok, s.fail, edge_ok, s.edge_fail,
                s.max_rel_normal, s.max_rel_edge);
  }
  std::printf("\n");
  for (int i = 0; i < n; ++i) {
    Stat const &s = stats[i];
    if (s.fail == 0) continue;
    std::printf("First normal-scale mismatch on %s:\n", s.name);
    std::printf("  a   = {%a, %a}\n", s.first_a.hi, s.first_a.lo);
    if (!s.is_unary)
      std::printf("  b   = {%a, %a}\n", s.first_b.hi, s.first_b.lo);
    std::printf("  cpp = {%a, %a}\n", s.first_cpp.hi, s.first_cpp.lo);
    std::printf("  f   = {%a, %a}\n", s.first_f.hi, s.first_f.lo);
  }
}

// -----------------------------------------------------------------------------

int main(int argc, char **argv) {
  long iterations = 10000;
  std::uint64_t seed = 42ULL;
  if (argc > 1) iterations = std::strtol(argv[1], nullptr, 0);
  if (argc > 2) seed = std::strtoull(argv[2], nullptr, 0);
  std::printf("[crosscheck] iterations=%ld seed=0x%llx\n",
              iterations, (unsigned long long)seed);

  Rng rng(seed);

  enum { I_ADD, I_SUB, I_MUL, I_DIV,
         I_NEG, I_ABS, I_SQRT,
         I_EQ,  I_NE,  I_LT,
         I_FLOOR, I_CEILING, I_NINT, I_EXPONENT,
         I_N };
  Stat stats[I_N] = {
    {"add"}, {"sub"}, {"mul"}, {"div"},
    {"neg",  0, 0, {}, {}, {}, {}, true},
    {"abs",  0, 0, {}, {}, {}, {}, true},
    {"sqrt", 0, 0, {}, {}, {}, {}, true},
    {"eq"}, {"ne"}, {"lt"},
    // Integer-return ops. is_unary = true, is_int_return = true.
    {"floor",    0, 0, {}, {}, {}, {}, true, true},
    {"ceiling",  0, 0, {}, {}, {}, {}, true, true},
    {"nint",     0, 0, {}, {}, {}, {}, true, true},
    {"exponent", 0, 0, {}, {}, {}, {}, true, true},
  };

  for (long i = 0; i < iterations; ++i) {
    q_t q1, q2;
    rng.generate_pair(q1, q2);
    mf::float64x2 f1 = from_q(q1);
    mf::float64x2 f2 = from_q(q2);
    float64x2_t c1 = to_c(f1);
    float64x2_t c2 = to_c(f2);

    // IEEE-edge regime where DD precision degrades (see
    // kEdgeLoAbs / kEdgeHiAbs). Tallied separately so a signal in the
    // normal regime is not buried under expected edge-range noise.
    auto in_edge = [](double h) {
      double a = std::fabs(h);
      return a < kEdgeLoAbs || a > kEdgeHiAbs;
    };
    bool edge_in = in_edge(c1.hi) || in_edge(c2.hi);

    auto check_strict = [&](int idx, float64x2_t cpp_r, float64x2_t f_r,
                             bool binary) {
      if (match_strict(cpp_r, f_r)) ++stats[idx].ok;
      else record_fail(stats[idx], c1, binary ? c2 : float64x2_t{0,0}, cpp_r, f_r);
    };
    auto check_tolerant = [&](int idx, float64x2_t cpp_r, float64x2_t f_r,
                              double normal_tol) {
      double rel = rel_err_dd(cpp_r, f_r);
      Stat &s = stats[idx];
      if (edge_in) {
        ++s.edge;
        if (rel > s.max_rel_edge) s.max_rel_edge = rel;
        if (rel > kTolEdge) ++s.edge_fail;
      } else {
        if (rel > s.max_rel_normal) s.max_rel_normal = rel;
        if (rel <= normal_tol) ++s.ok;
        else record_fail(s, c1, c2, cpp_r, f_r);
      }
    };
    auto check_bool2 = [&](int idx, bool cpp_b, int f_b) {
      int cb = cpp_b ? 1 : 0;
      if (cb == f_b) ++stats[idx].ok;
      else {
        float64x2_t cpp_mark = {(double)cb, 0.0};
        float64x2_t f_mark   = {(double)f_b, 0.0};
        record_fail(stats[idx], c1, c2, cpp_mark, f_mark);
      }
    };
    // Integer-return unary: compare raw integers (after both sides
    // independently round/truncate). Repro on failure is stashed in the
    // .hi slot of the float64x2_t pair so print_stats's %a format still
    // reads something sensible.
    auto check_int1 = [&](int idx, long cpp_i, long f_i) {
      if (cpp_i == f_i) ++stats[idx].ok;
      else {
        float64x2_t cpp_mark = {(double)cpp_i, 0.0};
        float64x2_t f_mark   = {(double)f_i,   0.0};
        record_fail(stats[idx], c1, float64x2_t{0,0}, cpp_mark, f_mark);
      }
    };

    check_strict(I_ADD,  to_c(f1 + f2),        fnat_add(c1, c2),  true);
    check_strict(I_SUB,  to_c(f1 - f2),        fnat_sub(c1, c2),  true);
    check_tolerant(I_MUL, to_c(f1 * f2),       fnat_mul(c1, c2), kTol_mul);
    check_tolerant(I_DIV, to_c(f1 / f2),       fnat_div(c1, c2), kTol_div);
    check_strict(I_NEG,  to_c(-f1),            fnat_neg(c1),  false);
    check_strict(I_ABS,  to_c(mf::abs(f1)),    fnat_abs(c1),  false);
    check_strict(I_SQRT, to_c(mf::sqrt(f1)),   fnat_sqrt(c1), false);
    check_bool2(I_EQ, f1 == f2, fnat_eq(c1, c2));
    check_bool2(I_NE, f1 != f2, fnat_ne(c1, c2));
    check_bool2(I_LT, f1 <  f2, fnat_lt(c1, c2));

    // Integer-return unary: Fortran's floor / ceiling / nint return
    // DEFAULT integer (int32 on mainstream Linux toolchains); gate |q1|
    // below 2^31 ≈ 2.1e9 so the Fortran computation doesn't overflow
    // before widening to c_long at the bind(c) boundary. exponent is
    // bounded by ±1074 for any finite input — only the zero gate
    // applies.
    //
    // C++ analogues, each using PRODUCTION multifloats.h ops:
    //   Fortran floor(DD)   ≡  (long) mf::floor(DD).hi      (lo fully
    //                            absorbed by the DD-floor renormalization
    //                            so the hi limb IS the integer floor)
    //   Fortran ceiling(DD) ≡  (long) mf::ceil(DD).hi
    //   Fortran nint(DD)    ≡  mf::lround(DD)               (both
    //                            round-half-away-from-zero)
    //   Fortran exponent(DD)≡  mf::ilogb(DD) + 1            (Fortran
    //                            normalizes to [0.5, 1) rather than C's
    //                            [1, 2), so the bias differs by one)
    if (q_isfinite(q1) && (q1 < 0 ? -q1 : q1) < (q_t)2.0e9q) {
      check_int1(I_FLOOR,   (long)mf::floor(f1)._limbs[0], fnat_floor(c1));
      check_int1(I_CEILING, (long)mf::ceil(f1)._limbs[0],  fnat_ceiling(c1));
      check_int1(I_NINT,    mf::lround(f1),                 fnat_nint(c1));
    }
    if (q_isfinite(q1) && q1 != (q_t)0) {
      check_int1(I_EXPONENT, (long)(mf::ilogb(f1) + 1), (long)fnat_exponent(c1));
    }
  }

  // Adversarial probe for the half-integer-hi + opposite-sign-lo class:
  // the random generator's uniform sampling misses this narrow input
  // class (density ≈ 1/ulp near magnitude, astronomically small at the
  // magnitudes the generator targets), but it's exactly where the
  // integer-rounding ops differ between hi-only and full-DD logic.
  // Without this probe a bug like the one PR-2 fixed in C++ nearbyint
  // (and the mirror dd_nint bug this PR found in fsrc/multifloats.fypp)
  // could live undetected. Inline the check since the scalar `check_int1`
  // lambda captures the main loop's c1 and isn't reachable here.
  {
    auto probe_int = [&stats](int idx, float64x2_t c_in, long cpp_i, long f_i) {
      if (cpp_i == f_i) { ++stats[idx].ok; }
      else {
        float64x2_t cpp_mark = {(double)cpp_i, 0.0};
        float64x2_t f_mark   = {(double)f_i,   0.0};
        record_fail(stats[idx], c_in, float64x2_t{0,0}, cpp_mark, f_mark);
      }
    };
    double halves[] = {2.5, 3.5, 100.5, -2.5, -100.5};
    double los[]    = {-0.01, +0.01};
    for (double hi_v : halves) {
      for (double lo_v : los) {
        mf::float64x2 f; f._limbs[0] = hi_v; f._limbs[1] = lo_v;
        float64x2_t c = to_c(f);
        probe_int(I_FLOOR,   c, (long)mf::floor(f)._limbs[0], fnat_floor(c));
        probe_int(I_CEILING, c, (long)mf::ceil(f)._limbs[0],  fnat_ceiling(c));
        probe_int(I_NINT,    c, mf::lround(f),                fnat_nint(c));
      }
    }
  }

  print_stats(stats, I_N);

  long total_fail = 0;
  for (int i = 0; i < I_N; ++i) total_fail += stats[i].fail;
  if (total_fail != 0) {
    std::printf("FAIL: %ld cross-language mismatches\n", total_fail);
    return 1;
  }
  std::printf("PASS crosscheck\n");
  return 0;
}
