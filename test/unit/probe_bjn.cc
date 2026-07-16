// Regime probe for jn(n, x): deterministic sweep comparing multifloats
// jn vs boost::math::cyl_bessel_j against the __float128 oracle. Bins
// max_rel by regime so we can see whether the 19× gap reported in the
// fuzz report is forward-recurrence (n <= x), Miller (n > x), near-root,
// or large-x.
//
// Build with -DMULTIFLOATS_BUILD_BOOST_COMPARE=ON; target name `bjn_probe`.

#include "test_common.hh"

#include <boost/multiprecision/cpp_double_fp.hpp>
#include <boost/math/special_functions/bessel.hpp>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>

#include <quadmath.h>

using boost::multiprecision::cpp_double_double;
using multifloats_test::q_t;
using multifloats_test::to_q;
using multifloats_test::from_q;
using multifloats_test::qstr;

static inline q_t bdd_to_q(cpp_double_double const &x) {
  auto const &r = x.backend().crep();
  return (q_t)r.first + (q_t)r.second;
}

static double q_rel_err(q_t got, q_t ref) {
  if (isnanq(got) || isnanq(ref)) return 0.0;
  q_t denom = ref;
  if (denom < 0) denom = -denom;
  if (denom == 0) denom = 1;
  q_t e = got - ref;
  if (e < 0) e = -e;
  return (double)(e / denom);
}

struct Bin {
  char const *label;
  long n = 0;
  double max_rel_mf = 0;
  double max_rel_bd = 0;
  double sum_rel_mf = 0;
  double sum_rel_bd = 0;
  double worst_x = 0;
  int worst_n = 0;
};

static void update(Bin &b, double rel_mf, double rel_bd, double x, int n) {
  ++b.n;
  b.sum_rel_mf += rel_mf;
  b.sum_rel_bd += rel_bd;
  if (rel_mf > b.max_rel_mf) {
    b.max_rel_mf = rel_mf;
    b.worst_x = x;
    b.worst_n = n;
  }
  if (rel_bd > b.max_rel_bd) b.max_rel_bd = rel_bd;
}

int main() {
  // Regime bins
  Bin forward_far    {"forward (n <= x/2)"};
  Bin forward_near   {"forward (x/2 < n <= x)"};
  Bin miller_near    {"miller (x < n <= 2x)"};
  Bin miller_far     {"miller (n > 2x)"};
  Bin near_root      {"|J_n(x)| < 1e-3 (any regime)"};

  static constexpr int kOrders[] = {2, 3, 5, 8};

  // Sweep x with linear + log spacing on (0, 200].
  // 1) Dense linear in [0.1, 50]: 5000 points
  // 2) Log in [50, 200]: 1000 points
  // 3) Targeted near-root: scan x near each j_n root via 100 offsets
  //    at +- 1e-3 / 1e-5 / 1e-8 around computed root via boost.

  auto probe_x = [&](double x) {
    if (!(x > 0)) return;
    q_t qx = (q_t)x;
    cpp_double_double bdx;
    {
      auto &r = bdx.backend().rep();
      r.first = x; r.second = 0.0;
    }
    multifloats::float64x2 mfx;
    mfx.limbs[0] = x; mfx.limbs[1] = 0.0;

    for (int n : kOrders) {
      q_t qref = jnq(n, qx);
      if (!finiteq(qref)) continue;
      double aref = (double)(qref < 0 ? -qref : qref);

      multifloats::float64x2 mfr = multifloats::jn(n, mfx);
      cpp_double_double bdr = boost::math::cyl_bessel_j(n, bdx);
      double rmf = q_rel_err(to_q(mfr), qref);
      double rbd = q_rel_err(bdd_to_q(bdr), qref);

      Bin *b;
      if (n <= x * 0.5)      b = &forward_far;
      else if (n <= x)       b = &forward_near;
      else if (n <= 2.0 * x) b = &miller_near;
      else                   b = &miller_far;
      update(*b, rmf, rbd, x, n);

      if (aref < 1e-3) update(near_root, rmf, rbd, x, n);
    }
  };

  for (int i = 0; i < 5000; ++i) {
    double x = 0.1 + (50.0 - 0.1) * (double)i / 4999.0;
    probe_x(x);
  }
  for (int i = 0; i < 1000; ++i) {
    double t = (double)i / 999.0;
    double x = std::pow(10.0, std::log10(50.0) + t * (std::log10(200.0) - std::log10(50.0)));
    probe_x(x);
  }

  auto report = [](Bin const &b) {
    if (b.n == 0) { std::printf("  %-32s   (empty)\n", b.label); return; }
    double mean_mf = b.sum_rel_mf / (double)b.n;
    double mean_bd = b.sum_rel_bd / (double)b.n;
    double ratio = b.max_rel_mf > 0 && b.max_rel_bd > 0
                   ? b.max_rel_mf / b.max_rel_bd : 0.0;
    std::printf("  %-32s n=%-6ld  max mf=%.2e  max bd=%.2e  mf/bd=%5.1fx  "
                "mean mf=%.2e mean bd=%.2e  worst x=%.4f n=%d\n",
                b.label, b.n, b.max_rel_mf, b.max_rel_bd, ratio,
                mean_mf, mean_bd, b.worst_x, b.worst_n);
  };

  std::printf("[bjn_probe] DETERMINISTIC SWEEP — regime breakdown for jn\n");
  report(forward_far);
  report(forward_near);
  report(miller_near);
  report(miller_far);
  report(near_root);

  // ===== Fuzz-style random sweep, mirrors fuzz.cc generate_pair distribution
  // =====================================================================
  std::mt19937_64 eng(42);
  std::uniform_real_distribution<double> u01(0.0, 1.0);

  Bin fz_forward_far  {"forward (n <= x/2)"};
  Bin fz_forward_near {"forward (x/2 < n <= x)"};
  Bin fz_miller_near  {"miller (x < n <= 2x)"};
  Bin fz_miller_far   {"miller (n > 2x)"};
  Bin fz_near_root    {"|J_n(x)| < 1e-3 (any regime)"};

  // Worst-case capture
  struct Worst { double x = 0; int n = 0; double rel_mf = 0; double rel_bd = 0;
                 q_t qref = 0; q_t qmf = 0; q_t qbd = 0; };
  Worst worst;

  long total = 0, kept = 0;
  long iters = 1000000;
  for (long it = 0; it < iters; ++it) {
    double r1 = u01(eng), r2 = u01(eng), r3 = u01(eng), r4 = u01(eng);
    int mode = (int)(r1 * 10.0);
    q_t q1;
    if (mode == 0) {
      // skip non-finite for this probe
      continue;
    } else if (mode == 1 || mode == 2) {
      int k = (int)(r3 * 20.0) - 10;
      q1 = (q_t)(r2 - 0.5) * powq((q_t)10.0q, (q_t)k);
    } else if (mode == 3) {
      int k = (int)(r3 * 10.0) + 20;
      q1 = (q_t)(r2 - 0.5) * 2.0q * powq((q_t)10.0q, (q_t)k);
    } else if (mode == 4) {
      int k = -((int)(r3 * 10.0) + 20);
      q1 = (q_t)(r2 - 0.5) * 2.0q * powq((q_t)10.0q, (q_t)k);
    } else {
      int k = (int)(r3 * 60.0) - 30;
      q1 = (q_t)(r2 - 0.5) * powq((q_t)10.0q, (q_t)k);
    }
    if (!finiteq(q1)) continue;
    q_t aq = q1 < 0 ? -q1 : q1;
    if (!(aq < (q_t)200)) continue;

    double x = (double)q1;
    double ax = std::fabs(x);
    if (!(ax > 0)) continue;

    cpp_double_double bdx; { auto &r = bdx.backend().rep(); r.first=x; r.second=(double)(q1 - (q_t)x); }
    multifloats::float64x2 mfx; mfx.limbs[0]=x; mfx.limbs[1]=(double)(q1 - (q_t)x);

    ++total;
    for (int n : kOrders) {
      q_t qref = jnq(n, q1);
      if (!finiteq(qref)) continue;
      double aref = (double)(qref < 0 ? -qref : qref);
      multifloats::float64x2 mfr = multifloats::jn(n, mfx);
      cpp_double_double bdr = boost::math::cyl_bessel_j(n, bdx);
      q_t qmf = to_q(mfr), qbd = bdd_to_q(bdr);
      double rmf = q_rel_err(qmf, qref);
      double rbd = q_rel_err(qbd, qref);

      Bin *b;
      if ((double)n <= ax * 0.5)      b = &fz_forward_far;
      else if ((double)n <= ax)       b = &fz_forward_near;
      else if ((double)n <= 2.0 * ax) b = &fz_miller_near;
      else                            b = &fz_miller_far;
      update(*b, rmf, rbd, x, n);
      if (aref < 1e-3) update(fz_near_root, rmf, rbd, x, n);
      ++kept;

      if (rmf > worst.rel_mf) {
        worst.rel_mf = rmf; worst.rel_bd = rbd;
        worst.x = x; worst.n = n;
        worst.qref = qref; worst.qmf = qmf; worst.qbd = qbd;
      }
    }
  }

  std::printf("\n[bjn_probe] FUZZ-STYLE SWEEP (seed=42, iters=%ld, draws=%ld, "
              "samples=%ld)\n", iters, total, kept);
  report(fz_forward_far);
  report(fz_forward_near);
  report(fz_miller_near);
  report(fz_miller_far);
  report(fz_near_root);

  std::printf("\nworst multifloats jn case:\n");
  std::printf("  n=%d  x=%.17g (= %a)\n", worst.n, worst.x, worst.x);
  std::printf("  qref = %s\n", qstr(worst.qref));
  std::printf("  mf   = %s   rel=%.3e\n", qstr(worst.qmf), worst.rel_mf);
  std::printf("  bd   = %s   rel=%.3e\n", qstr(worst.qbd), worst.rel_bd);
  return 0;
}
