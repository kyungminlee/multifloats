// A/B precision + perf probe: current 4-mul vs Karatsuba 3-mul for DD complex
// multiply. Run against a __float128 reference.
#include "../src/multifloats.hh"
#include <quadmath.h>
#include <cstdio>
#include <chrono>
#include <random>
namespace mf = multifloats;
using q = __float128;

static q to_q(mf::float64x2 const& x) { return (q)x._limbs[0] + (q)x._limbs[1]; }
static mf::float64x2 from_q(q v) {
  mf::float64x2 r;
  r._limbs[0] = (double)v;
  r._limbs[1] = (double)(v - (q)r._limbs[0]);
  return r;
}

struct cdd { mf::float64x2 re, im; };

static cdd mul_4(cdd const& a, cdd const& b) {
  return { a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re };
}
static cdd mul_k(cdd const& a, cdd const& b) {
  mf::float64x2 p = a.re * b.re;
  mf::float64x2 q = a.im * b.im;
  mf::float64x2 r = (a.re + a.im) * (b.re + b.im);
  return { p - q, r - p - q };
}

int main() {
  std::mt19937_64 rng(42);
  std::uniform_real_distribution<double> uni(-1.0, 1.0);

  int const N = 200000;
  double max_rel_4_re = 0, max_rel_4_im = 0;
  double max_rel_k_re = 0, max_rel_k_im = 0;
  for (int i = 0; i < N; ++i) {
    cdd a { from_q((q)uni(rng) + (q)uni(rng) * ldexpq(1, -52)),
            from_q((q)uni(rng) + (q)uni(rng) * ldexpq(1, -52)) };
    cdd b { from_q((q)uni(rng) + (q)uni(rng) * ldexpq(1, -52)),
            from_q((q)uni(rng) + (q)uni(rng) * ldexpq(1, -52)) };
    q ar=to_q(a.re), ai=to_q(a.im), br=to_q(b.re), bi=to_q(b.im);
    q ref_re = ar*br - ai*bi;
    q ref_im = ar*bi + ai*br;
    q mag = fabsq(ref_re) + fabsq(ref_im); if (mag < 1e-100q) continue;

    cdd g4 = mul_4(a, b);
    double r4re = (double)(fabsq(to_q(g4.re) - ref_re) / mag);
    double r4im = (double)(fabsq(to_q(g4.im) - ref_im) / mag);
    if (r4re > max_rel_4_re) max_rel_4_re = r4re;
    if (r4im > max_rel_4_im) max_rel_4_im = r4im;

    cdd gk = mul_k(a, b);
    double rkre = (double)(fabsq(to_q(gk.re) - ref_re) / mag);
    double rkim = (double)(fabsq(to_q(gk.im) - ref_im) / mag);
    if (rkre > max_rel_k_re) max_rel_k_re = rkre;
    if (rkim > max_rel_k_im) max_rel_k_im = rkim;
  }
  printf("precision (max_rel over %d random points):\n", N);
  printf("  4-mul   re %.3e  im %.3e\n", max_rel_4_re, max_rel_4_im);
  printf("  karats  re %.3e  im %.3e\n", max_rel_k_re, max_rel_k_im);

  // Cancellation stress test: a = (1, eps), b = (-1, eps). Im should be 0.
  double eps = 1e-18;
  cdd a = { mf::float64x2(1.0), from_q((q)eps) };
  cdd b = { mf::float64x2(-1.0), from_q((q)eps) };
  cdd g4 = mul_4(a, b);
  cdd gk = mul_k(a, b);
  printf("\ncancellation stress (true Im = 0):\n");
  printf("  4-mul   im = %.3e  (hi=%g lo=%g)\n",
         (double)to_q(g4.im), g4.im._limbs[0], g4.im._limbs[1]);
  printf("  karats  im = %.3e  (hi=%g lo=%g)\n",
         (double)to_q(gk.im), gk.im._limbs[0], gk.im._limbs[1]);

  // Perf: loop-carried dependency so the compiler can't hoist.
  int const M = 10'000'000;
  cdd a0 { mf::float64x2(0.7), mf::float64x2(0.3) };
  cdd b0 { mf::float64x2(0.9), mf::float64x2(-0.2) };
  volatile double sink = 0;

  auto t1 = std::chrono::steady_clock::now();
  cdd acc = a0;
  for (int i = 0; i < M; ++i) acc = mul_4(acc, b0);
  auto t2 = std::chrono::steady_clock::now();
  sink += acc.re._limbs[0];
  double dt4 = std::chrono::duration<double>(t2 - t1).count();

  auto t3 = std::chrono::steady_clock::now();
  acc = a0;
  for (int i = 0; i < M; ++i) acc = mul_k(acc, b0);
  auto t4 = std::chrono::steady_clock::now();
  sink += acc.re._limbs[0];
  double dtk = std::chrono::duration<double>(t4 - t3).count();

  printf("\nperf (%d iters, loop-carried):\n", M);
  printf("  4-mul   %.3f s  (%.2f ns/op)\n", dt4, dt4 * 1e9 / M);
  printf("  karats  %.3f s  (%.2f ns/op)\n", dtk, dtk * 1e9 / M);
  printf("  ratio karats/4-mul = %.3f\n", dtk / dt4);
  (void)sink;
  return 0;
}
