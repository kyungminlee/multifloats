// Sweep matmul behavior across (k, renorm_interval) pairs to characterize
// the precision vs. perf tradeoff that DD_FMA_RENORM_INTERVAL trades off.
// Reference is a libquadmath dot product (__float128).
#include "../src/multifloats.hh"
#include "../src/multifloats_c.h"
#include <quadmath.h>
#include <cstdio>
#include <chrono>
#include <random>
#include <vector>

namespace mf = multifloats;
using q = __float128;

static q to_q(mf::float64x2 const& x) { return (q)x._limbs[0] + (q)x._limbs[1]; }

int main() {
  int const M = 4;
  int const N = 4;
  int const trials = 40;
  int const ks[] = {16, 64, 256, 1024, 4096, 16384, 65536};
  int const ris[] = {0, 4, 8, 16, 32, 64};

  printf("  k     ri  max_rel    time_s\n");
  printf("--------------------------------\n");

  for (int k : ks) {
    // Build A (M×k) and B (k×N) once.
    std::mt19937_64 rng(12345);
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    std::vector<float64x2_t> A(M * k), B(k * N);
    std::vector<q> Aq(M * k), Bq(k * N);
    for (int i = 0; i < M * k; ++i) {
      double hi = u(rng);
      A[i].hi = hi; A[i].lo = 0.0;
      Aq[i] = (q)hi;
    }
    for (int i = 0; i < k * N; ++i) {
      double hi = u(rng);
      B[i].hi = hi; B[i].lo = 0.0;
      Bq[i] = (q)hi;
    }
    // Quad-precision reference C (column-major).
    std::vector<q> Cref(M * N, 0);
    for (int j = 0; j < N; ++j)
      for (int p = 0; p < k; ++p)
        for (int i = 0; i < M; ++i)
          Cref[i + j * M] += Aq[i + p * M] * Bq[p + j * k];

    for (int ri : ris) {
      std::vector<float64x2_t> C(M * N);
      auto t1 = std::chrono::steady_clock::now();
      for (int rep = 0; rep < trials; ++rep)
        matmuldd_mm(A.data(), B.data(), C.data(), M, k, N, ri);
      auto t2 = std::chrono::steady_clock::now();
      double dt = std::chrono::duration<double>(t2 - t1).count();

      double max_rel = 0;
      for (int i = 0; i < M * N; ++i) {
        q got = (q)C[i].hi + (q)C[i].lo;
        q ref = Cref[i];
        q mag = fabsq(ref);
        if (mag < 1e-30q) continue;
        double rel = (double)(fabsq(got - ref) / mag);
        if (rel > max_rel) max_rel = rel;
      }
      printf("%6d %5d  %.3e  %.5f\n", k, ri, max_rel, dt / trials);
    }
    printf("--------------------------------\n");
  }
  return 0;
}
