
// Inline matmul micro-ops used by the dispatched panel templates below.
// `always_inline` is a hint: the panel templates rely on both ends of
// the MAC fusing into one register-resident body; a stray function call
// here would spill the accumulator ladder and cost ~2× the hot loop.
__attribute__((always_inline))
static inline void mac_inl(double ah, double al, double bh, double bl,
                              double &s_hi, double &s_lo) {
  double p = ah * bh;
  double e = std::fma(ah, bh, -p);
  double cross = std::fma(ah, bl, al * bh);
  double t = s_hi + p;
  double bp = t - s_hi;
  double ap = t - bp;
  double aerr = s_hi - ap;
  double berr = p - bp;
  s_lo += (e + cross) + (aerr + berr);
  s_hi = t;
}

// Renormalize an accumulator pair. Uses full two_sum rather than
// fast_two_sum because |s_hi| >= |s_lo| is not guaranteed after
// cancellation in a long compensated dot-product.
__attribute__((always_inline))
static inline void renorm_inl(double &s_hi, double &s_lo) {
  double t = s_hi + s_lo;
  double bp = t - s_hi;
  double ap = t - bp;
  double aerr = s_hi - ap;
  double berr = s_lo - bp;
  s_lo = aerr + berr;
  s_hi = t;
}

__attribute__((always_inline))
static inline float64x2 finalize_inl(double s_hi, double s_lo) {
  renorm_inl(s_hi, s_lo);
  return {s_hi, s_lo};
}

// AXPY-style panel µkernel. Processes a row-panel of A of compile-time
// height MR against a contiguous x / b-column, writing MR finalized DD
// outputs. `lda` is the leading dim of A (so the panel can sit in a
// larger matrix); `renorm_interval > 0` triggers a two_sum of each
// accumulator pair every `renorm_interval` reductions (keeps s_lo
// bounded and precision ~DD for large k, matching dot_product).
template <int MR>
static inline void gaxpy_mv_panel(const float64x2 *__restrict__ a,
                                     const float64x2 *__restrict__ x,
                                     float64x2 *__restrict__ y, int64_t lda,
                                     int64_t k, int64_t renorm_interval) {
  double s_hi[MR] = {}, s_lo[MR] = {};
  // Fast path: if no intermediate renorm is needed, drop the chunking
  // scaffolding entirely — the p-loop is a single tight sweep.
  if (renorm_interval <= 0 || k <= renorm_interval) {
    for (int64_t p = 0; p < k; ++p) {
      const double xh = x[p].limbs[0];
      const double xl = x[p].limbs[1];
      const float64x2 *__restrict__ acol = a + p * lda;
      for (int i = 0; i < MR; ++i) {
        mac_inl(acol[i].limbs[0], acol[i].limbs[1], xh, xl, s_hi[i], s_lo[i]);
      }
    }
    for (int i = 0; i < MR; ++i) y[i] = finalize_inl(s_hi[i], s_lo[i]);
    return;
  }
  // Chunked path: for large k, do `renorm_interval` reductions then
  // two_sum each accumulator pair to keep s_lo bounded.
  const int64_t chunk = renorm_interval;
  int64_t p0 = 0;
  while (p0 < k) {
    int64_t pend = p0 + chunk;
    if (pend > k) pend = k;
    for (int64_t p = p0; p < pend; ++p) {
      const double xh = x[p].limbs[0];
      const double xl = x[p].limbs[1];
      const float64x2 *__restrict__ acol = a + p * lda;
      for (int i = 0; i < MR; ++i) {
        mac_inl(acol[i].limbs[0], acol[i].limbs[1], xh, xl, s_hi[i], s_lo[i]);
      }
    }
    p0 = pend;
    if (p0 < k) {
      for (int i = 0; i < MR; ++i) renorm_inl(s_hi[i], s_lo[i]);
    }
  }
  for (int i = 0; i < MR; ++i) y[i] = finalize_inl(s_hi[i], s_lo[i]);
}

// Cold tail handler: 1..7 leftover rows. Kept out-of-line so the hot
// dispatcher stays small enough for the compiler to inline panel<8>
// into matmuldd_mv — otherwise the 7 extra template instantiations
// would bloat the caller and push the accumulators out of registers.
__attribute__((noinline))
static void gaxpy_mv_tail(const float64x2 *__restrict__ a,
                             const float64x2 *__restrict__ x,
                             float64x2 *__restrict__ y, int tail,
                             int64_t lda, int64_t k,
                             int64_t renorm_interval) {
  switch (tail) {
    case 1: gaxpy_mv_panel<1>(a, x, y, lda, k, renorm_interval); break;
    case 2: gaxpy_mv_panel<2>(a, x, y, lda, k, renorm_interval); break;
    case 3: gaxpy_mv_panel<3>(a, x, y, lda, k, renorm_interval); break;
    case 4: gaxpy_mv_panel<4>(a, x, y, lda, k, renorm_interval); break;
    case 5: gaxpy_mv_panel<5>(a, x, y, lda, k, renorm_interval); break;
    case 6: gaxpy_mv_panel<6>(a, x, y, lda, k, renorm_interval); break;
    case 7: gaxpy_mv_panel<7>(a, x, y, lda, k, renorm_interval); break;
  }
}

// Tile any m into MR=8 row panels plus a 1..7 tail. Keep this small so
// the hot panel<8> body inlines — accumulators need to stay in registers.
static inline void gaxpy_mv_dispatch(const float64x2 *__restrict__ a,
                                        const float64x2 *__restrict__ x,
                                        float64x2 *__restrict__ y, int64_t m,
                                        int64_t k, int64_t lda,
                                        int64_t renorm_interval) {
  constexpr int MR = 8;
  int64_t i = 0;
  for (; i + MR <= m; i += MR) {
    gaxpy_mv_panel<MR>(a + i, x, y + i, lda, k, renorm_interval);
  }
  int tail = static_cast<int>(m - i);
  if (tail > 0) gaxpy_mv_tail(a + i, x, y + i, tail, lda, k, renorm_interval);
}

// MR×NR GEMM µkernel: loads one column-slice of A per p and reuses it
// across NR B-column entries. Amortizes A bandwidth vs a per-column mv
// dispatch at the cost of MR*NR accumulator pairs (stack-spilled on
// machines with < 32 FP regs, but the spill cost is paid once per p
// and dwarfed by the DD-mac FLOPs).
template <int MR, int NR>
static inline void gemm_panel(const float64x2 *__restrict__ a,
                                 const float64x2 *__restrict__ b,
                                 float64x2 *__restrict__ c,
                                 int64_t lda, int64_t ldb, int64_t ldc,
                                 int64_t k, int64_t renorm_interval) {
  double s_hi[MR][NR] = {}, s_lo[MR][NR] = {};
  auto kernel_block = [&](int64_t p0, int64_t pend) {
    for (int64_t p = p0; p < pend; ++p) {
      const float64x2 *__restrict__ acol = a + p * lda;
      double ah[MR], al[MR];
      for (int i = 0; i < MR; ++i) { ah[i] = acol[i].limbs[0]; al[i] = acol[i].limbs[1]; }
      for (int jj = 0; jj < NR; ++jj) {
        const double bh = b[jj * ldb + p].limbs[0];
        const double bl = b[jj * ldb + p].limbs[1];
        for (int i = 0; i < MR; ++i) {
          mac_inl(ah[i], al[i], bh, bl, s_hi[i][jj], s_lo[i][jj]);
        }
      }
    }
  };
  if (renorm_interval <= 0 || k <= renorm_interval) {
    kernel_block(0, k);
  } else {
    const int64_t chunk = renorm_interval;
    int64_t p0 = 0;
    while (p0 < k) {
      int64_t pend = p0 + chunk;
      if (pend > k) pend = k;
      kernel_block(p0, pend);
      p0 = pend;
      if (p0 < k) {
        for (int jj = 0; jj < NR; ++jj)
          for (int i = 0; i < MR; ++i)
            renorm_inl(s_hi[i][jj], s_lo[i][jj]);
      }
    }
  }
  for (int jj = 0; jj < NR; ++jj)
    for (int i = 0; i < MR; ++i)
      c[jj * ldc + i] = finalize_inl(s_hi[i][jj], s_lo[i][jj]);
}
