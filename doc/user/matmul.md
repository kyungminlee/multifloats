# Matmul

The `matmul` surface provides three shape-dispatched operations, for both
real (`real64x2`) and complex (`cmplx64x2`) operands, through the same Fortran
generic:

| Op             | Signature                       |
| -------------- | ------------------------------- |
| `matmul(A, B)` | `A(m,k) * B(k,n) → C(m,n)`      |
| `matmul(A, x)` | `A(m,k) * x(k)   → y(m)`        |
| `matmul(x, B)` | `x(k)   * B(k,n) → y(n)`        |

The C ABI exposes the three kernels directly as `matmuldd_mm` / `matmuldd_mv`
/ `matmuldd_vm` (complex matmul is composed in Fortran from four real matmuls).

## Scope vs GEMM

`multifloats::matmul` is **not** a drop-in replacement for a GEMM:

- **No `transa` / `transb` flags.** Inputs are always interpreted in storage
  order. Materialize the transpose first (`matmul(transpose(A), B)`).
- **No `alpha` / `beta` scaling.** The output is overwritten, not accumulated.
  Build `C := α·A·B + β·C` externally.
- **No leading dimension (LDA/LDB/LDC).** Storage is contiguous column-major
  with leading dimension equal to the first extent.

These constraints are deliberate: the compensated DD kernels use a
register-blocked panel design that assumes contiguous column-major storage.
Relaxing them (GEMM-style trans/alpha/beta/LDA) is tracked as deferred work in
[`doc/dev/architecture.md`](https://github.com/kyungminlee/multifloats/blob/main/doc/dev/architecture.md).

## Renormalization interval

Matmul and `dot_product` run a compensated fused-multiply-accumulate loop. The
low-limb accumulator grows by ~1 ULP per iteration, so for large inner
dimension `k` it must be re-normalized periodically. The module constant
`DD_FMA_RENORM_INTERVAL` (default `8`) controls the frequency; set `0` to
disable periodic renorm.

| k     | ri=0    | ri=8    | ri=64   |
| ----- | ------- | ------- | ------- |
| 64    | 2.9e-31 | 2.9e-31 | 2.9e-31 |
| 1024  | 3.2e-29 | 4.9e-30 | 3.3e-30 |
| 65536 | 1.5e-28 | 3.2e-30 | 1.5e-29 |

`ri=8` is near-optimal at ~3–4% overhead. Call
`dd_set_fma_renorm_interval(n)` from Fortran (or pass `n` to the `matmuldd_*`
C kernels) to override per call. For `k < 100`, `ri=0` is fine; for
`k > 10000`, keep `ri=8`; avoid `ri > 32` at very large `k` since `s_lo` can
alias off significant bits before the next renormalization fires.
