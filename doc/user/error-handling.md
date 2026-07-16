# Error handling

`multifloats` follows a strict **NaN-in-NaN-out** policy. Invalid inputs
propagate through the computation as IEEE 754 non-finite values rather than
interrupting control flow:

- **No `errno`.** `<cerrno>` is never read or written. `log(-1)` returns a
  DD NaN; `errno` is left untouched.
- **No `fenv` side effects.** The kernels do not call `feraiseexcept`,
  `feclearexcept`, or depend on the floating-point rounding mode.
- **No exceptions.** No C++ `throw`, no `std::terminate`, no `__builtin_trap`,
  no Fortran `error stop`.
- **No signaling-NaN handling.** Any NaN (quiet or signalling) is treated as
  data and propagated by the arithmetic kernels.
- **No input validation.** Domain checks are limited to what IEEE 754 naturally
  produces — `sqrt(-x)` returns NaN because the underlying `std::sqrt(hi)` does;
  there is no pre-check that raises an error.

This mirrors the behavior of `double` in C and `real(kind=16)` in Fortran, and
makes the kernels safe to call from hot loops, vectorized code, and parallel
regions without synchronizing on shared error state.

```{tip}
To detect non-finite results, use the classification functions on the leading
limb — `isnan`, `isinf`, `isfinite`, `signbit` (C++ ADL on `float64x2`) — after
the computation, the same way you would for `double`.
```
