# Missing C23 `<math.h>` Functions for `float64x2`

Audit of `include/multifloats/float64x2.h` against the C23 `<math.h>`
surface. These are candidates to implement later. Names below use the
existing `dd` suffix convention.

## C23-new functions (missing)

| Proposed name      | C23 name        | Notes                                                  |
| ------------------ | --------------- | ------------------------------------------------------ |
| `powndd`           | `pown`          | Integer exponent. We already ship `powidd`; adding the C23-spelled alias is cheap. |
| `powrdd`           | `powr`          | `exp(y*log(x))`, requires `x > 0`.                     |
| `rootndd`          | `rootn`         | `x^(1/n)`.                                             |
| `rsqrtdd`          | `rsqrt`         | `1/sqrt(x)`. Cheap, frequently useful.                 |
| `compoundndd`      | `compoundn`     | `(1+x)^n`.                                             |
| `roundevendd`      | `roundeven`     | Round to nearest, ties to even. Distinct from `rint`/`nearbyint` (those follow the current rounding mode). |
| `fmaximumdd`       | `fmaximum`      | C23 max: NaN propagates (unlike `fmax`).               |
| `fminimumdd`       | `fminimum`      | C23 min: NaN propagates.                               |
| `fmaximum_numdd`   | `fmaximum_num`  | NaN-suppressing variant.                               |
| `fminimum_numdd`   | `fminimum_num`  | NaN-suppressing variant.                               |
| `fmaximum_magdd`   | `fmaximum_mag`  | Max by magnitude.                                      |
| `fminimum_magdd`   | `fminimum_mag`  | Min by magnitude.                                      |
| `fmaximum_mag_numdd` | `fmaximum_mag_num` |                                                  |
| `fminimum_mag_numdd` | `fminimum_mag_num` |                                                  |
| `fmaxmagdd`        | `fmaxmag`       | Pre-C23 max-by-magnitude (kept for completeness).      |
| `fminmagdd`        | `fminmag`       | Pre-C23 min-by-magnitude.                              |
| `nextupdd`         | `nextup`        | Next representable toward `+inf`.                      |
| `nextdowndd`       | `nextdown`      | Next representable toward `-inf`.                      |
| `fromfpdd`         | `fromfp`        | Round to integer in integer type.                      |
| `ufromfpdd`        | `ufromfp`       |                                                        |
| `fromfpxdd`        | `fromfpx`       | `x` variants raise `inexact`.                          |
| `ufromfpxdd`       | `ufromfpx`      |                                                        |
| `canonicalizedd`   | `canonicalize`  | Trivial for IEEE-754 binary formats; semantics for DD need defining. |
| `totalorderdd`     | `totalorder`    | Total ordering predicate.                              |
| `totalordermagdd`  | `totalordermag` |                                                        |
| `getpayloaddd`     | `getpayload`    | NaN payload accessors.                                 |
| `setpayloaddd`     | `setpayload`    |                                                        |
| `setpayloadsigdd`  | `setpayloadsig` |                                                        |

## C99 functions (missing)

| Proposed name    | C99 name      | Notes                                                |
| ---------------- | ------------- | ---------------------------------------------------- |
| `remquodd`       | `remquo`      | `remainder` plus low quotient bits.                  |
| `scalblndd`      | `scalbln`     | `long` exponent variant of `scalbn`.                 |
| `nandd`          | `nan`         | Parse NaN payload string.                            |
| `fpclassifydd`   | `fpclassify`  | Returns `FP_NAN`/`FP_INFINITE`/`FP_NORMAL`/`FP_SUBNORMAL`/`FP_ZERO`. We have the individual predicates but no unified classifier. |
| `isnormaldd`     | `isnormal`    |                                                      |

## C23 classification macros (missing)

- `iscanonicaldd`
- `issignalingdd`
- `issubnormaldd`
- `iszerodd`

## Not applicable

- `nexttoward` — requires `long double`, no meaningful DD analogue.
- `isgreater` / `isgreaterequal` / `isless` / `islessequal` / `islessgreater` / `isunordered` — these are quiet-comparison wrappers; users can call the operators directly, no DD version needed.

## Suggested priority

High value, low cost:

1. `rsqrtdd` — useful, fast.
2. `roundevendd` — fills a real gap (mode-independent round-to-nearest-even).
3. `remquodd` — common in argument reduction.
4. `fpclassifydd`, `isnormaldd` — finishes the classification set.
5. `fmaximumdd` / `fminimumdd` family — C23 NaN-propagating semantics differ from existing `fmin`/`fmax`.

Medium value:

6. `nextupdd`, `nextdowndd`, `scalblndd`.
7. `powrdd`, `rootndd`, `compoundndd`, `powndd` (alias).

Low priority (rarely used; non-trivial semantics for DD):

8. `fromfp*`, `ufromfp*`, `canonicalize`, `totalorder*`, `getpayload`/`setpayload*`, `nandd`.
