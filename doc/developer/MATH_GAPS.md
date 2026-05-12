# Math Surface Gaps for `float64x2` / `real64x2`

Audit of the public math surface against:

- C23 `<math.h>` — `include/multifloats/float64x2.h`
- gfortran intrinsic procedures
  (https://gcc.gnu.org/onlinedocs/gfortran/Intrinsic-Procedures.html)
  and Intel Fortran intrinsic categories
  (https://www.intel.com/content/www/us/en/docs/fortran-compiler/developer-guide-reference/2024-2/categories-of-intrinsic-functions.html)
  — `fsrc/multifloats.fypp`

These are candidates to implement later.

---

# Part 1 — C23 `<math.h>` gaps

Names use the existing `dd` suffix convention.

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

---

# Part 2 — Fortran intrinsic gaps

Audited against gfortran's intrinsic procedure list and Intel Fortran's
intrinsic categories. Existing coverage lives in
`fsrc/multifloats.fypp` (search the `UNARY_DD_DD`, `UNARY_CDD_CDD`,
`UNARY_CDD_DD`, `BINARY_REAL`, `BINARY_MINMAX`, `DD_INT_TO_DD`,
`INQUIRY_*`, `ARRAY_REDUCE_*`, `BESSEL_N` sets).

## Trigonometric — degree variants (F2023)

F2023 standardized degree-argument trig. None are currently provided
for `real64x2` / `cmplx64x2`.

| Intrinsic   | Notes                                 |
| ----------- | ------------------------------------- |
| `sind`      | sin in degrees                        |
| `cosd`      | cos in degrees                        |
| `tand`      | tan in degrees                        |
| `asind`     | arcsin returning degrees              |
| `acosd`     | arccos returning degrees              |
| `atand`     | arctan returning degrees              |
| `atan2d`    | atan2 returning degrees               |

## Trigonometric — cotangent family (extension)

gfortran/Intel extension; non-standard but widely used.

| Intrinsic   | Notes                                 |
| ----------- | ------------------------------------- |
| `cotan`     | cot(x) (radians)                      |
| `cotand`    | cot(x) (degrees)                      |

(Intel additionally documents `acot`, `acotd`, `acosh`/`asinh`/`atanh`
which we already provide.)

## Numeric inquiry / manipulation

| Intrinsic       | Status                                                                     |
| --------------- | -------------------------------------------------------------------------- |
| `out_of_range`  | F2018; tests if a real value would overflow target integer/real kind. Useful for guarding `int(dd)` / `nint(dd)` casts. Missing. |
| `selected_real_kind` | N/A — kind selector, not a value function.                            |
| `kind`          | N/A on derived type.                                                       |

## IEEE arithmetic (`ieee_arithmetic` module)

We use `ieee_is_finite` / `ieee_is_nan` internally but expose nothing
for DD. A DD-aware `ieee_arithmetic`-style interface would cover:

- `ieee_is_finite`, `ieee_is_nan`, `ieee_is_normal`, `ieee_is_negative`
- `ieee_class` (returns `ieee_class_type`)
- `ieee_value` (build NaN / Inf / signed zero of given class)
- `ieee_copy_sign`, `ieee_signbit`
- `ieee_logb`, `ieee_scalb`, `ieee_rint`, `ieee_rem`, `ieee_next_after`
- `ieee_unordered`, `ieee_signaling_compare`
- `ieee_max_num`, `ieee_min_num`, `ieee_max_num_mag`, `ieee_min_num_mag`
  (mirror of the C23 `fmaximum`/`fminimum` family)
- `ieee_get_underflow_mode` / set — N/A (mode is hardware).
- `ieee_get_rounding_mode` / set — N/A.

These are largely thin wrappers; many overlap with the C23 list above
(e.g. `ieee_max_num` ↔ `fmaximum_num`).

## Bit / integer-only intrinsics

`ishft`, `ishftc`, `iand`, `ior`, `ieor`, `not`, `popcnt`, `poppar`,
`leadz`, `trailz`, `bge`, `bgt`, `ble`, `blt`, `btest`, `ibclr`,
`ibset`, `ibits`, `mvbits`, `dshiftl`, `dshiftr`, `maskl`, `maskr` —
all integer-only, **N/A** for DD.

## Array intrinsics

These are generic over derived types and should already work without
explicit overloads:

`all`, `any`, `count`, `parity`, `merge`, `pack`, `unpack`, `cshift`,
`eoshift`, `reshape`, `spread`, `transpose`, `lbound`, `ubound`,
`shape`, `size`, `rank`.

If we observe miscompiles or surprising behavior with derived types in
any of these, they become candidates for explicit overloads — but
nothing is known-broken today.

## String / I/O / system

`achar`, `iachar`, `char`, `ichar`, `len`, `len_trim`, `trim`,
`adjustl`, `adjustr`, `index`, `scan`, `verify`, `repeat`,
`new_line`, `command_argument_count`, `get_command*`,
`get_environment_variable`, `system_clock`, `cpu_time`, `date_and_time`
— **N/A**, not numeric.

## Suggested priority (Fortran)

High value:

1. `sind`, `cosd`, `tand`, `asind`, `acosd`, `atand`, `atan2d` —
   F2023 standard, easy via the existing `*pi` infrastructure scaled by
   `pi/180`.
2. `out_of_range` for `real64x2` — guards integer conversions.
3. A small `ieee_arithmetic`-style interface: `ieee_is_finite`,
   `ieee_is_nan`, `ieee_is_normal`, `ieee_class`, `ieee_value`,
   `ieee_copy_sign`. These are cheap and cover the common
   classification/construction needs.

Medium value:

4. `cotan`, `cotand` — gfortran/Intel extension; trivial.
5. Remaining `ieee_*` (`logb`, `scalb`, `rint`, `rem`, `next_after`,
   `unordered`, `max_num`/`min_num` family).

Low priority:

6. `ieee_signaling_compare` and other rarely-used IEEE entries.
