"""Op metadata driving BENCHMARK.md rendering.

Each op describes:
  - ``bench_key``   : the label printed by ``fortran_bench`` / ``cpp_bench``
                     (used to look up the speedup in the per-system JSON).
  - ``display``     : the markdown-escaped label shown in the BENCHMARK.md row.
  - ``approach``    : prose description of the algorithm (one cell).
  - ``prec``        : short precision label (e.g. ``"full DD"``, ``"exact"``).
                     Only the Fortran tables include a ``prec`` column.
  - ``fuzz``        : how to assemble the err cell from the fuzz stats.
                       * ``None``     — no fuzz data (shows "—")
                       * ``"exact"``  — always shows "exact"
                       * ``"<key>"``  — pulls one ``max_rel`` value
                       * list of ``("re"|"im", "<key>")`` — shows as
                         ``"<v_re> (re) / <v_im> (im)"``

Sections are the H3 headings inside each H2 ("Fortran" / "C++").
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence, Tuple, Union


FuzzSpec = Union[None, str, Sequence[Tuple[str, str]]]


@dataclass
class Op:
    bench_key: str
    display: str
    approach: str
    prec: str = ""
    fuzz: FuzzSpec = None


@dataclass
class Section:
    title: str
    ops: list[Op] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Fortran: float64x2 vs real(kind=16)
# ---------------------------------------------------------------------------

FORTRAN_SECTIONS: list[Section] = [
    Section("Arithmetic", [
        Op("add",         "add",         "Julia: two\\_sum EFT",                                   "full DD", "add"),
        Op("sub",         "sub",         "Julia: two\\_sum EFT (negate + add)",                    "full DD", "sub"),
        Op("mul",         "mul",         "Julia: two\\_prod EFT via FMA",                          "full DD", "mul"),
        Op("div",         "div",         "original: Newton refinement (1/y seed, one step)",       "full DD", "div"),
        Op("sqrt",        "sqrt",        "Julia: Karp\u2013Markstein (reciprocal sqrt seed + Newton)", "full DD", "sqrt"),
        Op("add (dd+dp)", "add (dd+dp)", "Julia: two\\_sum EFT",                                   "exact",   "exact"),
        Op("mul (dp*dd)", "mul (dp\\*dd)", "Julia: two\\_prod EFT via FMA",                        "full DD", "mul_df"),
    ]),
    Section("Unary", [
        Op("abs",          "abs",          "original: sign-check + negate limbs",               "exact", "abs"),
        Op("neg",          "neg",          "original: negate both limbs",                       "exact", "neg"),
        Op("aint",         "aint",         "original: truncate hi, check DD fractional part",   "exact", "trunc"),
        Op("anint",        "anint",        "original: truncate hi, DD fractional part vs \u00b10.5", "exact", "round"),
        Op("fraction",     "fraction",     "original: scale both limbs by \u2212exponent",      "exact", "exact"),
        Op("scale",        "scale",        "original: ldexp on both limbs",                     "exact", "exact"),
        Op("set_exponent", "set\\_exponent", "original: scale + set\\_exponent on hi",           "exact", "exact"),
    ]),
    Section("Binary", [
        Op("min",    "min",    "original: DD comparison + select",                         "full DD", "fmin"),
        Op("max",    "max",    "original: DD comparison + select",                         "full DD", "fmax"),
        Op("min3",   "min3",   "original: chained min",                                    "full DD", "min3"),
        Op("max3",   "max3",   "original: chained max",                                    "full DD", "max3"),
        Op("sign",   "sign",   "original: sign-check + negate",                            "exact",   "copysign"),
        Op("dim",    "dim",    "original: DD comparison, then subtract or zero",           "full DD", "fdim"),
        Op("hypot",  "hypot",  "original: scaled sqrt(x\u00b2+y\u00b2)",                   "full DD", "hypot"),
        Op("mod",    "mod",    "sample: floor-multiple reduction loop; fallback to div chain", "full DD", "fmod"),
        Op("modulo", "modulo", "original: mod + sign adjustment",                          "full DD", "fmod"),
    ]),
    Section("Exponential / logarithmic", [
        Op("exp",     "exp",     "Julia: exp2 polynomial (14-term Horner) + ldexp reconstruction",  "full DD", "exp"),
        Op("log",     "log",     "Julia: log2 table lookup (32 centers) + polynomial (7-term Horner)", "full DD", "log"),
        Op("log10",   "log10",   "Julia: log2 kernel \u00d7 DD log10(2)",                           "full DD", "log10"),
        Op("pow",     "pow",     "Julia: exp(y \u00d7 log(x))",                                     "full DD", "pow"),
        Op("pow_int", "pow\\_int", "original: repeated squaring via DD mul",                        "full DD", "pow_int"),
    ]),
    Section("Trigonometric", [
        Op("sin",    "sin",    "original: 13-term Taylor Horner + 3-part Cody\u2013Waite \u03c0/2 + \u03c0/8 split", "full DD", "sin"),
        Op("cos",    "cos",    "original: 13-term Taylor Horner + 3-part Cody\u2013Waite \u03c0/2 + \u03c0/8 split", "full DD", "cos"),
        Op("sincos", "sincos", "fused sin/cos: one range-reduction, two Taylor outputs",                             "full DD",
            fuzz=[("sin", "sincos_s"), ("cos", "sincos_c")]),
        Op("sinpi",  "sinpi",  "Julia: sinpi Horner polynomial, direct",                                              "full DD", "sinpi"),
        Op("cospi",  "cospi",  "Julia: cospi Horner polynomial, direct",                                              "full DD", "cospi"),
        Op("tan",    "tan",    "original: sin/cos Taylor kernels + DD divide",                                        "full DD", "tan"),
        Op("tanpi",  "tanpi",  "original: sinpi/cospi ratio",                                                         "full DD", "tanpi"),
        Op("asin",   "asin",   "original: piecewise rational P/Q (3 regions, from libquadmath asinq.c)",              "full DD", "asin"),
        Op("asinpi", "asinpi", "original: asin(x)/\u03c0 with exact-DD \u03c0 division",                               "full DD", "asinpi"),
        Op("acos",   "acos",   "original: asin polynomial + half-angle identity",                                     "full DD", "acos"),
        Op("acospi", "acospi", "original: acos(x)/\u03c0 with exact-DD \u03c0 division",                               "full DD", "acospi"),
        Op("atan",   "atan",   "original: 84-entry table lookup + rational P(t\u00b2)/Q(t\u00b2) (from libquadmath atanq.c)", "full DD", "atan"),
        Op("atanpi", "atanpi", "original: atan(x)/\u03c0 with exact-DD \u03c0 division",                               "full DD", "atanpi"),
        Op("atan2",  "atan2",  "original: table-based atan + quadrant correction",                                    "full DD", "atan2"),
        Op("atan2pi","atan2pi","original: atan2(y,x)/\u03c0 with exact-DD \u03c0 division",                            "full DD", "atan2pi"),
    ]),
    Section("Hyperbolic", [
        Op("sinh",  "sinh",  "original: Taylor series (\\|x\\|<0.1) or (exp\u2212exp\u207b\u00b9)/2", "full DD", "sinh"),
        Op("cosh",  "cosh",  "original: (exp+exp\u207b\u00b9)/2",                                    "full DD", "cosh"),
        Op("sinhcosh", "sinhcosh", "fused sinh/cosh: one range-reduction, two outputs",              "full DD",
            fuzz=[("sinh", "sinhcosh_s"), ("cosh", "sinhcosh_c")]),
        Op("tanh",  "tanh",  "original: sinh/cosh (\\|x\\|<0.5) or (1\u2212e\u207b\u00b2\u02e3)/(1+e\u207b\u00b2\u02e3)", "full DD", "tanh"),
        Op("asinh", "asinh", "original: Taylor series (\\|x\\|<0.01) or log(x+\u221a(x\u00b2+1)) with Newton", "full DD", "asinh"),
        Op("acosh", "acosh", "original: log(x+\u221a(x\u00b2\u22121)) with Newton correction",       "full DD", "acosh"),
        Op("atanh", "atanh", "original: Taylor series (\\|x\\|<0.01) or \u00bd\u00b7log((1+x)/(1\u2212x))", "full DD", "atanh"),
    ]),
    Section("Error / special functions", [
        Op("erf",           "erf",           "piecewise rational approx (libquadmath erfq.c)",                "full DD", "erf"),
        Op("erfc",          "erfc",          "piecewise rational approx + split exp(-x^2)",                   "full DD", "erfc"),
        Op("erfc_scaled",   "erfc\\_scaled", "exp(x^2)\u00b7erfc(x) with asymptotic cancellation",            "full DD", "erfcx"),
        Op("gamma",         "gamma",         "piecewise rational approx + Stirling + reflection",             "full DD", "gamma"),
        Op("log_gamma",     "log\\_gamma",   "piecewise rational approx + Stirling asymptotic",               "full DD", "lngamma"),
        Op("bessel_j0",     "bessel\\_j0",   "piecewise rational + Hankel asymptotic (j0q.c) via C++",        "full DD", "bj0"),
        Op("bessel_j1",     "bessel\\_j1",   "piecewise rational + Hankel asymptotic (j1q.c) via C++",        "full DD", "bj1"),
        Op("bessel_jn(3,.)","bessel\\_jn(3,.)", "forward/backward recurrence from j0/j1",                     "full DD", "bjn"),
        Op("bessel_y0",     "bessel\\_y0",   "piecewise rational + Hankel asymptotic (j0q.c) via C++",        "full DD", "by0"),
        Op("bessel_y1",     "bessel\\_y1",   "piecewise rational + Hankel asymptotic (j1q.c) via C++",        "full DD", "by1"),
        Op("bessel_yn(3,.)","bessel\\_yn(3,.)", "forward recurrence from y0/y1",                              "full DD", "byn"),
    ]),
    Section("Complex arithmetic", [
        # Bench key stays `cdd_*` (matches bench.cc output); fuzz keys use the
        # libquadmath-style `c*` naming emitted by fuzz.cc's CHK_C1 / CHK_COP.
        Op("cdd_add",   "cdd\\_add",   "original: component-wise DD add",      "full DD",
            fuzz=[("re", "cadd_re"), ("im", "cadd_im")]),
        Op("cdd_sub",   "cdd\\_sub",   "original: component-wise DD sub",      "full DD",
            fuzz=[("re", "csub_re"), ("im", "csub_im")]),
        Op("cdd_mul",   "cdd\\_mul",   "original: (ac\u2212bd, ad+bc) via DD ops", "full DD",
            fuzz=[("re", "cmul_re"), ("im", "cmul_im")]),
        Op("cdd_div",   "cdd\\_div",   "original: (ac+bd, bc\u2212ad)/(c\u00b2+d\u00b2)", "full DD / deriv",
            fuzz=[("re", "cdiv_re"), ("im", "cdiv_im")]),
        Op("cdd_conjg", "cdd\\_conjg", "original: negate im limbs",            "exact", "exact"),
        Op("cdd_abs",   "cdd\\_abs",   "original: hypot(re, im)",              "full DD", "cabs"),
    ]),
    Section("Complex transcendentals", [
        Op("cdd_sqrt",  "cdd\\_sqrt",  "original: Kahan-style (\\|z\\|+\\|a\\|)/2 with scaling", "full DD",
            fuzz=[("re", "csqrt_re"), ("im", "csqrt_im")]),
        Op("cdd_exp",   "cdd\\_exp",   "original: exp(re)\u00b7(cos(im), sin(im))", "full DD",
            fuzz=[("re", "cexp_re"), ("im", "cexp_im")]),
        Op("cdd_expm1", "cdd\\_expm1", "original: expm1 + complex rotation with cancellation-safe Re", "reduced DD",
            fuzz=[("re", "cexpm1_re"), ("im", "cexpm1_im")]),
        Op("cdd_log",   "cdd\\_log",   "original: (log(\\|z\\|), atan2(im,re))", "full DD",
            fuzz=[("re", "clog_re"), ("im", "clog_im")]),
        Op("cdd_log2",  "cdd\\_log2",  "original: clog / log(2) component-wise",            "full DD",
            fuzz=[("re", "clog2_re"), ("im", "clog2_im")]),
        Op("cdd_log10", "cdd\\_log10", "original: clog / log(10) component-wise",           "full DD",
            fuzz=[("re", "clog10_re"), ("im", "clog10_im")]),
        Op("cdd_log1p", "cdd\\_log1p", "original: cancellation-safe log(1+z) near z=0",     "full DD",
            fuzz=[("re", "clog1p_re"), ("im", "clog1p_im")]),
        Op("cdd_pow",   "cdd\\_pow",   "original: exp(w\u00b7log(z))",                      "full DD",
            fuzz=[("re", "cpow_re"), ("im", "cpow_im")]),
        Op("cdd_sin",   "cdd\\_sin",   "original: sin(re)cosh(im), cos(re)sinh(im)", "full DD",
            fuzz=[("re", "csin_re"), ("im", "csin_im")]),
        Op("cdd_cos",   "cdd\\_cos",   "original: cos(re)cosh(im), \u2212sin(re)sinh(im)", "full DD",
            fuzz=[("re", "ccos_re"), ("im", "ccos_im")]),
        Op("cdd_tan",   "cdd\\_tan",   "original: complex sin/cos ratio", "full DD",
            fuzz=[("re", "ctan_re"), ("im", "ctan_im")]),
        Op("cdd_sinh",  "cdd\\_sinh",  "original: sinh(re)cos(im), cosh(re)sin(im)", "full DD",
            fuzz=[("re", "csinh_re"), ("im", "csinh_im")]),
        Op("cdd_cosh",  "cdd\\_cosh",  "original: cosh(re)cos(im), sinh(re)sin(im)", "full DD",
            fuzz=[("re", "ccosh_re"), ("im", "ccosh_im")]),
        Op("cdd_tanh",  "cdd\\_tanh",  "original: complex tanh via sinh/cosh", "full DD",
            fuzz=[("re", "ctanh_re"), ("im", "ctanh_im")]),
        Op("cdd_asin",  "cdd\\_asin",  "original: \u2212i\u00b7log(iz+\u221a(1\u2212z\u00b2))", "deriv / full DD",
            fuzz=[("re", "casin_re"), ("im", "casin_im")]),
        Op("cdd_acos",  "cdd\\_acos",  "original: \u03c0/2 \u2212 asin(z)", "full DD",
            fuzz=[("re", "cacos_re"), ("im", "cacos_im")]),
        Op("cdd_atan",  "cdd\\_atan",  "original: (i/2)\u00b7log((i+z)/(i\u2212z))", "full DD",
            fuzz=[("re", "catan_re"), ("im", "catan_im")]),
        Op("cdd_asinh", "cdd\\_asinh", "original: log(z+\u221a(z\u00b2+1))", "deriv / full DD",
            fuzz=[("re", "casinh_re"), ("im", "casinh_im")]),
        Op("cdd_acosh", "cdd\\_acosh", "original: log(z+\u221a(z\u00b2\u22121))", "full DD",
            fuzz=[("re", "cacosh_re"), ("im", "cacosh_im")]),
        Op("cdd_atanh", "cdd\\_atanh", "original: \u00bd\u00b7log((1+z)/(1\u2212z))", "deriv / full DD",
            fuzz=[("re", "catanh_re"), ("im", "catanh_im")]),
    ]),
    Section("Array reductions", [
        Op("arr_sum (n=8)",       "arr\\_sum (n=8)",       "original: chained DD add",                                                                       "full DD", "arr_sum"),
        Op("arr_product (n=8)",   "arr\\_product (n=8)",   "original: chained DD mul",                                                                       "full DD", "arr_prod"),
        Op("arr_maxval (n=8)",    "arr\\_maxval (n=8)",    "original: chained DD compare",                                                                   "full DD", "arr_max"),
        Op("arr_minval (n=8)",    "arr\\_minval (n=8)",    "original: chained DD compare",                                                                   "full DD", "arr_min"),
        Op("arr_dot (n=8)",       "arr\\_dot (n=8)",       "original: fused multiply-accumulate with periodic renormalization",                              "full DD", "arr_dot"),
        Op("arr_norm2 (n=8)",     "arr\\_norm2 (n=8)",     "original: sqrt(dot(x,x))",                                                                       "full DD", "arr_norm2"),
        Op("arr_matmul (8x8*8)",  "arr\\_matmul (8\u00d78\u00b78)", "original: AXPY-order C kernel, MR=8 register-blocked panels + 1..7 tail, periodic renorm", "full DD", "arr_matmul"),
    ]),
]


# ---------------------------------------------------------------------------
# C: C ABI (sindd / cdd_muldd / j0dd / ...) vs __float128
#
# Scope: every function exported from include/multifloats.h. Shared ops
# (sin, cdd_mul, …) are duplicated between C_SECTIONS and FORTRAN_SECTIONS
# so the Fortran-elemental ABI overhead is visible per-op. C-section-only
# ops: none — everything here is also callable from Fortran via the
# iso_c_binding bridge. Fortran-section-only ops (not in C ABI): array
# reductions, numeric inquiry — see FORTRAN_SECTIONS.
# ---------------------------------------------------------------------------

C_SECTIONS: list[Section] = [
    Section("Arithmetic", [
        Op("add",   "add",   "Julia: two\\_sum EFT",                                        fuzz="add"),
        Op("sub",   "sub",   "Julia: two\\_sum EFT (negate + add)",                         fuzz="sub"),
        Op("mul",   "mul",   "Julia: two\\_prod EFT via FMA",                               fuzz="mul"),
        Op("div",   "div",   "original: Newton refinement (1/y seed, one step)",            fuzz="div"),
        Op("sqrt",  "sqrt",  "Julia: Karp\u2013Markstein (reciprocal sqrt seed + Newton)",  fuzz="sqrt"),
        Op("fma",   "fma",   "original: x\\*y + z via DD ops",                              fuzz="fmadd"),
        Op("abs",   "abs",   "original: sign-check + negate limbs",                         fuzz="abs"),
        Op("neg",   "neg",   "original: negate both limbs",                                 fuzz="neg"),
    ]),
    Section("Rounding", [
        Op("trunc",     "trunc",     "original: signbit ? \u2212floor(\u2212x) : floor(x)", fuzz="exact"),
        Op("round",     "round",     "original: trunc(x + \u00bd\u00b7sign(x))",            fuzz="exact"),
    ]),
    Section("Binary", [
        Op("fmin",     "fmin",     "original: DD comparison + select",                           fuzz="fmin"),
        Op("fmax",     "fmax",     "original: DD comparison + select",                           fuzz="fmax"),
        Op("fdim",     "fdim",     "original: DD comparison, then subtract or zero",             fuzz="fdim"),
        Op("copysign", "copysign", "original: sign-bit copy to hi, propagate to lo",             fuzz="copysign"),
        Op("fmod",     "fmod",     "sample: floor-multiple reduction loop; fallback to div chain", fuzz="fmod"),
        Op("hypot",    "hypot",    "original: scaled sqrt(x\u00b2+y\u00b2)",                     fuzz="hypot"),
    ]),
    Section("Exponential / logarithmic", [
        Op("exp",    "exp",    "Julia: exp2 polynomial (14-term Horner) + ldexp reconstruction",   fuzz="exp"),
        Op("exp2",   "exp2",   "Julia: exp2 polynomial (14-term Horner)",                          fuzz="exp2"),
        Op("expm1",  "expm1",  "original: exp(x) \u2212 1 via DD sub",                             fuzz="expm1"),
        Op("log",    "log",    "Julia: log2 table lookup (32 centers) + polynomial (7-term Horner)", fuzz="log"),
        Op("log10",  "log10",  "Julia: log2 kernel \u00d7 DD log10(2)",                            fuzz="log10"),
        Op("log2",   "log2",   "Julia: log2 table lookup + polynomial",                            fuzz="log2"),
        Op("log1p",  "log1p",  "original: log(1 + x) via DD add",                                  fuzz="log1p"),
        Op("pow",    "pow",    "Julia: exp(y \u00d7 log(x))",                                      fuzz="pow"),
    ]),
    Section("Trigonometric", [
        Op("sin",    "sin",    "original: 13-term Taylor Horner + 3-part Cody\u2013Waite \u03c0/2 + \u03c0/8 split", fuzz="sin"),
        Op("cos",    "cos",    "original: 13-term Taylor Horner + 3-part Cody\u2013Waite \u03c0/2 + \u03c0/8 split", fuzz="cos"),
        Op("sincos", "sincos", "fused sin/cos: one range-reduction, two Taylor outputs",
            fuzz=[("sin", "sincos_s"), ("cos", "sincos_c")]),
        Op("sinpi",  "sinpi",  "Julia: sinpi Horner polynomial, direct",                            fuzz="sinpi"),
        Op("cospi",  "cospi",  "Julia: cospi Horner polynomial, direct",                            fuzz="cospi"),
        Op("tan",    "tan",    "original: sin/cos Taylor kernels + DD divide",                      fuzz="tan"),
        Op("tanpi",  "tanpi",  "original: sinpi/cospi ratio",                                       fuzz="tanpi"),
        Op("asin",   "asin",   "original: piecewise rational P/Q (3 regions, from libquadmath asinq.c)", fuzz="asin"),
        Op("asinpi", "asinpi", "original: asin(x)/\u03c0 with exact-DD \u03c0 division",            fuzz="asinpi"),
        Op("acos",   "acos",   "original: asin polynomial + half-angle identity",                   fuzz="acos"),
        Op("acospi", "acospi", "original: acos(x)/\u03c0 with exact-DD \u03c0 division",            fuzz="acospi"),
        Op("atan",   "atan",   "original: 84-entry table lookup + rational P(t\u00b2)/Q(t\u00b2) (from libquadmath atanq.c)", fuzz="atan"),
        Op("atanpi", "atanpi", "original: atan(x)/\u03c0 with exact-DD \u03c0 division",            fuzz="atanpi"),
        Op("atan2",  "atan2",  "original: table-based atan + quadrant correction",                  fuzz="atan2"),
        Op("atan2pi","atan2pi","original: atan2(y,x)/\u03c0 with exact-DD \u03c0 division",         fuzz="atan2pi"),
    ]),
    Section("Hyperbolic", [
        Op("sinh",  "sinh",  "original: Taylor series (\\|x\\|<0.1) or (exp\u2212exp\u207b\u00b9)/2",  fuzz="sinh"),
        Op("cosh",  "cosh",  "original: (exp+exp\u207b\u00b9)/2",                                     fuzz="cosh"),
        Op("sinhcosh", "sinhcosh", "fused sinh/cosh: one range-reduction, two outputs",
            fuzz=[("sinh", "sinhcosh_s"), ("cosh", "sinhcosh_c")]),
        Op("tanh",  "tanh",  "original: sinh/cosh (\\|x\\|<0.5) or (1\u2212e\u207b\u00b2\u02e3)/(1+e\u207b\u00b2\u02e3)", fuzz="tanh"),
        Op("asinh", "asinh", "original: Taylor series (\\|x\\|<0.01) or log(x+\u221a(x\u00b2+1)) with Newton", fuzz="asinh"),
        Op("acosh", "acosh", "original: log(x+\u221a(x\u00b2\u22121)) with Newton correction",        fuzz="acosh"),
        Op("atanh", "atanh", "original: Taylor series (\\|x\\|<0.01) or \u00bd\u00b7log((1+x)/(1\u2212x))", fuzz="atanh"),
    ]),
    Section("Error / special functions", [
        Op("erf",           "erf",           "piecewise rational approx (ported from libquadmath erfq.c)",        fuzz="erf"),
        Op("erfc",          "erfc",          "piecewise rational approx + split exp(-x^2)",                       fuzz="erfc"),
        Op("erfcx",         "erfcx",         "exp(x\u00b2)\u00b7erfc(x); scaled form avoiding tail cancellation", fuzz="erfcx"),
        Op("tgamma",        "tgamma",        "piecewise rational approx + Stirling + reflection, exp(lgamma)",    fuzz="gamma"),
        Op("lgamma",        "lgamma",        "piecewise rational approx + Stirling asymptotic",                   fuzz="lngamma"),
        Op("bj0",           "bj0",           "piecewise rational + Hankel asymptotic (j0q.c)",                    fuzz="bj0"),
        Op("bj1",           "bj1",           "piecewise rational + Hankel asymptotic (j1q.c)",                    fuzz="bj1"),
        Op("bjn",           "bjn(3,.)",      "forward/backward recurrence from j0/j1",                            fuzz="bjn"),
        Op("by0",           "by0",           "piecewise rational + Hankel asymptotic (j0q.c)",                    fuzz="by0"),
        Op("by1",           "by1",           "piecewise rational + Hankel asymptotic (j1q.c)",                    fuzz="by1"),
        Op("byn",           "byn(3,.)",      "forward recurrence from y0/y1",                                     fuzz="byn"),
        Op("yn_range(0..5)","yn\\_range(0..5)", "single forward-recurrence sweep, 6 outputs / call",              fuzz="yn_range"),
    ]),
    Section("Complex arithmetic", [
        Op("cdd_add",   "cdd\\_add",   "original: component-wise DD add",
            fuzz=[("re", "cadd_re"), ("im", "cadd_im")]),
        Op("cdd_sub",   "cdd\\_sub",   "original: component-wise DD sub",
            fuzz=[("re", "csub_re"), ("im", "csub_im")]),
        Op("cdd_mul",   "cdd\\_mul",   "original: (ac\u2212bd, ad+bc) via DD ops",
            fuzz=[("re", "cmul_re"), ("im", "cmul_im")]),
        Op("cdd_div",   "cdd\\_div",   "original: (ac+bd, bc\u2212ad)/(c\u00b2+d\u00b2)",
            fuzz=[("re", "cdiv_re"), ("im", "cdiv_im")]),
        Op("cdd_conjg", "cdd\\_conjg", "original: negate im limbs",                fuzz="exact"),
        Op("cdd_proj",  "cdd\\_proj",  "C99 Annex G Riemann-sphere projection (identity for finite z)",
            fuzz=[("re", "cproj_re"), ("im", "cproj_im")]),
        Op("cdd_abs",   "cdd\\_abs",   "original: hypot(re, im)",                  fuzz="cabs"),
        Op("cdd_arg",   "cdd\\_arg",   "original: atan2(im, re)",                  fuzz="carg"),
    ]),
    Section("Complex transcendentals", [
        Op("cdd_sqrt",  "cdd\\_sqrt",  "original: Kahan-style (\\|z\\|+\\|a\\|)/2 with scaling",
            fuzz=[("re", "csqrt_re"), ("im", "csqrt_im")]),
        Op("cdd_exp",   "cdd\\_exp",   "original: exp(re)\u00b7(cos(im), sin(im))",
            fuzz=[("re", "cexp_re"), ("im", "cexp_im")]),
        Op("cdd_expm1", "cdd\\_expm1", "original: expm1(re)\u00b7cos(im) + (cos(im)\u22121) + i\u00b7exp(re)\u00b7sin(im)",
            fuzz=[("re", "cexpm1_re"), ("im", "cexpm1_im")]),
        Op("cdd_log",   "cdd\\_log",   "original: (log(\\|z\\|), atan2(im,re))",
            fuzz=[("re", "clog_re"), ("im", "clog_im")]),
        Op("cdd_log2",  "cdd\\_log2",  "original: clog / log(2) component-wise",
            fuzz=[("re", "clog2_re"), ("im", "clog2_im")]),
        Op("cdd_log10", "cdd\\_log10", "original: clog / log(10) component-wise",
            fuzz=[("re", "clog10_re"), ("im", "clog10_im")]),
        Op("cdd_log1p", "cdd\\_log1p", "original: cancellation-safe log(1+z) near z=0",
            fuzz=[("re", "clog1p_re"), ("im", "clog1p_im")]),
        Op("cdd_pow",   "cdd\\_pow",   "original: exp(w\u00b7log(z))",
            fuzz=[("re", "cpow_re"), ("im", "cpow_im")]),
        Op("cdd_sin",   "cdd\\_sin",   "original: sin(re)cosh(im), cos(re)sinh(im)",
            fuzz=[("re", "csin_re"), ("im", "csin_im")]),
        Op("cdd_cos",   "cdd\\_cos",   "original: cos(re)cosh(im), \u2212sin(re)sinh(im)",
            fuzz=[("re", "ccos_re"), ("im", "ccos_im")]),
        Op("cdd_tan",   "cdd\\_tan",   "original: complex sin/cos ratio",
            fuzz=[("re", "ctan_re"), ("im", "ctan_im")]),
        Op("cdd_sinpi", "cdd\\_sinpi", "original: csin(\u03c0\u00b7z) via π-scaled trig kernels",
            fuzz=[("re", "csinpi_re"), ("im", "csinpi_im")]),
        Op("cdd_cospi", "cdd\\_cospi", "original: ccos(\u03c0\u00b7z) via π-scaled trig kernels",
            fuzz=[("re", "ccospi_re"), ("im", "ccospi_im")]),
        Op("cdd_sinh",  "cdd\\_sinh",  "original: sinh(re)cos(im), cosh(re)sin(im)",
            fuzz=[("re", "csinh_re"), ("im", "csinh_im")]),
        Op("cdd_cosh",  "cdd\\_cosh",  "original: cosh(re)cos(im), sinh(re)sin(im)",
            fuzz=[("re", "ccosh_re"), ("im", "ccosh_im")]),
        Op("cdd_tanh",  "cdd\\_tanh",  "original: complex tanh via sinh/cosh",
            fuzz=[("re", "ctanh_re"), ("im", "ctanh_im")]),
        Op("cdd_asin",  "cdd\\_asin",  "original: \u2212i\u00b7log(iz+\u221a(1\u2212z\u00b2))",
            fuzz=[("re", "casin_re"), ("im", "casin_im")]),
        Op("cdd_acos",  "cdd\\_acos",  "original: \u03c0/2 \u2212 asin(z)",
            fuzz=[("re", "cacos_re"), ("im", "cacos_im")]),
        Op("cdd_atan",  "cdd\\_atan",  "original: (\u2212i/2)\u00b7log((1+iz)/(1\u2212iz))",
            fuzz=[("re", "catan_re"), ("im", "catan_im")]),
        Op("cdd_asinh", "cdd\\_asinh", "original: log(z+\u221a(z\u00b2+1))",
            fuzz=[("re", "casinh_re"), ("im", "casinh_im")]),
        Op("cdd_acosh", "cdd\\_acosh", "original: log(z+\u221a(z\u00b2\u22121))",
            fuzz=[("re", "cacosh_re"), ("im", "cacosh_im")]),
        Op("cdd_atanh", "cdd\\_atanh", "original: \u00bd\u00b7log((1+z)/(1\u2212z))",
            fuzz=[("re", "catanh_re"), ("im", "catanh_im")]),
    ]),
]


@dataclass
class UnifiedRow:
    """One row of the single-arch BENCHMARK.md table.

    Built by ``unified_rows()`` by walking ``C_SECTIONS`` and
    ``FORTRAN_SECTIONS`` and merging entries that share a ``bench_key``.
    For shared ops (most), the C and Fortran ``Op``s carry identical
    ``approach``/``fuzz``; ``prec`` is only set on the Fortran side, so
    we copy it from there.
    """
    display: str
    scope: str           # "C" | "Fortran" | "both"
    prec: str            # "" when only the C-side op is present
    approach: str
    fuzz: FuzzSpec       # for err_dd / err_qp lookup
    c_bench_key: str | None  # None iff scope == "Fortran"
    f_bench_key: str | None  # None iff scope == "C"


def unified_rows() -> list[UnifiedRow]:
    """Flat row list for the single-arch BENCHMARK.md table.

    Walks ``C_SECTIONS`` first (in declaration order), merging each op
    with its same-``bench_key`` counterpart in ``FORTRAN_SECTIONS`` if
    one exists; then appends Fortran-only ops in their original section
    order. No reordering — the table mirrors the section sequence in
    ``ops.py``.
    """
    f_remaining: dict[str, Op] = {
        op.bench_key: op
        for sect in FORTRAN_SECTIONS
        for op in sect.ops
    }
    out: list[UnifiedRow] = []

    for sect in C_SECTIONS:
        for c_op in sect.ops:
            f_op = f_remaining.pop(c_op.bench_key, None)
            if f_op is not None:
                out.append(UnifiedRow(
                    display=f_op.display,
                    scope="both",
                    prec=f_op.prec,
                    approach=f_op.approach,
                    fuzz=f_op.fuzz,
                    c_bench_key=c_op.bench_key,
                    f_bench_key=f_op.bench_key,
                ))
            else:
                out.append(UnifiedRow(
                    display=c_op.display,
                    scope="C",
                    prec=c_op.prec,
                    approach=c_op.approach,
                    fuzz=c_op.fuzz,
                    c_bench_key=c_op.bench_key,
                    f_bench_key=None,
                ))

    for sect in FORTRAN_SECTIONS:
        for f_op in sect.ops:
            if f_op.bench_key not in f_remaining:
                continue
            del f_remaining[f_op.bench_key]
            out.append(UnifiedRow(
                display=f_op.display,
                scope="Fortran",
                prec=f_op.prec,
                approach=f_op.approach,
                fuzz=f_op.fuzz,
                c_bench_key=None,
                f_bench_key=f_op.bench_key,
            ))

    return out
