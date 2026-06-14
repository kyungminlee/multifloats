# C / C++ API reference

This reference is generated from `include/multifloats/float64x2.h` by Doxygen
and rendered with [Breathe](https://breathe.readthedocs.io/). It covers both
the C ABI (the `float64x2` / `complex64x2` structs and the `extern "C"` `*dd`
entry points) and the C++ `multifloats::float64x2` class API.

```{note}
The header does not yet carry Doxygen-style (`@brief`) comments, so entries
currently show signatures without prose. As `///` / `@brief` comments are
added to the header they will appear here automatically on the next build.
```

The two public types are:

- **`multifloats::float64x2`** — the core type: a struct of two IEEE doubles
  that also carries the full C++ class API (constructors, operators, and the
  `<cmath>`-style member surface). In C it is the plain `float64x2` struct.
- **`multifloats::complex64x2`** — the double-double complex type, an identical
  POD in C and C++.

Below is the complete header surface as extracted by Doxygen — both structs,
the `extern "C"` `*dd` entry points (`muldd`, `sindd`, `logdd`,
`matmuldd_mm`, …), the `<cmath>`-style free functions, and the macros
(`MULTIFLOATS_ABI_VERSION`):

```{eval-rst}
.. doxygenfile:: float64x2.h
   :project: multifloats
```
