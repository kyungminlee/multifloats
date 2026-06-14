# Getting started

This page gets you from a fresh checkout to a working `float64x2` program in
each supported language. For the full set of build options see
[Building](building.md).

## Install / build the libraries

```sh
cmake -B build -S .
cmake --build build
```

A default build produces only the two installable libraries:

- **`libmultifloats.a`** — the C / C++ kernels. The API is header-only via
  `include/multifloats/float64x2.h`; this archive holds the out-of-line math
  bodies and the `extern "C"` `*dd` entry points.
- **`libmultifloatsf-<compiler>.a`** — the Fortran module library.

Both are exported as CMake packages. Prefer `find_package(multifloats)` /
`find_package(multifloatsf)` over hard-coded `-l` flags so the
compiler-tagged archive name is resolved automatically.

## C++

```cpp
#include <multifloats/float64x2.h>   // header-only class API
using namespace multifloats;

int main() {
    float64x2 x(1.0);
    float64x2 pi = float64x2(4.0) * atan(x);
    std::cout << to_string(pi, 32) << "\n";   // scientific, 32 digits
}
```

Compile and link against `libmultifloats.a`:

```sh
c++ -std=c++17 main.cc -lmultifloats
```

The full `<cmath>` double-name surface (`sqrt`, `exp`, `log`, `sin`, `atan2`,
`erf`, `fma`, `hypot`, the `isnan`/`signbit` classifiers, …) is available via
argument-dependent lookup on `float64x2`. See the [API reference](../api/index.md).

## C (via the C ABI)

```c
#include <multifloats/float64x2.h>   /* -lmultifloats */
#include <stdio.h>

int main(void) {
    float64x2 x   = {1.0, 0.0};
    float64x2 pi4 = atandd(x);                       /* pi/4 */
    float64x2 pi  = muldd((float64x2){4.0, 0.0}, pi4);
    printf("pi.hi = %.17g  pi.lo = %.17g\n", pi.limbs[0], pi.limbs[1]);
}
```

The C functions follow the libc / libquadmath naming convention: a type
suffix on the name (`sindd`, `logdd`, `cexpdd`) mirroring `sinf` / `sinq`.

## Fortran

```fortran
program demo
    use multifloats          ! generic sqrt/atan overloads live here
    integer, parameter :: dp = 8
    type(real64x2) :: x, y
    x = real64x2(1.0_dp)     ! from a dp literal
    y = sqrt(atan(x) * 4)    ! sqrt(pi) to full DD
    print *, y               ! defined I/O: ~32 digits
end program
```

Link against `libmultifloatsf-<compiler>.a`. The `real64x2` interface is
designed to mirror `REAL(KIND=16)`, so porting existing quad-precision code is
mostly a matter of changing the type declaration.

## Next steps

- [Building](building.md) — every CMake option, tests, optional comparison harnesses.
- [Precision](precision.md) — what is full double-double and what is not, and why.
- [Error handling](error-handling.md) — the NaN-in-NaN-out policy.
