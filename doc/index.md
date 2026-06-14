# multifloats

**Double-double arithmetic for Fortran, C and C++.**

`multifloats` provides two derived/template types — `real64x2` / `float64x2`
(~106-bit, ~32 decimal digits) and a complex counterpart — implemented
natively from error-free transformations on regular IEEE doubles. There are
no quad-precision (`real(16)` / `__float128`) temporaries in any arithmetic,
transcendental, or array kernel.

| Fortran type | C / C++ type  | Storage         | Precision                      |
| ------------ | ------------- | --------------- | ------------------------------ |
| `real64x2`   | `float64x2`   | 2 × `real(dp)`  | ~106 bits (~32 decimal digits) |
| `cmplx64x2`  | `complex64x2` | 2 × `float64x2` | ~106 bits per component        |

The C++ class API mirrors `<cmath>`; the Fortran `real64x2` interface mirrors
`REAL(KIND=16)` so existing quad-precision code can switch with minimal
changes; and a stable C ABI (the `*dd` entry points) bridges both.

```{toctree}
:maxdepth: 2
:caption: Guides

guides/getting-started
guides/building
guides/precision
guides/matmul
guides/error-handling
```

```{toctree}
:maxdepth: 2
:caption: Reference

api/index
BENCHMARK
CHANGELOG
```

## Quick taste

::::{tab-set}
:::{tab-item} C++
```cpp
#include <multifloats/float64x2.h>
using namespace multifloats;
float64x2 x(1.0);
float64x2 pi = float64x2(4.0) * atan(x);
std::cout << to_string(pi, 32) << "\n";   // 32 digits
```
:::
:::{tab-item} C
```c
#include <multifloats/float64x2.h>   /* -lmultifloats */
float64x2 x = {1.0, 0.0};
float64x2 pi4 = atandd(x);            /* pi/4 as a DD */
float64x2 pi  = muldd((float64x2){4.0, 0.0}, pi4);
```
:::
:::{tab-item} Fortran
```fortran
use multifloats
integer, parameter :: dp = 8
type(real64x2) :: x, y
x = real64x2(1.0_dp)
y = sqrt(atan(x) * 4)   ! sqrt(pi) to full DD
print *, y
```
:::
::::

```{note}
The Fortran API reference is not yet part of this site. It will be added
later via [FORD](https://forddocs.readthedocs.io/). For now the Fortran
surface is documented narratively in the guides and in the project README.
```
