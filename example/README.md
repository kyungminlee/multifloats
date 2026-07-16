# Examples

Minimal, self-contained programs that compute π to full double-double
precision (`4·atan(1)`) — one via the C++ API, one via the Fortran module.
They double as a template for consuming an installed `multifloats` /
`multifloatsf` package with `find_package()`.

| File | Language | Links |
|---|---|---|
| `pi.cpp` | C++ | `multifloats::multifloats` |
| `pi.f90` | Fortran | `multifloats::fortran_module` |

## Build

The examples build against an **installed** package, not the source tree:

```sh
# From the repo root: build and stage the libraries into a scratch prefix.
cmake -B build -S .
cmake --build build
cmake --install build --prefix /tmp/multifloats-prefix

# Configure and build the examples against that prefix.
cmake -S example -B build/example -DCMAKE_PREFIX_PATH=/tmp/multifloats-prefix
cmake --build build/example

./build/example/pi_cxx    # pi = 3.1415926535897932384626433832795...
./build/example/pi_f
```

CXX is enabled in the example project even for the Fortran program:
`multifloats::fortran_module` links `libmultifloats` (a C++-compiled
archive), and its package selector probes the C++ compiler tag.

See the top-level [`README.md`](../README.md) for the full API surface and
the [developer docs](../doc/dev/) for internals.
