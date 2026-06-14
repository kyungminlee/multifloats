# multifloats documentation

User documentation, built with **Sphinx + Breathe + MyST** and the **Furo**
theme. Doxygen parses the C/C++ header into XML; Breathe pulls that into
Sphinx; the narrative pages are authored in Markdown (MyST).

## One-time setup

Python toolchain (via [uv](https://docs.astral.sh/uv/)):

```sh
uv venv doc/.venv
uv pip install --python doc/.venv -r doc/requirements.txt
```

Plus the **doxygen** binary as a system package:

```sh
sudo apt-get install -y doxygen      # Debian/Ubuntu
# brew install doxygen               # macOS
```

## Build

```sh
./doc/build.sh
# or, from doc/:  make html
```

Output lands in `doc/_build/html/index.html`.

## Layout

```
doc/
  conf.py            Sphinx configuration
  Doxyfile           Doxygen config (header -> _doxygen/xml, XML only)
  index.md           landing page + toctree
  guides/            narrative user guides (MyST Markdown)
  api/               Breathe-extracted C/C++ API reference
  requirements.txt   Python toolchain (sphinx, breathe, myst-parser, furo, ...)
  build.sh, Makefile build drivers
  developer/         developer-only notes (excluded from the built site)
```

## Notes

- The Fortran API is not yet in the site; it is planned via
  [FORD](https://forddocs.readthedocs.io/) and noted on the landing page.
- The header has no `@brief` comments yet, so API entries currently show
  signatures only. Adding `///` / `@brief` comments to
  `include/multifloats/float64x2.h` enriches the reference automatically.
