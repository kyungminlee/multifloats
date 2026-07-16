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
  conf.py            Sphinx configuration (site root = user/index)
  Doxyfile           Doxygen config (header -> _doxygen/xml, XML only)
  user/              the published user manual (MyST Markdown)
    index.md         landing page + toctree
    api/             Breathe-extracted C/C++ API reference
  dev/               contributor docs (excluded from the built site)
  CHANGELOG.md       release history (linked from user/index)
  requirements.txt   Python toolchain (sphinx, breathe, myst-parser, furo, ...)
  build.sh, Makefile build drivers
```

## Notes

- The Fortran API is not yet in the site; it is planned via
  [FORD](https://forddocs.readthedocs.io/) and noted on the landing page.
- The header has no `@brief` comments yet, so API entries currently show
  signatures only. Adding `///` / `@brief` comments to
  `include/multifloats/float64x2.h` enriches the reference automatically.
