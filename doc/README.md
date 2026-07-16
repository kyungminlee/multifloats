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
  conf.py.in         Sphinx config template (@MULTIFLOATS_VERSION@ from VERSION)
  Doxyfile.in        Doxygen config template (header -> _doxygen/xml, XML only)
  user/              the published user manual (MyST Markdown)
    index.md         landing page + toctree (site root = user/index)
    api/             Breathe-extracted C/C++ API reference
  dev/               contributor docs (excluded from the built site)
  changelog.md       stub that {include}s the root CHANGELOG.md into the site
  _static/           static assets copied verbatim into the site
  _templates/        Jinja template overrides for the theme
  requirements.txt   Python toolchain (sphinx, breathe, myst-parser, furo, ...)
  build.sh, Makefile build drivers
```

`conf.py` and `Doxyfile` are **generated** from the `.in` templates by
`build.sh` / `make` (the version string is read from the root `VERSION` file)
and are git-ignored — edit the `.in` files, not the generated ones.

## Notes

- The Fortran API is not yet in the site; it is planned via
  [FORD](https://forddocs.readthedocs.io/) and noted on the landing page.
- The header has no `@brief` comments yet, so API entries currently show
  signatures only. Adding `///` / `@brief` comments to
  `include/multifloats/float64x2.h` enriches the reference automatically.
