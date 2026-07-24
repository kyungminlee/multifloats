# multifloats documentation

- [user/](user/index.md) — using the library: install, link, call from C, C++, and Fortran. Includes the Breathe-extracted [API reference](user/api/index.md).
- [dev/](dev/index.md) — developing multifloats: configure, build, codegen, test, debug, benchmark, release, internals.
- [changelog.md](changelog.md) — include-stub pulling the root `CHANGELOG.md` into the rendered site.

## Building the HTML site

Built with **Sphinx + MyST** (Furo theme) plus **Doxygen + Breathe** for the
C/C++ API; published to GitHub Pages by `.github/workflows/docs.yml` on
pushes to `main`.

One-time setup (via [uv](https://docs.astral.sh/uv/); `doxygen` must be on
the PATH):

```sh
uv venv doc/.venv
uv pip install --python doc/.venv -r doc/requirements.txt
sudo apt-get install -y doxygen      # Debian/Ubuntu
# brew install doxygen               # macOS
```

Then:

```sh
./doc/build.sh
```

Output lands in `doc/_build/html/index.html`. `conf.py` and `Doxyfile` are
**generated** from the `.in` templates by `build.sh` (the version string
comes from the root `VERSION` file) and are git-ignored — edit the `.in`
files, not the generated ones.

## Notes

- The Fortran API is not yet in the site; it is planned via
  [FORD](https://forddocs.readthedocs.io/) and noted on the landing page.
- The header has no `@brief` comments yet, so API entries currently show
  signatures only. Adding `///` / `@brief` comments to
  `include/multifloats/float64x2.h` enriches the reference automatically.
