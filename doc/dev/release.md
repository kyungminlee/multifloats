# Releasing

How a versioned release is cut and published. The whole pipeline lives in
`.github/workflows/release.yml`; this is the operator's view.

## Versioning

Semantic versioning ([semver](https://semver.org/)). Two independent axes:

- **Release version** (`vMAJOR.MINOR.PATCH`) lives in exactly two places:
  the top entry of `CHANGELOG.md` and the git tag. There is no version
  string in `CMakeLists.txt` or any header — nothing else to bump.
- **ABI version** is `#define MULTIFLOATS_ABI_VERSION` in
  `include/multifloats/float64x2.h`. It is the single source of truth for the
  shared library `SOVERSION` (parsed at configure time in `CMakeLists.txt` and
  `src/CMakeLists.txt`). Bump it **only** on a breaking C-ABI change — a new
  release version does not imply a new ABI version.

Pick the bump by what changed: build/CI/packaging only → PATCH; new
backward-compatible API/behavior → MINOR; breaking ABI or API → MAJOR (and
bump `MULTIFLOATS_ABI_VERSION`).

## Cut a release

1. **Update `CHANGELOG.md`.** Add a `## [X.Y.Z] — YYYY-MM-DD` (ISO-8601
   UTC) entry above the previous one, grouped `Added` / `Changed` / `Fixed` /
   `Removed`. Lead with one line stating whether numeric behavior or the ABI
   changed. Commit to `main` (`docs(changelog): cut X.Y.Z`).

2. **(Optional) Dry-run the matrix** without publishing — see below. Recommended
   when the release touches the build or the workflow itself.

3. **Tag and push:**

   ```sh
   git tag -a vX.Y.Z -m "vX.Y.Z"
   git push origin vX.Y.Z
   ```

   The tag push triggers the workflow. Everything after this is automatic.

## What the workflow does

Trigger: push of a `v*` tag (publishes) or manual `workflow_dispatch`
(dry-run, does **not** publish — the publish step is gated on
`startsWith(github.ref, 'refs/tags/v')`).

Jobs, all gating the release:

- **`build-lto`** — per-compiler LTO matrix: GCC 12–15, LLVM 19–21, Intel
  2023–2025 (Linux x86-64). GCC entries build + **run the full test suite**
  (`BUILD_TESTING=ON`); non-GCC entries build the two libraries only (their C++
  tests need libquadmath / `-fext-numeric-literals`, GCC-only).
- **`build-lto-mac`** — Homebrew GCC 13–15 (macOS arm64); the newest runs the
  suite.
- **`build-non-lto`** / **`build-non-lto-mac`** — the portable compiler-agnostic
  archive + the shared library; run the suite, and `build-non-lto` additionally
  audits the shipped `.so`'s dynamic symbol table (`check_shared_exports.sh`).
- **`release`** — needs all builds; merges artifacts into per-platform combined
  tarballs and publishes a GitHub Release with auto-generated notes.

A test failure or a symbol-audit failure (`check_exported_symbols`,
`check_shared_exports`) fails the build and blocks the release. Before this was
wired up (≤ v0.9.0) the jobs configured with `BUILD_TESTING` off and ran zero
tests — see the v0.9.1 CHANGELOG entry.

## Artifacts

Per release, on the GitHub Release page:

- **`libmultifloats-vX.Y.Z-{linux-x86_64,darwin-arm64}.tar.gz`** — the combined
  per-platform bundle a consumer wants: every per-compiler archive + the
  portable non-LTO archive + shared library, under one install prefix whose
  `multifloatsConfig.cmake` routes `find_package` to the best-fit archive.
- **`multifloats-vX.Y.Z-<platform>-<lto|nolto>-<compiler>.tar.gz`** — the
  individual per-compiler archives that compose the bundle.

## Dry-run

To verify the full matrix (especially macOS) without cutting a tag:

```sh
gh workflow run release.yml --ref <branch>
gh run watch $(gh run list --workflow release.yml --branch <branch> \
  --limit 1 --json databaseId -q '.[0].databaseId') --exit-status
```

It builds and tests everything and uploads artifacts, but skips publishing.
(Tarball names are derived from the ref with `/` sanitized to `-`, so a
slashed branch name packages cleanly — a no-op on real tags.)

## If a release goes wrong

The tag is the release. To redo, fix on `main`, then either move the tag
(`git tag -f`, force-push — only safe if no one has pulled it) or, preferably,
cut the next PATCH. Deleting a published release + its assets is a manual step
on the GitHub Release page.
