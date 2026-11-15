# Versioning

IGLOO requires its libraries as MOSE and ICE do: ORION and OSLO by a range of releases, FiNeR at
one release, and every one is checked at configure time, so that a library IGLOO was not tested
with stops the configure before anything is built.

## How IGLOO requires its libraries

`CMakeLists.txt` names the releases IGLOO is tested with:

| Library | Requirement | How |
|---|---|---|
| ORION | `IGLOO_ORION_VERSION 2.2`: 2.2.0 and later, up to 3.0.0 excluded | `semver_dependency` of `cmake/SemVerDeps.cmake` |
| OSLO | `IGLOO_OSLO_VERSION 1.1`: 1.1.0 and later, up to 2.0.0 excluded | `semver_dependency` of `cmake/SemVerDeps.cmake` |
| FiNeR | `IGLOO_FINER_TAG v2.3.0` | the submodule `lib/FiNeR`, checked against its `VERSION` file |

`cmake/SemVerDeps.cmake` is FLINT's file, copied unmodified: update it from FLINT, never here.
A requirement accepts the releases from its own up to the next change of its leftmost nonzero
number (`2.2` accepts 2.2.0 to 3.0.0 excluded, `0.3` accepts 0.3.x only). Only `vX.Y.Z` tags are
releases; a tag such as `v2.2.0-rc1` is never selected.

IGLOO requires the oldest release that has everything it uses. ORION `2.2`: the property table
is read with `tec_read_structured_multiblock`, whose Tecplot reader 2.1.0 fixed and 2.2.0 made
refuse data that a zone header does not announce. OSLO `1.1`: `run_odesolver` with the
Jacobian arguments IGLOO passes.

FiNeR's own dependencies are fetched as FiNeR's author does (`fobis fetch`: their default
branches, `IGLOO_FINER_DEPS=FREE`). If an upstream change breaks IGLOO, `IGLOO_FINER_DEPS=BACKUP`
checks out the releases of `IGLOO_FINER_BACKUP_PINS`, the ones ATLAS, MOSE and ICE share.

## The checks

| Library | What stops the configure |
|---|---|
| `lib/FiNeR` (submodule) | checked out at another commit than the gitlink of this IGLOO commit (after a `git checkout` without `git submodule update`), or its `VERSION` differs from `IGLOO_FINER_TAG` |
| FiNeR clone (`--include-finer`) | its `VERSION` differs from `IGLOO_FINER_TAG`; with `BACKUP`, a dependency not at its release |
| ORION, OSLO clone (`--include-orion`, `--include-oslo`) | its newest `vX.Y.Z` tag is not in the required range |

To build against a work-in-progress checkout on purpose (an ORION or FiNeR branch under test),
configure with `-DIGLOO_PIN_CHECK=WARNING`: each check then warns instead.

## Which release a build gets

The release of ORION and of OSLO is, in order:

1. the clone given with `-DFETCHCONTENT_SOURCE_DIR_<NAME>=<folder>` (`install.sh --include-orion`,
   `--include-oslo`), as it is; the configure stops if it is not a release that IGLOO accepts;
2. the release in `deps.lock`, if it satisfies the requirement;
3. the newest release on GitHub that satisfies it.

The configure checks the releases out in `lib/ORION` and `lib/OSLO` and writes them to
`deps.lock`; all are git ignored. `-DSEMVER_UPDATE=ALL` (`install.sh build --update-deps`)
selects the newest ones again.

## The requirements in hydra

Hydra builds IGLOO with its own ORION, OSLO and FiNeR. When hydra selects ORION or OSLO with
SemVerDeps, IGLOO's requirement is checked against that release, so a new major ORION or OSLO
in IGLOO's requirement, or another FiNeR release, must move in hydra too.

## Bumping a library

- ORION or OSLO: raise `IGLOO_<NAME>_VERSION` in `CMakeLists.txt`, configure with
  `--update-deps`, run the suite.
- FiNeR: check the new release out in `lib/FiNeR`, `git add lib/FiNeR`, set `IGLOO_FINER_TAG` to
  it in the same commit, run the suite.
