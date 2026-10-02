# Installation

This page describes how to build IGLOO. Two build modes exist: **standalone** (IGLOO owns its own copies of ORION, OSlo, and FiNeR) and **hydra-submodule** (IGLOO reuses the dependency tree already present in `$HYDRADIR/lib/`).

---

## Prerequisites

| Requirement | Details |
|-------------|---------|
| **CMake** | ≥ 3.23 |
| **Fortran compiler** | GNU (`gfortran`) or Intel/oneAPI (`ifx`) |
| **C/C++ compiler** | Required only for optional TecIO support |
| **OpenMP** | Optional; required for multi-threaded execution |
| **MPI** | Optional; hybrid MPI + OpenMP execution (`--use-mpi`) |

---

## Build modes

### Standalone (default)

IGLOO is built against the in-tree copies of ORION, OSlo, and FiNeR under `lib/`. This is the standard mode for users who do not have the hydra suite installed.

A plain `git clone https://github.com/open-hydra/IGLOO.git` is enough: `install.sh build` initialises the `lib/ORION`, `lib/OSLO` and `lib/third_party/FiNeR` submodules, and the CMake configure clones FiNeR's own dependencies, all over https (a bare `cmake -B build` does the same). OSLO's nested SUNDIALS submodule is not fetched, since IGLOO builds OSLO without SUNDIALS.

```bash
# Intel compilers with OpenMP (recommended for production)
./install.sh build --compilers=intel --use-openmp

# GNU compilers with OpenMP
./install.sh build --compilers=gnu --use-openmp

# With optional Tecplot binary I/O (requires a C++ compiler)
./install.sh build --compilers=intel --use-openmp --use-tecio
```

### Hydra-submodule (`--include-*`)

When IGLOO is a submodule of hydra, point the build at `$HYDRADIR/lib/{ORION,OSLO,third_party/FiNeR}` instead of the in-tree copies:

```bash
./install.sh build --compilers=intel --use-openmp \
  --include-orion=$HYDRADIR/lib/ORION \
  --include-oslo=$HYDRADIR/lib/OSLO \
  --include-finer=$HYDRADIR/lib/third_party/FiNeR
```

The three options are independent: any dependency left out comes from the in-tree copy. A relative path is taken from the repository root.

!!! warning "`compile` keeps the dependency paths of the last `build`"
    `./install.sh build` records the three paths in `CMakePresets.json`, and `compile` reuses them. Always run `./install.sh build` (not `compile`) when changing a dependency path.

---

## Build options

| Flag | Description |
|------|-------------|
| `--include-orion=<path>` | ORION checkout to build against. Default: the in-tree `lib/ORION`. |
| `--include-oslo=<path>` | OSLO checkout to build against. Default: the in-tree `lib/OSLO`. |
| `--include-finer=<path>` | FiNeR checkout to build against. Default: the in-tree `lib/third_party/FiNeR`. |
| `--compilers=intel` or `--compilers=gnu` | Selects the compiler family. When omitted, CMake decides. |
| `--use-openmp` | Enables OpenMP parallelization. |
| `--use-tecio` | Enables Tecplot binary I/O (requires C++ compiler). |
| `--use-mpi` | Enables the hybrid MPI + OpenMP build (see [Parallelization](../development/parallelization.md)). Cannot be combined with `--use-tecio`. |

The `build` command wipes `build/`, runs a clean CMake configure and build, then writes `CMakePresets.json` from the populated `CMakeCache.txt`. The executable is placed at `bin/IGLOO`.

---

## Incremental rebuild

After an initial `build`, use `compile` for fast iteration when only source files have changed:

```bash
./install.sh compile
```

`compile` runs `cmake --preset default && cmake --build build` without wiping the build directory. It depends on a valid `CMakePresets.json` written by the last `build` run.

!!! note
    `./install.sh build` **wipes** `build/` — never use it for incremental work. Use `compile` between source-only changes.

---

## CMake presets

After a successful `build`, `CMakePresets.json` at the repo root records the compiler paths and cache variables used. You can rebuild directly with CMake:

```bash
cmake --preset default
cmake --build build
```

This is equivalent to `./install.sh compile`.

---

## Next steps

- **[Quick Start](quick-start.md)** — run the first verification case and confirm the solver is working.
