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

---

## Build modes

### Standalone (`--master=None`)

IGLOO is built against the in-tree copies of ORION, OSlo, and FiNeR under `lib/`. This is the standard mode for users who do not have the hydra suite installed.

```bash
# Intel compilers with OpenMP (recommended for production)
./install.sh build --master=None --compilers=intel --use-openmp

# GNU compilers with OpenMP
./install.sh build --master=None --compilers=gnu --use-openmp

# With optional Tecplot binary I/O (requires a C++ compiler)
./install.sh build --master=None --compilers=intel --use-openmp --use-tecio
```

### Hydra-submodule (`--master=hydra`)

When IGLOO is a submodule of hydra, the build reuses `$HYDRADIR/lib/{ORION,OSlo,third_party/FiNeR}` instead of the in-tree copies. The `HYDRADIR` environment variable must be set (see the hydra `bootstrap.sh`).

```bash
./install.sh build --master=hydra --compilers=intel --use-openmp
```

!!! warning "Do not mix build modes"
    Switching between `--master=None` and `--master=hydra` without running a full `build` leaves a stale `lib/` checkout that the build silently ignores. Always run `./install.sh build` (not `compile`) when changing the master setting.

---

## Build options

| Flag | Description |
|------|-------------|
| `--master=None` or `--master=hydra` | Required. Selects the dependency source. |
| `--compilers=intel` or `--compilers=gnu` | Selects the compiler family. When omitted, CMake decides. |
| `--use-openmp` | Enables OpenMP parallelization. |
| `--use-tecio` | Enables Tecplot binary I/O (requires C++ compiler). |

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
