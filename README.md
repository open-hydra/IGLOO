<p align="center">
  <h1 align="center">IGLOO</h1>
  <p align="center"><b>Integration of a General Lagrangian One-way ODE set</b></p>
</p>

<p align="center">
  <a href="https://open-hydra.github.io/IGLOO/"><img src="https://img.shields.io/badge/docs-online-brightgreen.svg" alt="Documentation"></a>
  <img src="https://img.shields.io/badge/language-Fortran-734f96.svg" alt="Language: Fortran">
  <a href="https://github.com/open-hydra/IGLOO/blob/main/LICENSE"><img src="https://img.shields.io/badge/license-GPLv3-blue.svg" alt="License: GPLv3"></a>
</p>

---

IGLOO is an open-source Lagrangian particle solver written in modern Fortran. It integrates each particle's ODE state — position, velocity, temperature, diameter — through a *steady* background gas field, from injection to exit of the domain. The gas field is either read from a Tecplot file or injected in memory by a parent solver (the hydra suite). Drag and heat transfer are always integrated; evaporation, secondary breakup (with optional child particles), and eulerian feedback are opt-in physics.

## Features

- **Drag & heat transfer** — momentum and thermal coupling with the carrier gas through 13 literature drag laws (incompressible and compressible) and six Nusselt-number correlations.
- **Evaporation** — opt-in phase change (d²-law, CEM, CEM-B, Abramzon–Sirignano, Tonini–Cossali) with equilibrium or Langmuir–Knudsen interface and coupled diameter/temperature evolution.
- **Metal combustion** — Beckstead d^n aluminium burn law with ignition gate, heat release and burnout.
- **Secondary breakup** — five breakup models (Pilch–Erdman, Reitz–Diwakar, Reitz-KHRT, TAB, ETAB) with KHRT child-particle generation; the particle array is compacted and grown on the fly.
- **Injection** — assigned particle positions or boundary-patch injection with per-cell spacing, 3D advancing-front packing, and stochastic diameter sampling (Dirac, Normal, LogNormal, Rosin–Rammler).
- **Eulerian feedback** — particle statistics deposited back onto the gas mesh as eulerian and source fields, smoothed by a volume-weighted binomial mollifier.
- **Adaptive time integration** — explicit (DOPRI5) or stiff implicit (SDIRK4) ODE stepping via the OSlo library, with geometric cell tracking by ray/face intersection.
- **Parallel execution** — the hot loop is OpenMP-parallel over particles; an optional hybrid MPI + OpenMP build distributes particles over ranks with rank-count-invariant output.
- **Flexible I/O** — Tecplot input/output via ORION; trajectory, eulerian-field, and scatter-cloud output.

## Quick Start

### Prerequisites

| Requirement | Details |
|---|---|
| **CMake** | ≥ 3.23 |
| **Fortran compiler** | GNU (`gfortran`) or Intel/oneAPI (`ifort` / `ifx`) |
| **C/C++ compiler** | Required only for optional TecIO support |

### Build

```bash
git clone --recurse-submodules https://github.com/open-hydra/IGLOO.git
cd IGLOO

# Standalone build with GNU compilers and OpenMP
./install.sh build --master=None --compilers=gnu --use-openmp

# — or with Intel compilers and Tecplot binary I/O —
./install.sh build --master=None --compilers=intel --use-openmp --use-tecio

# As a hydra submodule (reuses $HYDRADIR's dependency tree)
./install.sh build --master=hydra --compilers=intel --use-openmp

# Hybrid MPI + OpenMP (not combinable with --use-tecio)
./install.sh build --master=None --compilers=intel --use-openmp --use-mpi
```

The executable is placed in `bin/IGLOO`. `./install.sh compile` performs an incremental rebuild from the existing CMake preset.

See the [Installation Guide](https://open-hydra.github.io/IGLOO/getting-started/installation/) for all build options and troubleshooting.

### Run the test suite

```bash
./tests/test.sh all
```

`tests/` is a single model-first tree (`standard/`, `evaporation/`, `combustion/`, `breakup/`, `infrastructure/`, plus `repeatability/` and `mpi/`); each category holds oracle-checked end-to-end cases and literature-grounded model tests, all registered in ctest. See the [Quick Start](https://open-hydra.github.io/IGLOO/getting-started/quick-start/) for a full walkthrough.

## Dependencies

IGLOO is built on top of companion libraries, included as Git submodules (standalone build) or shared with the hydra suite (`--master=hydra`):

| Library | Role |
|---|---|
| [ORION](https://github.com/MarcoGrossi92/ORION) | Multi-format I/O (Tecplot ASCII/binary) |
| [OSlo](https://github.com/MarcoGrossi92/OSlo) | ODE solver library (DOPRI5, SDIRK4) |
| [FiNeR](https://github.com/szaghi/FiNeR) | INI configuration file parser |

Optional external libraries: **OpenMP**, **TecIO** (pulled in transitively by ORION).

## Project Structure

```
IGLOO/
├── src/
│   ├── app/           # IGLOO executable (driver) + DocGen
│   └── lib/           # Solver library (IGLOOL)
├── lib/               # Git submodule dependencies (standalone build)
├── tests/             # Verification & validation suite (ctest, model-first)
│   ├── standard/      # Drag + heat: unit families & e2e cases
│   ├── evaporation/   # Evaporation unit families + d²-law, LK, TC and MHB98 e2e
│   ├── combustion/    # Beckstead burn-law unit family + burn-box e2e
│   ├── breakup/       # TAB, ETAB, Pilch–Erdman, Reitz–Diwakar, Reitz-KHRT
│   ├── infrastructure/# Gas reconstruction, INI pipeline, injection/BC/wedge cases, refusals
│   ├── repeatability/ # Two-sweep state-leak gates
│   └── mpi/           # Rank-count gates (USE_MPI builds)
├── docs/              # MkDocs documentation source
├── install.sh         # Build helper script
└── CMakeLists.txt
```

## Documentation

Full documentation is available at **[open-hydra.github.io/IGLOO](https://open-hydra.github.io/IGLOO/)**, covering:

- Installation & quick start
- User guide & input file reference
- Theory guide (Lagrangian formulation, drag, heat, evaporation, combustion, breakup, injection, Eulerian feedback, mollification, time integration, geometry)
- Verification & validation cases

## License

IGLOO is free and open-source software released under the [GNU General Public License v3.0](LICENSE).
