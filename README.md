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

- **Drag & heat transfer** — momentum and thermal coupling with the carrier gas through a library of literature drag laws and Nusselt-number correlations.
- **Evaporation** — opt-in phase change with vapor-pressure driven mass transfer and coupled diameter/temperature evolution.
- **Secondary breakup** — five breakup models with optional child-particle generation; the particle array is compacted and grown on the fly.
- **Injection** — assigned particle positions or boundary-patch injection with per-cell spacing, 3D advancing-front packing, and stochastic diameter sampling (Dirac, Normal, LogNormal, Rosin–Rammler).
- **Eulerian feedback** — particle statistics deposited back onto the gas mesh as eulerian and source fields, smoothed by a volume-weighted binomial mollifier.
- **Adaptive time integration** — explicit (DOPRI5) or stiff implicit (SDIRK4) ODE stepping via the OSlo library, with geometric cell tracking by ray/face intersection.
- **Parallel execution** — the hot loop is OpenMP-parallel over particles.
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
```

The executable is placed in `bin/IGLOO`. `./install.sh compile` performs an incremental rebuild from the existing CMake preset.

See the [Installation Guide](https://open-hydra.github.io/IGLOO/getting-started/installation/) for all build options and troubleshooting.

### Run the test suite

```bash
./tests/test.sh all
```

`tests/` is a single model-first tree (`standard/`, `evaporation/`, `breakup/`, `infrastructure/`); each category holds oracle-checked end-to-end cases and literature-grounded model tests, all registered in ctest. The legacy `test/` cases are deprecated and will be removed once `tests/` is complete. See the [Quick Start](https://open-hydra.github.io/IGLOO/getting-started/quick-start/) for a full walkthrough.

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
│   ├── app/           # IGLOO executable (8-line driver) + DocGen
│   └── lib/           # Solver library (IGLOOL)
├── lib/               # Git submodule dependencies (standalone build)
├── tests/             # Verification & validation suite (ctest, model-first)
│   ├── standard/      # Drag + heat: unit families & e2e cases
│   ├── evaporation/   # Evaporation unit family + d²-law e2e
│   ├── breakup/       # TAB, Pilch–Erdman, Reitz–Diwakar, ETAB
│   └── infrastructure/# Gas reconstruction & INI pipeline
├── test/              # Legacy cases (deprecated, superseded by tests/)
├── docs/              # MkDocs documentation source
├── install.sh         # Build helper script
└── CMakeLists.txt
```

## Documentation

Full documentation is available at **[open-hydra.github.io/IGLOO](https://open-hydra.github.io/IGLOO/)**, covering:

- Installation & quick start
- User guide & input file reference
- Theory guide (Lagrangian formulation, drag, heat, evaporation, breakup, mollification)
- Verification & validation cases

## License

IGLOO is free and open-source software released under the [GNU General Public License v3.0](LICENSE).
