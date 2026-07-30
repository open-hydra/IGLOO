---
title: Overview
---

# Overview

IGLOO (**I**ntegration of a **G**eneral **L**agrangian **O**ne-way **O**DE set) is an open-source Lagrangian particle solver written in modern Fortran. It integrates each particle's ODE state — position, velocity, temperature, diameter — through a *steady* background gas field, sweeping every particle from injection to exit of the domain in a single pass. The gas field is either read from a Tecplot file or passed in memory by a parent solver. Drag and heat transfer are always integrated; evaporation, secondary breakup (with optional child-particle generation), and Eulerian feedback are opt-in physics.

---

## Hydra CFD Suite

IGLOO is the **Lagrangian particle solver** of the **Hydra** CFD ecosystem — an integrated suite of tools for multi-physics simulation of complex systems involving compressible flows and dispersed phases.

| Component | Role | Status |
|-----------|------|--------|
| [**ATLAS**](https://github.com/open-hydra/ATLAS) | Pre-processor: mesh preparation, initial and boundary conditions, particle BC definition | Separate package |
| [**MOSE**](https://github.com/open-hydra/MOSE) | Compressible Euler/Navier–Stokes gas solver on multi-block structured grids | Companion solver |
| **IGLOO** | Lagrangian particle solver: integrates dispersed-phase ODE through a steady gas field | This package |

A typical workflow is:

1. **ATLAS** generates the mesh, preprocesses boundary conditions, and writes `INPUT/bc.txt`, `INPUT/phase.txt`, and `INPUT/properties.dat`.
2. **MOSE** (or another CFD solver) computes the steady gas field and writes a Tecplot solution file.
3. **IGLOO** reads that solution as a frozen background, injects particles, and integrates them through the domain.

IGLOO can also run standalone — without MOSE — when the gas field is provided as a Tecplot file directly.

!!! info "Using IGLOO without ATLAS"
    ATLAS is the standard way to produce the BC and material files that IGLOO requires. If ATLAS is not available, all files can be prepared manually; see the [User Guide](user/using.md) for the expected formats.

---

## Capabilities

| Feature | Details |
|---------|---------|
| **Drag** | Stokes, Morsi–Alexander, Crowe, Hermsen, Henderson, Putnam models |
| **Heat transfer** | Kavanau–Drake, Ranz–Marshall Nusselt correlations |
| **Evaporation** | Vapor-pressure driven mass transfer; coupled diameter and temperature evolution |
| **Secondary breakup** | Five models (Pilch–Erdman, Reitz, KHRT, TAB, Beale–Reitz); optional child-particle generation with dynamic array compaction/growth |
| **Injection** | Boundary-patch injection (2D sweep or 3D advancing-front hexagonal packing); assigned position injection; stochastic diameter distributions (Dirac, Normal, LogNormal, Rosin–Rammler) |
| **Eulerian feedback** | Particle statistics deposited onto the gas mesh as Eulerian and source fields, smoothed by a volume-weighted binomial mollifier |
| **ODE integration** | Explicit (DOPRI5) or stiff implicit (H-SDIRK4) stepping via the OSlo library; geometric cell tracking by ray/face intersection |
| **Parallelism** | OpenMP-parallel hot loop over particles |
| **I/O** | Tecplot ASCII/binary via ORION; trajectory, Eulerian-field, scatter-cloud, and exit-location output |

---

## Dependencies

IGLOO is built on three companion libraries, included as Git submodules in a standalone build or shared with the hydra suite when built as a submodule.

| Library | Role |
|---------|------|
| [ORION](https://github.com/MarcoGrossi92/ORION) | Multi-format I/O (Tecplot ASCII/binary) |
| [OSlo](https://github.com/MarcoGrossi92/OSlo) | ODE solver library (DOPRI5, H-SDIRK4) |
| [FiNeR](https://github.com/szaghi/FiNeR) | INI configuration file parser |

Optional: **OpenMP** (parallelism), **TecIO** (Tecplot binary output, pulled in transitively by ORION).
