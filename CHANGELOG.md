# Changelog

All notable changes to IGLOO are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.0.0] - 2026-11-16

### Added
- Initial public release of IGLOO, an open-source Lagrangian particle solver written in modern Fortran.
- Lagrangian particle tracking through a steady background gas field, integrating particle position, velocity, temperature, and diameter.
- Aerodynamic drag and heat transfer models, including multiple drag laws, Nusselt-number correlations, and optional radiative heat exchange.
- Evaporation models, including the d²-law, CEM, CEM-B, Abramzon–Sirignano, and Tonini–Cossali formulations.
- Aluminium particle combustion based on the Beckstead burn law.
- Secondary breakup models: Pilch–Erdman, Reitz–Diwakar, Reitz-KHRT, TAB, and ETAB, including child-particle generation.
- Particle injection at assigned positions and through boundary patches, with support for stochastic particle-size distributions.
- Eulerian feedback through the deposition of particle statistics and source fields onto the gas mesh.
- Adaptive ODE integration using explicit DOPRI5 and stiff implicit SDIRK4 schemes.
- Parallel particle tracking using OpenMP, with optional hybrid MPI and OpenMP execution.
- Tecplot-based input/output, including particle trajectories, Eulerian fields, and scatter-cloud data.
- Build and installation scripts, project documentation, and a verification and validation test suite.