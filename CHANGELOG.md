# Changelog

All notable changes to IGLOO are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Combustion: `d-smoke`, a constant of both burning models in the models file (refused in `[IGLOO-Models]`, on the
  material line and below 0). A material with `d-smoke > 0` gets five rows of the mass source for its smoke: the
  weight W of what a burning parcel burns in a cell (the metal consumed and the oxide its cap gained), W u and W T,
  u and T weighted by the metal burnt over the cell residence. A coupler reads them to give the smoke its birth
  state; `source.tec` keeps the net rows only.

### Changed
- Combustion: the `combustion` row of a burning material in the models file must write `q-comb`, the heat of
  reaction per kg burnt (`q-comb=0` is allowed); a row without it stops the run.
- Combustion: `combustion=direct` without `T-burn` in its row takes `[IGLOO-Properties] boiling-temperature` as its
  burning temperature only when that is given at one atmosphere; with another `boiling-pressure` the run stops and
  asks for `T-burn`.
- Combustion: a value of the models file is one finite real number; `NaN`, `Inf` and a value holding `,` or `/`
  stop the run.
- Combustion: a `cap=<law>` word on a material line needs its `cap` row in the models file; without it the run stops.
- Combustion: `q-comb`, `d-stop` and `d-burnout` given in `[IGLOO-Models]` stop the run by name; they are constants of
  the material's row in the models file.

### Fixed
- Combustion: the burnout removal of a burning parcel (ODE models 5 and 8) requires the parcel to burn, at or above
  `T-ign` under `Beckstead` and `T-burn` under `direct`; a parcel at the removal mass that does not burn stays an
  inert parcel instead of being handed to the gas. The evaporating models 2 and 4 are unchanged.

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