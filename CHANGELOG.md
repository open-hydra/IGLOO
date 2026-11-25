# Changelog

All notable changes to IGLOO are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Coupling: `IGLOO_Mod_Exchange` (`src/lib/coupling/`), the one module a parent code needs to exchange with IGLOO:
  `exchange_rate` and `gas_rate` (what the parcels give the gas, deposit and absorb), `track_processed` (the state the
  processed mass leaves at), `add_absorption_rate` and `absorption_packet` (mass offered to the parcels, from the gas
  or any other phase), `add_birth_rate` (new parcels of a 0-group line), `clear_import_rates`, `material_id`. A
  refusal comes back in plain words. See docs/development/coupling.md.
- Coupling: `add_removal_rate` takes mass off the parcels of a material in a cell, for the gas or any other phase:
  whole droplets leave at their own state, in proportion to the mass the parcels held there in the last sweep, and
  `exchange_rate` reports what left (`removed`). Every model can lose droplets this way.
- Combustion: `d-gen` in a burning material's `combustion` row, the diameter a coupler gives the drops of the mass
  the material generates in another phase (`gen_diameter`); IGLOO only holds it.
- Coupling: the parcels of a material absorb the mass offered in a cell in proportion to the mass they held there in
  the last sweep, mixing its momentum and energy; what they absorbed is reported per cell, apart from the smoothed
  source fields.
- `[IGLOO-General]`: `max-reflections` (default 10000) drops a parcel after that many reflections; `graze-standoff`
  and `corner-clear` set where a parcel stands after a grazing impact on a symmetry face (300). All three are read
  again at every sweep.

### Changed
- Combustion: the `combustion` row of a burning material in the models file must write `q-comb`, the heat of
  reaction per kg burnt (`q-comb=0` is allowed); a row without it stops the run.
- Combustion: `combustion=direct` without `T-burn` in its row takes `[IGLOO-Properties] boiling-temperature` as its
  burning temperature only when that is given at one atmosphere; with another `boiling-pressure` the run stops and
  asks for `T-burn`.
- Combustion: a value of the models file is one finite real number; `NaN`, `Inf` and a value holding `,` or `/`
  stop the run.
- Combustion: the oxide cap is a generic **carried material**: the line word `carry=constant|Gallier21` and the row
  `<material> carry <law> share= yield= carried= carried-0= residue= d-residue=` replace `cap=`, `xi-cap`, `xi-ox0`
  and `product`. The yield (kg formed per kg burnt, 1.889441 for aluminium and alumina) is an input, no longer built
  in; a carrying material without it stops the run. A `carry=<law>` word needs its `carry` row.
- A line of the phase file with 0 groups holds born parcels only (burnout residues, a coupler's births) and always has
  one family, named by `carried=` or not.
- Combustion: `q-comb`, `d-stop`, `d-burnout` and `d-gen` given in `[IGLOO-Models]` stop the run by name; they are constants of
  the material's row in the models file.
- Breakup: the constants of a material's breakup model are its row of the models file,
  `<material> breakup <model> key=value ...`, plain keys (`Comega`, `B1`, ...), each optional (absent: the published
  default; `WeCrit` as written, stored doubled). The 20 `[IGLOO-Models]` constants (`Cd` ... `WeTrans`) and the 23
  `<tag>-<key>` tokens of the material line (`TAB-Comega`, ...) stop the run by name; a row key its model does not
  read stops it too. `[IGLOO-Models] breakup` and the line's `breakup=<model>` are unchanged.
- Reflection: a parcel that grazes a symmetry face (below 2 % incidence) is set back 1 % of its diameter
  (`graze-standoff`), no longer a fixed 1e-6 m; with `corner-clear` it stands off every symmetry face of its cell.
- Docs: the user guide is split by topic (Run Settings, Materials & Models, Material Properties, Integrator, Wall
  Erosion), the theory pages hold the physics only, a Code Map serves developers, and Getting Started lists the
  whole Hydra suite.

### Fixed
- Geometry: a cell's volume is the integral of det J of the trilinear map of its nodes, the volume its bilinear faces
  bound, as MOSE and ICE compute it; it was a sum of five tetrahedra, exact only on plane faces, whose volumes did not
  tile a mesh with twisted faces (a 6×6×6 block with its inner nodes moved by up to 5 % of a cell summed to
  1 + 3.7e-5 of its box, a cell 3.9e-2 off). A cell's eight octants (gas order 2) add up to its volume, and the
  chimera overlap's 24 tetrahedra enclose the same volume. Plane-faced meshes change by round-off.
- Combustion: the burnout removal of a burning parcel (ODE models 5 and 8) requires the parcel to burn, at or above
  `T-ign` under `Beckstead` and `T-burn` under `direct`; a parcel at the removal mass that does not burn stays an
  inert parcel instead of being handed to the gas. The evaporating models 2 and 4 are unchanged.
- Combustion: a `combustion=direct` drop at its removal mass is removed from the temperature its burning branch
  burns from, `T-burn` within 1e-10 (`directOnset`, shared by both tests). A drop whose step into `T-burn` ended in
  that margin burned on past `d-burnout` to a solver failure instead of being removed (18 of the 25 drops of
  `tests/combustion/burnout-direct`).
- Injection (3-D): a parcel pinned on a warped injecting face is moved inside its cell; where the face bulged
  outwards it was never injected (2.4 % of the injecting area of a star grain).
- Breakup at gas order 2 with `H-sdirk4`: the first step of a parcel read the gas before it was set, so a drop could
  get a non-finite Weber number and never break up.

### Removed
- The Monolith timing harness (`tests/bench/monolith/`) and the `IGLOO_IO_PROBE` configure option: the harness, its
  I/O probe and its measured-performance table are now the igloo-bench repository, which builds the probe through
  a CMake project hook.

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