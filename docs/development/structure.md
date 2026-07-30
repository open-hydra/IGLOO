# Code Structure

---

## Library/executable split

`src/lib/` compiles to a static library `IGLOOL`; `src/app/IGLOO.f90` is an
eight-line driver that links against it:

```fortran
program IGLOO
  use IGLOO_module, only: obj_IGLOO
  implicit none
  type(obj_IGLOO) :: IGLOOsolver

  call IGLOOsolver%setup()
  call IGLOOsolver%solve()
  call IGLOOsolver%writeout()
end program IGLOO
```

All simulation state lives in the single `obj_IGLOO` instance.

---

## Top-level object: `obj_IGLOO`

Defined in `src/lib/obj_IGLOO.f90` (`module IGLOO_module`):

| Field | Type | Role |
|-------|------|------|
| `geoblock(:)` | `obj_block` | Geometry arrays (nodes, face normals, cell volumes) |
| `gasblock(:)` | `obj_flowblock` | Frozen gas field interpolated onto the mesh |
| `source(:)` | `obj_sourceblock` | Gas-coupling source terms (mass, momentum, energy) |
| `euler(:,:)` | `obj_eulerblock` | Equivalent Eulerian fields; 2-D: `(block, family)` |
| `material(:)` | `obj_material` | Per-material data; each `obj_material` owns a `group(:)` array of particle groups |
| `eulSwitch`, `srcSwitch` | `logical` | Enable Eulerian/source output paths |

Methods:

| Method | Purpose |
|--------|---------|
| `setup([external_gas])` | Read config and gas field, pin particles at injection points |
| `solve()` | OpenMP-parallel sweep of all particles |
| `writeout()` | Write Eulerian and source fields to `OUTPUT/` |
| `getSourceTerms(...)` | `pure` function; returns drag force and heat-transfer rate for given local state (used by hydra) |

---

## Hydra integration hook

`obj_IGLOO%setup` accepts an *optional* `external_gas` argument of
`type(orion_data)` (from the ORION library). When present, IGLOO copies the gas
field from memory via `copyORION` and skips reading a Tecplot file; `gas-file` in
`[IGLOO-General]` is ignored.

This is the **only** entry point hydra calls. IGLOO does not participate in the
gas time loop; it integrates all particles in a single sweep over the frozen field.

---

## Hot path: `Lib_Integration.f90::integrate`

Called from `obj_IGLOO%solve` inside `!$OMP PARALLEL DO`. For each particle it:

1. Selects one of four ODE RHS functions via `Lib_RHS.f90::determineModel(phaseChange, brkupEqOde)`:

   | `phaseChange` | `brkupEqOde` | RHS |
   |:---:|:---:|-----|
   | `.false.` | `.false.` | `rhsStandard` |
   | `.true.`  | `.false.` | `rhsEvaporation` |
   | `.false.` | `.true.`  | `rhsBreakupOnly` |
   | `.true.`  | `.true.`  | `rhsEvapBreakup` |

2. Advances the ODE state via `oslo::Run_ODESolver`.
3. Tracks geometric cell crossings via `IGLOO_RayFaceIntersection3D` (ray/face
   intersection on hexahedral cells).
4. On breakup events, records child-particle state in `gr%child(ip)`. The outer
   `maxLoop` iteration in `solve` compacts scattered children and grows the
   particle array geometrically (2×) as needed.

---

## Module naming

Most modules carry the `IGLOO_` prefix: `IGLOO_module`, `IGLOO_particles`,
`IGLOO_data_block`, `IGLOO_IO`, `IGLOO_IO_INI`, `IGLOO_bcBox`, etc.

The integration core is deliberately unprefixed: `Lib_Integration`, `Lib_RHS`,
`Lib_Equations`.

---

## Shared runtime state: `IGLOO_variables`

`src/lib/variables.f90` declares all shared module-level variables: `nb`
(number of blocks), `nm` (number of materials), `eulerSwitch`, `phaseChange`,
`rtol`, `atol`, `ode_word`, tolerances, output unit numbers, body-acceleration
vector, mollification flags, and more. Populated once at startup by
`IGLOO_IO_INI::read_IGLOO_input` reading `input.ini`.

!!! warning "OpenMP and module-level state"
    Any module-level scratch used inside the OMP parallel region must be declared
    `!$omp threadprivate`. See [Contributing](contributing.md) for the rules.

---

## CMake layout

| Target | Links | Notes |
|--------|-------|-------|
| `IGLOOL` (static lib) | `FiNeR::FiNeR OSlo ORION` | `src/lib/*.f90` + `src/lib/config/*.f90` |
| `IGLOO` (exe) | `IGLOOL FiNeR::FiNeR` (+ `-pthread` when `USE_TECIO=ON`) | `src/app/IGLOO.f90` |
| `DocGen` (exe) | `IGLOOL` | `src/app/docgen.f90`; standalone builds only (see below) |

---

## Two build modes

| Mode | Dependency source | Switch |
|------|-------------------|--------|
| `--master=None` | In-tree `lib/{ORION,OSlo,third_party/FiNeR}` | Standalone |
| `--master=hydra` | `$HYDRADIR/lib/{ORION,OSlo,third_party/FiNeR}` | Submodule of hydra |

When IGLOO is configured *by* hydra, only `src/lib/CMakeLists.txt` runs — the
top-level `if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)` block is
skipped and `src/app/CMakeLists.txt` (which adds the `DocGen` target) is not
included. DocGen is therefore a **standalone-only** tool.

!!! warning "Stale lib/ after mode switch"
    Switching between `--master=None` and `--master=hydra` without a full
    `./install.sh build` leaves a stale `lib/` checkout that the build silently
    ignores. Always run `build`, not `compile`, when changing the master setting.

---

## DocGen

`src/app/docgen.f90` is a doc-only executable: it calls
`src/lib/config/Register_IGLOO.f90::Register_IGLOO_Params()` and then
`reg%generate_markdown(...)` to write `docs/user/registry.md`.

The registry in `src/lib/config/Register_IGLOO.f90` **must be kept in sync**
with the runtime reader `src/lib/Lib_INI.f90`. When adding or changing an
`input.ini` key, update both files; see [Contributing](contributing.md).
