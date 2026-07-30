# Running IGLOO

This page covers the complete workflow for running an IGLOO simulation: preparing the case directory, launching the solver, and inspecting the output. It also describes how IGLOO is embedded in hydra as a one-way particle sub-solver.

---

## The IGLOO Driver

The IGLOO executable (`src/app/IGLOO.f90`) is an eight-line driver:

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

All simulation state is owned by a single `obj_IGLOO` instance. The three calls correspond to three distinct phases:

| Phase | Method | What it does |
|-------|--------|--------------|
| Setup | `obj_IGLOO%setup` | Reads `input.ini`, loads the gas field, reads phase/BC files, pins particles at injection points. |
| Solve | `obj_IGLOO%solve` | OpenMP-parallel integration loop: sweeps every particle from injection to exit of the domain. |
| Write | `obj_IGLOO%writeout` | Writes Eulerian and source fields to `OUTPUT/`. |

IGLOO does not own a time loop. Each call to `solve` integrates all particles from injection to exit in a single pass over the frozen gas field.

---

## Case Directory Structure

Every IGLOO case follows this layout:

```
my_case/
├── input.ini             ← solver configuration (IGLOO sections + ATLAS BC/material sections)
├── INPUT/
│   ├── solfile.tec       ← frozen gas field (Tecplot ASCII; name set by gas-file in input.ini)
│   ├── bc.txt            ← boundary condition table (numeric codes, ATLAS-generated)
│   ├── phase.txt         ← material name and group count per material
│   └── properties.dat    ← thermodynamic property tables (Tecplot format: cp, rho, h vs T)
└── OUTPUT/               ← created at run time
    ├── trajectories-<mat>.dat   ← per-cell-crossing state (X Y Z U V W T d_p m_p ID)
    ├── outloc-<mat>.dat         ← exit location, speed, angle, area per particle
    ├── scatter-<mat>.dat        ← number-density scatter cloud
    ├── source.tec               ← gas-coupling source terms on the mesh
    └── eulerian-<mat>.tec       ← equivalent Eulerian fields on the mesh
```

The material name `<mat>` comes from `INPUT/phase.txt` (first column). With a single material called `A`, files are named `trajectories-A.dat`, etc.

!!! warning
    IGLOO reads the gas mesh directly from the Tecplot solution file (`gas-file`). There is no separate mesh file; `MESH/` directories seen in legacy test cases are not used by IGLOO.

---

## Running the Solver

```bash
ulimit -s unlimited           # IGLOO requires a large stack
export KMP_STACKSIZE=100M     # Intel compilers: explicit stack for OpenMP threads
export OMP_NUM_THREADS=4      # number of OpenMP threads (default: system decides)

cd my_case/
mkdir -p OUTPUT
/path/to/bin/IGLOO
```

!!! warning "Stack size"
    Always set `ulimit -s unlimited` (or `KMP_STACKSIZE=100M` for Intel compilers) before launching the solver. Deep recursion in the ray/face intersection and ODE integrator requires more stack than the default system limit.

Output is written to `OUTPUT/` on completion. The solver prints progress to stdout; any fatal errors go to stderr.

### Running a verification case

Every end-to-end case under `tests/` is a runnable simulation directory:

```bash
cd tests/infrastructure/db-2daxi/
rm -rf OUTPUT && mkdir OUTPUT
OMP_NUM_THREADS=4 ../../../bin/IGLOO      # run
python3 check.py                          # oracle gate
```

(The legacy `test/` tree with its `IGLOO.sh` wrapper was retired on 2026-07-16;
a copy is parked out-of-git at `~/Desktop/Software/toBeRemoved/IGLOO-legacy-test/`.)

---

## Output Files

### Trajectories

`OUTPUT/trajectories-<mat>.dat` records particle state at every cell crossing (one row per crossing). Columns (Tecplot ASCII variables header):

```
"X" "Y" "Z" "U" "V" "W" "T" "d_p" "m_p" "ID"
```

Output is controlled by `out-traj` in `[IGLOO-General]` (default: `on`). Sampling frequency is set by `fsample-traj` (default: 1 — every crossing; increase to reduce file size).

!!! warning "Non-deterministic ordering"
    With OpenMP enabled, the row order within a zone is non-deterministic across runs (different thread scheduling). Verification scripts must sort by particle ID and position before comparing, not rely on byte-identical output.

### Scatter cloud

`OUTPUT/scatter-<mat>.dat` is a number-density point cloud: each point represents `dNscat` real droplets, where `dNscat` is auto-sized so the cloud holds approximately `fsample-traj` points per injection stream. Controlled by `out-scatter` in `[IGLOO-General]` (default: `on`).

### Exit locations

`OUTPUT/outloc-<mat>.dat` records, for every particle that exits the domain, its exit position, speed, impact angle, face area, and particle ID.

### Field output

After `solve`, `obj_IGLOO%writeout` writes the mesh-projected fields:

| File | Content |
|------|---------|
| `source.tec` | Gas-coupling source terms: mass, momentum, energy deposition per cell |
| `eulerian-<mat>.tec` | Equivalent Eulerian fields: number density, velocity, temperature |

Both outputs are mollified by a local binomial smoother (controlled by `mollify` and `mollify-passes` in `[IGLOO-General]`) to reduce deposition noise.

---

## Embedding in Hydra

When IGLOO is used as a hydra submodule, `obj_IGLOO%setup` accepts an optional `external_gas` argument of `type(orion_data)` (from the ORION library):

```fortran
use IGLOO_module,  only: obj_IGLOO
use Lib_ORION_data, only: orion_data
type(obj_IGLOO)  :: particles
type(orion_data) :: gas_field   ! provided by the hydra gas solver

call particles%setup(external_gas=gas_field)
call particles%solve()
call particles%writeout()
```

When `external_gas` is present, IGLOO copies the gas field directly from memory and skips reading a Tecplot file. The `gas-file` key in `[IGLOO-General]` is then ignored. This is the only entry point hydra calls — IGLOO does not participate in the gas time loop.

After `solve`, hydra can retrieve the gas-coupling source terms via `obj_IGLOO%getSourceTerms`, a pure function that evaluates drag force and heat transfer for a given local gas/particle state.
