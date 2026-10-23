# Running IGLOO

This page covers the complete workflow for running an IGLOO simulation: preparing the case directory, launching the solver, and inspecting the output. It also describes how IGLOO is embedded in hydra as a one-way particle sub-solver.

---

## The IGLOO Driver

The IGLOO executable (`src/app/IGLOO.f90`) is a short driver:

```fortran
program IGLOO
  use IGLOO_module,   only: obj_IGLOO
  use IGLOO_Mod_MPI,  only: mpi_init_env, mpi_finalize_env
  implicit none
  type(obj_IGLOO) :: IGLOOsolver

  call mpi_init_env()          ! no-op without USE_MPI
  call IGLOOsolver%setup()
  call IGLOOsolver%solve()
  call IGLOOsolver%writeout()
  call mpi_finalize_env()
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
    ├── outloc-<mat>.dat         ← exit or stop location, speed, angle, area, boundary code, diameter per particle
    ├── scatter-<mat>.dat        ← number-density scatter cloud
    ├── snapshot-<mat>.dat       ← each particle's state at time-end (only when time-end is set)
    ├── source.tec               ← gas-coupling source terms on the mesh
    └── euler<fam>.tec           ← equivalent Eulerian fields on the mesh, one per particle family
```

The material name `<mat>` comes from `INPUT/phase.txt` (first column). With a single material called `A`, files are named `trajectories-A.dat`, etc. Families are numbered across materials in phase-file order (one per material × injection group), so a single-group, single-material run writes `euler1.tec`. When `[IGLOO-General] phase` is set, every output name carries the `<phase>-` prefix. A material that solidifies ends every row of its four particle files (trajectories, exit records, scatter cloud, snapshot) with its frozen fraction `f` ([Frozen fraction](#frozen-fraction)).

!!! warning
    IGLOO reads the gas mesh directly from the Tecplot solution file (`gas-file`). There is no separate mesh file; a `MESH/` directory in a case is not used by IGLOO.

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
    Always set `ulimit -s unlimited` (or `KMP_STACKSIZE=100M` for Intel compilers) before launching the solver. Large automatic arrays in the per-particle integration (cell-crossing ray tests and the ODE-solver work arrays) requires more stack than the default system limit.

Output is written to `OUTPUT/` on completion. The solver prints progress to stdout; any fatal errors go to stderr.

### Running a verification case

Every end-to-end case under `tests/` is a runnable simulation directory:

```bash
cd tests/infrastructure/db-2daxi/
rm -rf OUTPUT && mkdir OUTPUT
OMP_NUM_THREADS=4 ../../../bin/IGLOO      # run
python3 check.py                          # oracle gate
```

---

## Output Files

### Trajectories

`OUTPUT/trajectories-<mat>.dat` records particle state at every cell crossing (one row per crossing). Columns (Tecplot ASCII variables header):

```
"X" "Y" "Z" "U" "V" "W" "T" "d_p" "m_p" "ID"
```

Output is controlled by `out-traj` in `[IGLOO-General]` (default: `on`). Rows are written every `print-dcell` cell crossings (default 1) and, when `print-dtime` is positive, also at the first segment end at or past each multiple of that interval (for time rows only, set `print-dcell` above the number of crossings). A parcel that leaves the domain through a boundary face also writes a row at its exit point, whatever the cadence, so its trajectory ends where its exit record is. With `out-time = on` every row carries the particle time [s] as an eleventh column, `"t"`. The rows of a solidifying material end with `"f"`.

!!! warning "Non-deterministic ordering"
    With OpenMP enabled, the row order within a zone is non-deterministic across runs (different thread scheduling). Verification scripts must sort by particle ID and position before comparing, not rely on byte-identical output.

### Scatter cloud

`OUTPUT/scatter-<mat>.dat` is a number-density point cloud: each point represents `dNscat` real droplets, where `dNscat` is auto-sized so the cloud holds approximately `fsample-traj` points per injection stream. Controlled by `out-scatter` in `[IGLOO-General]` (default: `on`).

### Exit locations

`OUTPUT/outloc-<mat>.dat` holds one record per particle, written where its integration stopped: its position, temperature, speed, impact angle, mass flow, face area, and particle ID (`X Y Z T |u_p| alpha mdot Af ID`), followed by the time `t` with `out-time = on`, and then by the boundary code of the face the particle left through and its diameter (`bc d_p`), and by `f` for a solidifying material. `bc` is the face's code in the boundary file (`400`, `301`, `102`, `100`, …); it is `0` when the integration stopped inside the domain (`time-end`, evaporated or burnt out, stuck in a cell or without progress, the outer iteration cap, an ODE-solver failure, a non-finite state, a failure on another MPI rank), and the record is then the particle's state at the stop. A particle that is never integrated, because no cell holds its injection point, records its injection point, speed 0 and t = 0; its temperature and mass flow are the stream's for an assigned-position stream and 0 for a stream injected from a boundary face. With `probe-ids` only the listed particles are integrated and recorded. The mass flow `mdot` depends on the ODE model: for models 1, 3 and 6, whose stream mass changes only at a breakup event, it is the flow the parcel carries at exit; for model 4 (evaporation with ODE breakup) it is the flow at exit `ṅ_p·m`, after evaporation and stripping (a drop consumed in the domain leaves with its remnant); for models 2 and 5, whose drops lose mass to evaporation or combustion, it is the injected flow (a breakup child's birth flow).

### Snapshot at an end time

By default a particle is integrated until it leaves the domain. `time-end` [s] in `[IGLOO-General]` stops every particle at that time instead: a particle still inside gets a last trajectory row at `time-end` exactly, its state goes to `OUTPUT/snapshot-<mat>.dat`, one row per particle with the trajectory columns and the time (`X Y Z U V W T d_p m_p ID t`, then `f` for a solidifying material), and its exit record, code 0, has the same position, speed and time; a particle that leaves earlier writes its exit record as usual. The sources and the Eulerian fields then cover the interval from injection to `time-end`.

### Frozen fraction

A material that solidifies (`solidification=on` on its phase-file line) ends every row of its trajectory, exit, scatter and snapshot files with one more column, `"f"`: the parcel's frozen mass fraction at the row's state, 0 all liquid, 1 all solid, written with six decimals. A material that does not solidify writes no `f`. The files carry no phase column: the phase follows from `f`, the row's temperature `T` and the material's `T-melt` $T_m$:

| Phase | Where |
| :---: | :--- |
| 3, solid | `f = 1` |
| 2, freezing or melting plateau | `0 < f < 1` |
| 1, undercooled liquid | `f = 0` and `T < T_m` |
| 0, liquid | `f = 0` and `T ≥ T_m` |

At six decimals the rule misreads only a plateau state within 5e-7 of `f = 0` or `f = 1`, which is the boundary state itself (all liquid or all solid at $T_m$). The test tools carry the rule as `solid_phase(f, T, Tmelt)` in `tests/tools/datcols.py`.

### Field output

After `solve`, `obj_IGLOO%writeout` writes the mesh-projected fields:

| File | Content |
|------|---------|
| `source.tec` | Gas-coupling source terms: mass, momentum, energy deposition per cell |
| `euler<fam>.tec` | Equivalent Eulerian fields per particle family: density, velocity, temperature, number density (`rho_p u_p v_p w_p T_p n_p`) |

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

When `external_gas` is present, IGLOO copies the gas field directly from memory and skips reading a Tecplot file. The `gas-file` key in `[IGLOO-General]` is then ignored. IGLOO does not participate in the gas time loop: each `solve` is one complete sweep of every particle from injection to exit through the frozen field.

A parent whose gas field evolves runs **repeated sweeps**. `setup` is a wrapper over two calls that the parent then makes separately:

```fortran
call particles%setup_static(external_gas=gas_field)   ! once: INI, mesh, BCs, pinning
do while (coupling)
  call particles%reset_state(external_gas=gas_field)  ! per sweep: re-import the gas, restore the pinned population
  call particles%solve()
  call particles%writeout()
enddo
```

`reset_state` restores every pinned particle to its injection state (same stochastic diameter draw every sweep) and re-imports the gas; it is not optional — `solve` refuses to run on a stale state. The source and eulerian blocks (`particles%source`, `particles%euler`) are valid between a `solve` and the next `reset_state`, which deallocates them. Output files of sweep $N > 0$ carry a `-sweep<N>` suffix; sweep 0 is untagged, so single-sweep runs keep the plain names. Under MPI the parent owns `MPI_Init`/`MPI_Finalize` and must pass `external_gas` identically on every rank (see [Parallelization](../development/parallelization.md)).

After `solve`, hydra can also retrieve the instantaneous gas-coupling terms via `obj_IGLOO%getSourceTerms`, a pure function that evaluates drag force and heat transfer for a given local gas/particle state.

Under `gas-order = 2` the gas is interpolated on a dual mesh whose boundary cells reach one node past every
block face, and IGLOO fills that ring of ghost nodes itself at the start of each `solve` (see
[Boundary conditions](boundary-conditions.md#connected-faces-the-gas-ghost-ring)). A parent that holds better values
for a face — the ghost cells of its own solver — hands them over between `reset_state` and `solve`:

```fortran
call particles%reset_state(external_gas=gas_field)
call particles%set_ghost_face(b, f, density, velocity, temperature, mil, kl, gam, R)   ! per face, per sweep
call particles%solve()
```

The arrays hold the values at the face's ghost nodes over its interior range: `density(nspecies, na, nb)`,
`velocity(3, na, nb)` and the scalars `(na, nb)`, with `(na, nb)` = (Ny, Nz) on faces 1–2, (Nx, Nz) on faces 3–4
and (Nx, Ny) on faces 5–6 (a one-layer block has no ring on faces 5 and 6). They replace IGLOO's own fill on that
face — the linear extrapolation, a wall's mirror, the axis rule and, on a connected face, the interpolated values;
a node with a non-positive scalar falls back to the interior cell's values. The edge and corner nodes beside the
face are extrapolated from the provided values as from IGLOO's own, and on a connected face they are still
interpolated from the connected block unless their other face is provided too. `reset_state` drops every provided
face, so a parent provides them again for each sweep (under MPI identically on every rank, as the gas); a run
that provides none is unchanged.
