# Performance Notes

---

## OpenMP parallelism

The hot loop in `obj_IGLOO%solve` is `!$OMP PARALLEL DO SCHEDULE(DYNAMIC)` over
particles within each material group. Particles are independent — no communication
between them — so the parallel efficiency is limited only by load imbalance
(particles with different trajectory lengths) and the overhead of the cell-search
and ODE stepping per particle.

Set thread count via `OMP_NUM_THREADS` (and `KMP_STACKSIZE=100M` for Intel
compilers):

```bash
export OMP_NUM_THREADS=8
export KMP_STACKSIZE=100M
ulimit -s unlimited
```

---

## Per-particle cost drivers

| Cost driver | Notes |
|-------------|-------|
| **ODE tolerance** | `relative-tol`/`absolute-tol` in `[IGLOO-ODE]`; tighter tolerances increase step count and cost. |
| **Cell-crossing search** | `searchInBlock` in `obj_particles.f90`: guided walk from the previous cell, O(1) per crossing in regular grids. Cost rises in highly distorted meshes. |
| **Breakup child generation** | Adds outer `maxLoop` iterations over `solve`; each loop re-enters the OMP region. |
| **Mollifier passes** | `mollify-passes` in `[IGLOO-General]`; default 8. Each pass is O(N_cells). Increase sparingly. |
| **Scatter cloud** | `fsample-traj` in `[IGLOO-General]` sets the points per injection stream; the cloud is the largest output stream on breakup-heavy cases. |

---

## Cell locator: `searchInBlock`

`searchInBlock` (in `src/lib/obj_particles.f90`) is a guided walk: it seeds the
search from the particle's previous cell and steps neighbour-by-neighbour using
the face-intersection logic. In practice this is near-O(1) for particles with
smooth trajectories. The locator is called at every cell crossing, not every ODE
step.

---

## ODE solver choice: DOPRI5 vs H-SDIRK4

Set via `ode-solver` in `[IGLOO-ODE]`.

| Solver | Type | Use when |
|--------|------|----------|
| `H-dopri5` | Explicit Dormand–Prince 5(4) | Non-stiff; fast per step; typical drag/heat cases with large particles |
| `H-sdirk4` | Stiff implicit (default) | Evaporation or cases with very small particle inertia; more expensive per step but stable with large `dt` |

For stiff cases (e.g. small droplets at high gas temperature), SDIRK4 will take
fewer but costlier steps and converge where DOPRI5 would hit `maxInnerIter`.

!!! warning "Verify that particles advance"
    A solver stall (step size collapsing towards `dtMin`) ends a particle through the
    stuck-particle guards with a diagnostic line rather than a crash, and a run can
    finish "cleanly" with every particle ended at injection. When touching the
    integration path or the boundary conditions, check the trajectory files — and the
    oracle-gated cases in `tests/` — not just the exit status.

---

## Stack requirements

IGLOO requires a large stack:

```bash
ulimit -s unlimited        # system limit
export KMP_STACKSIZE=100M  # Intel OpenMP per-thread stack
```

The ray/face intersection logic (`IGLOO_RayFaceIntersection3D`) and the ODE
integrator use significant local arrays. Without these settings, silent stack
overflows may produce incorrect results rather than a crash.
