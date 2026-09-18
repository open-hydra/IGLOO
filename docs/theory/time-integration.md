# Time Integration

Each particle is integrated by `src/lib/Lib_Integration.f90::integrate`, called from
`obj_IGLOO%solve` inside an `!$OMP PARALLEL DO` over the particle array.  The function is
recursive-free: it advances a single particle from injection to domain exit in a single
call, with no inter-particle communication.

---

## Model selection

`determineModel` in `src/lib/Lib_RHS.f90` maps two boolean flags to one of four RHS
routines (plus the combustion family, model 5):

| `phaseChange` | `brkupEqOde` | Model | RHS | State dim |
| :---: | :---: | :---: | :--- | :---: |
| false | false | 1 | `rhsStandard` | 7 (+ 5 euler) |
| true | false | 2 | `rhsEvaporation` | 8 (+ 6 euler) |
| false | true | 3 | `rhsBreakupOnly` | 8 (+ 6 euler) |
| true | true | 4 | `rhsEvapBreakup` | 9 (+ 7 euler) |

A material with `combustion = Beckstead` is routed to model 5 (`rhsAlCombustion`,
model-2 layout) instead.  The euler extra equations (arc-length $\ell$ and the momentum,
energy and mass moments) are appended only when `eulerSwitch = .true.`; they do not
change the model index.  See [Governing equations](governing-equations.md) for the full
state vector.

---

## ODE solver

The time-stepping backend is the OSlo library, accessed via `oslo::Run_ODESolver`.  Two
adaptive solvers are selectable from `input.ini`:

```
[IGLOO-ODE]
ode-solver   = H-sdirk4    ; default; or H-dopri5
max-steps-ode = 100000     ; iopt(1); default 100000
relative-tol  = 1e-10      ; rtol; default 1e-10
absolute-tol  = 1e-10      ; atol; default 1e-10
```

`H-sdirk4` (L-stable, stiffly accurate, 4th-order DIRK) is the default because the
Stokes drag relaxation $\tau_p = \rho_p d^2/(18\mu_g)$ can be very short relative to
the cell-transit time, making the ODE stiff.  `H-dopri5` (explicit Dormand-Prince 5th
order) is appropriate for low-Re/large-particle problems where stiffness is absent.

---

## Outer loop

```
do while (part%time < tlimit .and. iter < maxIter)   ! maxIter = 500 000
```

Each outer iteration covers one cell transit:

1. `initializeCell` / `computeSource` — record entry state.
2. **`deltat`** budgeted from `computeDeltat(vert)` (cell diagonal / |v|, on the gas dual
   cell when `gas-order = 2`); a non-finite estimate falls back to
   $C_\tau\,\tau_p$ with $\tau_p$ the Stokes drag relaxation time (`tauFactor = 100`).
3. `t2 = tStart + 20 · deltat` — generous XEND horizon; the solver stops early at cell
   crossings or events.
4. **Inner loop** (`maxInnerIter = 10`): call `ODEsystem` repeatedly until `exitLoop`; if
   the limit is reached the particle is flagged `gone`.
5. Source / Euler accumulation.
6. Cell-tracking update (`IamOut`, `newGas`, `sectorOut`); BC handling, including the
   wedge fold when the segment ended on a sector plane.
7. Trajectory and scatter-cloud output.

---

## `ODEsystem` inner procedure

`ODEsystem` wraps a single `Run_ODESolver` call and interprets the result:

- `safety = 0.99` — pre-scales `deltat` before each ODE call to avoid exact-boundary
  roundoff.
- `err < 0` (solver failure, e.g. SDIRK4 B2 null-interval) → particle marked `gone`.
- `exitLoop` is set true when the solver reaches XEND without an interrupt, or when the
  step is infinitesimally small (`dout < eps`), `deltat < dtMin = 1e-14`, or a NaN
  appears in `y`.
- On interrupt (`IamOut`, `newGas`, `sectorOut` or `eventFlag`): `deltat` is shrunk to
  `max(deltat · min(1/nStep, din/dout), dtMin)` using the Möller-Trumbore face-distance
  `din` to target the boundary (`nStep = 10`; for a sector plane `din` is the exact
  point-to-plane distance).  The refined step is also floored at a displacement-based
  `dtSafe` (the time to move a fraction of the containment tolerance `eps` at the current
  speed) so refinement cannot stall on steps that no longer move the particle.

---

## `solout` callback

The OSlo solver calls `solout` at every accepted step.  It performs:

1. **Scatter cloud** — accumulates the npdot-weight of the accepted interval; emits a
   scatter marker whenever the running total exceeds `dNscat`.
2. **Cell-crossing detection** — `isPointInsideCell(y(1:3), vert, ...)` sets `IamOut`
   (and `sectorOut` when the point left the wedge sector); for `ord2` mode a separate
   gas-dual-cell test sets `newGas`.
3. **Burnout test** — for mass-consuming models (evaporation, combustion) the droplet is
   declared consumed once its mass drops to `mBurnTol = 1e-15` kg.
4. **NaN / tiny-step guard** — `exitLoop` is set if `any(y/=y)`, `deltat < dtMin`, or the
   step moved the particle by less than the containment tolerance.
5. **Event breakup** — `breakupEvent` (TAB/ETAB/KHRT) detects and applies discrete
   breakup events; `addChild` flag triggers child-parcel spawning.
6. **State snapshot** — on a clean interior step, updates `oldLocal`, `oldStLocal`,
   `oldEvLocal`.

`IRTRN = -2724` on any interrupt; the solver returns immediately.

---

## Robustness guards

| Location | Guard |
| :--- | :--- |
| outer loop | a non-finite `deltat` estimate falls back to `tauFactor · taup` |
| `solout::exitLoop` | `any(y /= y)` or `deltat /= deltat` interrupts the solver |
| `ODEsystem` | a non-finite state after the solver returns is reverted to the last good accepted step and the particle marked `gone` (a consuming droplet hands its remnant to the gas) |
| `ODEsystem` | a `deltat` update outside `[dtMin, huge)` is a hard error |
| RHS routines (models 2–5) | any non-finite RHS entry (unphysical Newton trial: $m \le 0$, $\dot n_p \le 0$) is replaced by a `1e30` penalty so SDIRK4 rejects the step |

!!! note "Stuck-particle detection"
    A particle that stays in the same cell for more than `nMaxCell = 10` outer
    iterations without a crossing, or makes no net displacement for `nMaxStall = 100`
    iterations, is flagged `gone` with a diagnostic message; a particle that reaches
    `maxIter` is flagged `gone` with its exit record written.

---

## Parameters reference

| INI key | Default | Effect |
| :--- | :---: | :--- |
| `ode-solver` | `H-sdirk4` | ODE integrator |
| `max-steps-ode` | 100000 | `iopt(1)`; internal step limit per `Run_ODESolver` call |
| `relative-tol` | 1e-10 | `rtol` |
| `absolute-tol` | 1e-10 | `atol` |

Full registry: [../user/registry.md](../user/registry.md).
