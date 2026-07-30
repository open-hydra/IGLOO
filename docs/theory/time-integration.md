# Time Integration

Each particle is integrated by `src/lib/Lib_Integration.f90::integrate`, called from
`obj_IGLOO%solve` inside an `!$OMP PARALLEL DO` over the particle array.  The function is
recursive-free: it advances a single particle from injection to domain exit in a single
call, with no inter-particle communication.

---

## Model selection

`determineModel` in `src/lib/Lib_RHS.f90` maps two boolean flags to one of four RHS
routines:

| `phaseChange` | `brkupEqOde` | Model | RHS | State dim |
| :---: | :---: | :---: | :--- | :---: |
| false | false | 1 | `rhsStandard` | 7 (+ 5 euler) |
| true | false | 2 | `rhsEvaporation` | 8 (+ 5 euler) |
| false | true | 3 | `rhsBreakupOnly` | 8 (+ 5 euler) |
| true | true | 4 | `rhsEvapBreakup` | 9 (+ 5 euler) |

The euler extra equations (arc-length $\ell$ and momentum integrand $\mathbf{v}_p\ell$)
are appended only when `eulerSwitch = .true.`; they do not change the model index.
See [Governing equations](governing-equations.md) for the full state vector.

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
2. **`deltat`** budgeted from `computeDeltat(vert)` (cell diagonal / |v|); capped by the
   Stokes drag relaxation time $\tau_p \cdot C_\tau$ (`tauFactor = 100`).
3. `t2 = tStart + 20 · deltat` — generous XEND horizon; the solver stops early at cell
   crossings or events.
4. **Inner loop** (`maxInnerIter = 10`): call `ODEsystem` repeatedly until `exitLoop`; if
   the limit is reached the particle is flagged `gone`.
5. Source / Euler accumulation.
6. Cell-tracking update (`IamOut`, `newGas`); BC handling.
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
- On interrupt (`IamOut` or `newGas` or `eventFlag`): `deltat` is shrunk to
  `max(deltat · min(1/nStep, din/dout), dtMin)` using the Möller-Trumbore face-distance
  `din` to target the boundary.

---

## `solout` callback

The OSlo solver calls `solout` at every accepted step.  It performs:

1. **Scatter cloud** — accumulates the npdot-weight of the accepted interval; emits a
   scatter marker whenever the running total exceeds `dNscat`.
2. **Cell-crossing detection** — `isPointInsideCell(y(1:3), vert, ...)` sets `IamOut`;
   for `ord2` mode a separate gas-dual-cell test sets `newGas`.
3. **NaN / tiny-step guard** — `exitLoop` is set if `any(y/=y)` or `deltat < dtMin`.
4. **Event breakup** — `breakupEvent` (TAB/ETAB/KHRT) detects and applies discrete
   breakup events; `addChild` flag triggers child-parcel spawning.
5. **State snapshot** — on a clean interior step, updates `oldLocal`, `oldStLocal`,
   `oldEvLocal`.

`IRTRN = -2724` on any interrupt; the solver returns immediately.

---

## NaN guards

| Location | Guard |
| :--- | :--- |
| `solout::exitLoop` | `any(y /= y)` or `deltat /= deltat` |
| `ODEsystem` | `any(y(1:6) /= y(1:6))` after solver return (DBGNANSOLV) |
| post-crossing | `any(part%stateVar(1:6) /= ...)` (DBGNANBC) |
| `rhsBreakupOnly`, `rhsEvapBreakup` | npdot or mass ≤ 0 → penalty return (SDIRK4 step rejection) |

!!! note "Stuck-particle detection"
    If a particle stays in the same cell for more than `nMaxCell = 10` outer iterations
    without a crossing (`part%Ncell > 10`), it is flagged `gone` with a diagnostic
    `DBGSTUCK` message.

---

## Parameters reference

| INI key | Default | Effect |
| :--- | :---: | :--- |
| `ode-solver` | `H-sdirk4` | ODE integrator |
| `max-steps-ode` | 100000 | `iopt(1)`; internal step limit per `Run_ODESolver` call |
| `relative-tol` | 1e-10 | `rtol` |
| `absolute-tol` | 1e-10 | `atol` |

Full registry: [../user/registry.md](../user/registry.md).
