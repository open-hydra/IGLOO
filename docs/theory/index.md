# Theoretical Guide

IGLOO integrates a set of Lagrangian particle ODEs through a steady background gas field.
Each computational particle is a parcel representing $\dot{n}_p$ real droplets per unit time.
The four governing physics layers — drag, heat transfer, evaporation, and breakup — are
selectable per material group; the active combination determines which ODE system is assembled
at run time.

---

## Model selection

Two boolean flags, set during initialisation from `input.ini`, control which right-hand-side
system is assembled by `determineModel` in `Lib_RHS.f90`:

- **`phaseChange`** — evaporation is active (requires an evaporation model and species data).
- **`brkupEqOde`** — at least one continuous-rate breakup model is active (Pilch-Erdman,
  Reitz-Diawakar, or Reitz-KHRT).  TAB and ETAB do *not* set this flag; they inject a discrete
  breakup event through the `solout` callback while the particle advances under model 1 or 2.

| `phaseChange` | `brkupEqOde` | Model | RHS routine | State dimension |
| :---: | :---: | :---: | :--- | :---: |
| false | false | 1 | `rhsStandard` | 7 (+ 5 euler) |
| true | false | 2 | `rhsEvaporation` | 8 (+ 5 euler) |
| false | true | 3 | `rhsBreakupOnly` | 8 (+ 5 euler) |
| true | true | 4 | `rhsEvapBreakup` | 9 (+ 5 euler) |

A material with `combustion = Beckstead` overrides this matrix and is routed to
model 5, `rhsAlCombustion` (state dimension 8, model-2 layout) — see
[Metal Combustion](combustion.md).

The euler extra equations (arc length $\ell$ and momentum integrand $\mathbf{v}\ell$) are appended
only when `eulerSwitch` is enabled; they do not change the model index.

!!! note "TAB and ETAB use model 1 or 2"
    The Taylor Analogy Breakup (TAB) and Enhanced TAB (ETAB) models treat breakup as a
    discrete event detected in the `solout` callback rather than adding an ODE equation.
    They therefore run under model 1 (no evaporation) or model 2 (with evaporation) and
    do not appear in the `brkupEqOde` column of the table above.

---

## Pages in this guide

<div class="grid cards" markdown>

-   :material-math-integral: **[Governing equations](governing-equations.md)**

    Per-particle state vector, the four RHS systems, and the body-acceleration term.

-   :material-arrow-expand-all: **[Drag](drag.md)**

    The 13 drag-coefficient correlations and their selection keyword.

-   :material-thermometer: **[Heat transfer](heat.md)**

    Nusselt-number correlations and the particle-property lookup tables.

-   :material-water-percent: **[Evaporation](evaporation.md)**

    d²-law, CEM, CEM-B, and Abramzon-Sirignano mass-transfer models.

-   :material-fire: **[Metal combustion](combustion.md)**

    Beckstead d^n aluminum burn law, ignition gate, and heat release.

-   :material-explosion: **[Breakup](breakup.md)**

    Pilch-Erdman, Reitz-Diawakar, Reitz-KHRT, TAB, and ETAB; child-particle appending.

-   :material-spray: **[Injection](injection.md)**

    BC-pinned face injection, advancing-front 3D hex packing, and stochastic diameter sampling.

-   :material-transfer: **[Eulerian feedback](eulerian-feedback.md)**

    Source and Euler field accumulation, Favre averaging, and body-force correction.

-   :material-blur: **[Mollification](mollification.md)**

    Dimension-split binomial smoother for deposited fields.

-   :material-clock-fast: **[Time integration](time-integration.md)**

    The `integrate` loop, OSlo ODE solver interface, geometric cell-crossing sweep, and NaN guards.

-   :material-cube-outline: **[Geometry](geometry.md)**

    Hex face conventions, point-in-cell tests, and ray-polyhedron intersection.

</div>

---

## Summary table

| Page | Key source files |
| :--- | :--- |
| [Governing equations](governing-equations.md) | `Lib_RHS.f90`, `Lib_Equations.f90`, `obj_particles.f90` |
| [Drag](drag.md) | `Lib_Drag.f90` |
| [Heat transfer](heat.md) | `Lib_Heat.f90`, `Lib_Thermodynamics.f90` |
| [Evaporation](evaporation.md) | `Lib_Evaporation.f90`, `Lib_RHS.f90` |
| [Metal combustion](combustion.md) | `Lib_Combustion.f90`, `Lib_RHS.f90` |
| [Breakup](breakup.md) | `Lib_Breakup.f90`, `Lib_Integration.f90` |
| [Injection](injection.md) | `initialization.f90`, `Lib_Statistics.f90` |
| [Eulerian feedback](eulerian-feedback.md) | `Lib_Integration.f90`, `obj_gas.f90`, `obj_condensed.f90`, `obj_block.f90` |
| [Mollification](mollification.md) | `Lib_Mollify.f90`, `obj_block.f90` |
| [Time integration](time-integration.md) | `Lib_Integration.f90` |
| [Geometry](geometry.md) | `geometry.f90`, `obj_block.f90` |
