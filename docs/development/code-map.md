# Code Map

The [Theory Guide](../theory/index.md) describes the physics and the algorithms in words and equations. This page
says where each one lives in the code, under the names the code uses. Paths are under `src/lib/`.

---

## Where each topic lives

| Topic | Files |
| :--- | :--- |
| [Governing equations](../theory/governing-equations.md) | `Lib_RHS.f90`, `Lib_Equations.f90`, `obj_particles.f90` |
| [Drag](../theory/drag.md) | `Lib_Drag.f90` |
| [Heat transfer](../theory/heat.md) | `Lib_Heat.f90`, `Lib_Thermodynamics.f90` |
| [Evaporation](../theory/evaporation.md) | `Lib_Evaporation.f90`, `Lib_RHS.f90` |
| [Metal combustion](../theory/combustion.md) | `Lib_Combustion.f90`, `Lib_RHS.f90` |
| [Solidification](../theory/solidification.md) | `Lib_Solidification.f90`, `Lib_RHS.f90` |
| [Breakup](../theory/breakup.md) | `Lib_Breakup.f90`, `Lib_Integration.f90` |
| [Wall erosion](../theory/erosion.md) | `Lib_Erosion.f90` (laws, input checks, wall sums), `Lib_Integration.f90::depositErosion`, `IO.f90::write_wallfield` |
| [Injection](../theory/injection.md) | `initialization.f90`, `Lib_Statistics.f90` |
| [Eulerian feedback](../theory/eulerian-feedback.md) | `Lib_Integration.f90`, `obj_block.f90`, `obj_condensed.f90` |
| [Mollification](../theory/mollification.md) | `Lib_Mollify.f90`, `obj_block.f90` |
| [Time integration](../theory/time-integration.md) | `Lib_Integration.f90` |
| [Geometry](../theory/geometry.md) | `geometry.f90`, `obj_particles.f90`, `obj_bc.f90` |
| [Coupling](coupling.md) | `coupling/IGLOO_Mod_Exchange.f90` |

---

## The ODE systems

`Lib_RHS.f90::determineModel` maps a material's flags to its system, the `model` index of the code. `phaseChange`
is set by an evaporation model; `brkupEqOde` by a continuous-rate breakup model (Pilch–Erdman, Reitz–Diwakar,
Reitz-KHRT). TAB and ETAB set `brkupEvent` only, and act through the `solout` callback.

| Model | Selected by | Right-hand side | State size (+ Euler slots) |
| :---: | :--- | :--- | :---: |
| 1 | neither flag | `rhsStandard` | 7 (+5) |
| 2 | `phaseChange` | `rhsEvaporation` | 8 (+6) |
| 3 | `brkupEqOde` | `rhsBreakupOnly` | 8 (+6) |
| 4 | both | `rhsEvapBreakup` | 9 (+7) |
| 5 | `combustion=Beckstead` or `direct` | `rhsBurn` | 8 (+6) |
| 6 | `solidification=on` | `rhsSolidification` | 8 (+5) |
| 7 | solidification and `brkupEqOde` | `rhsSolidification` | 9 (+6) |
| 8 | Beckstead with a carried material (`share > 0` or `carry=Gallier21`) | `rhsBurnCarry` | 9 (+6) |

### State layout (`stateVar`)

| Index | Content | Models |
| :---: | :--- | :--- |
| 1–3 | position | all |
| 4–6 | velocity | all |
| 7 | temperature, or enthalpy with a $c_p(T)$ table (not under solidification) | all |
| 8 | drop mass (the metal's under model 8) | 2, 4, 5, 8 |
| 8 | number rate | 3 |
| 8 | frozen fraction | 6, 7 |
| 9 | number rate | 4, 7 |
| 9 | carried mass $m_\mathrm{ox}$ | 8 |

With `eulerSwitch` (`out-file` holds `E`) the moments follow, at `nOde+k`:

| Offset | Content |
| :---: | :--- |
| +1 | arc length $\int|\mathbf v_p|$ |
| +2..4 | momentum integrand, weighted by $m$ (models 2, 5, 8) or $m\,\dot n_p$ (model 4), unweighted otherwise; meridian frame on a wedge |
| +5 | energy integrand ($T_p\ell$ or $h_p\ell$), same weights |
| +6 | mass moment $m$ (2, 5, 8), number moment $\dot n_p$ (3, 7), $m\,\dot n_p$ (4) |
| +7 | number moment (model 4 only) |

The body-force accumulators $J$, $W$ (`srcBodyForce`) sit at the tail for models 2, 4, 5, 8 with a body force and
source output; $J$ reuses the mass moment when both are on. `obj_particle` keeps beside the state `m`, `d`, `npdot`,
`tp`, `mdot`, the geometry cell `i(4)` and the gas dual cell `igas(3)`.

---

## Shared kernels

- `Lib_Equations.f90::interphase`: the drag force and the convective heat of every model. `cpFactor` is $1/c_p$ for
  the temperature state and 1 for the enthalpy state; it turns the heat into the rate of index 7.
- `Lib_Drag.f90::assign_drag`: the drag word to the correlation; each correlation guards its denominators with the
  local `toll = 1e-20`.
- `Lib_Heat.f90::assign_heat`, `heat(Re, Pr, Ma, heatSelect)`; `heat_Kavanau_Drake`; `radiativeHeat` with the CODATA
  constant `sigmaSB`.
- `obj_IGLOO.f90::getSourceTerms`: the drag and heat per unit volume of a cell's particle state, for hydra's Eulerian
  coupling (`IGLOO-iter = 2`).
- `Lib_Equations.f90::interp2ndOrder` (3D) and `sampleGas2D` (2D, axisymmetric): the second-order gas sample;
  `toMeridian` rotates the momentum deposit on a wedge. `obj_block.f90::fillGhostGradient` and `fillGhostPartners`
  fill the ghost ring.

### Evaporation

`Lib_Evaporation.f90::evaporation` computes the shared quantities and dispatches. The species properties come from
the array `ep(nep)`:

| Key | Slot |
| :--- | :--- |
| `Mv` | `ep(iMv)` |
| `Lv` | `ep(iLv)` |
| `cpv` | `ep(icpv)` |
| `Le` | `ep(iLe)` |
| `Yinf` | `ep(iYinf)` |
| `boiling-temperature` | `ep(iinvTboil)` = $1/T_\mathrm{boil}$ |
| Clausius–Clapeyron slope | `ep(iLvMvOverRu)` = $L_v M_v / R_u$ |

ASM and TC set `override_Qdot = .true.` to replace the convective heat of `interphase`. `lkCorrection` solves the
Langmuir–Knudsen fixed point; a parcel whose bracket fails is counted into
`[WARNING] LK-interface: material <name>: the Langmuir-Knudsen fixed point not bracketed in N parcel(s)`.

### Combustion

`Lib_Combustion.f90` holds the Beckstead closure (the oxidizer exponent is the parameter `aXeff`), `directRate` and
the onset factor `directOnset` of the direct law, the lobe table (`buildLobeTable`, `lobeF`: 513 nodes in $\alpha$,
bisection on $f_\mathrm{ox}$, linear in $(1 - f_\mathrm{ox})^{1/3}$) and the stretches: `carryStretch`,
`carryStretchX0` (initial oxide), `carryLawPath` (Gallier's law, per parcel). The yield is the auxiliary slot `imYld`.
The digitized points of Gallier's figures and their provenance are in `tests/combustion/burn-box-cap-gallier/reference/`.

### Breakup

`Lib_Breakup.f90::assign_breakup` sets `brkupEqOde`, `brkupEvent` and `brkupHasChild`; `breakupOde` is the continuous
rate inside the right-hand side, `breakupEvent` the events in `solout`. A KH-RT shed fills the parent's `shedList`
(`gr%shed(:)`) with `childState = [v_p, d_stable, n_child]`; `obj_IGLOO%solve` drains the lists into new
`obj_particle`s and integrates them in the next pass of its outer `maxLoop`. Every parcel's random stream is seeded
from `(seed, family, ID)` (`rngSeedFor`). The KH-RT validity count prints
`[WARNING] Reitz-KHRT outside its hypotheses (Oh > 1 or rho_g/rho_p > 0.1) in N model evaluations`.

---

## The integration loop

`Lib_Integration.f90::integrate` advances one parcel from injection to exit, called from `obj_IGLOO%solve` inside an
`!$OMP PARALLEL DO` over the parcels. The ODE library is OSlo (`oslo::Run_ODESolver`).

**Outer loop** (one iteration per segment, at most `maxIter = 500000`):

1. `initializeCell` / `computeSource`: the entry state.
2. `deltat` from `computeDeltat(vert)`: the distance `computeDs` returns along the velocity over the speed (on the gas
   dual cell at `gas-order = 2`); a non-finite value falls back to `tauFactor = 100` Stokes times.
3. The horizon `t2 = tStart + 20 deltat`.
4. The inner loop: `ODEsystem` up to `max-inner-iter` times, until `exitLoop`.
5. `computeSrcField`, `computeEulField`: the deposits (each an `!$OMP ATOMIC UPDATE`).
6. `updateCell`, the boundary (`bcDef`), the wedge fold.
7. Trajectory and scatter output.

**`ODEsystem`** wraps one `Run_ODESolver` call:

- `safety = 0.99` scales `deltat` before each call.
- `err < 0` (a solver failure): the parcel returns to its last accepted step, is deposited up to it and is `gone`.
- `exitLoop` when the solver reaches XEND, when the step is below `dtMin = 1e-14` or `dout < eps`, or on a NaN.
- On an interrupt (`IamOut`, `newGas`, `sectorOut`, `eventFlag`):
  `deltat = max(deltat * min(1/nStep, din/dout), dtMin)`, `nStep = 10`, `din` the Möller–Trumbore distance to the
  face (exact for a sector plane), floored at `dtSafe`, the time to move a fraction of the containment tolerance.
- A solidification event is applied here, by absolute assignment, never inside `solout` (whose `Y` is a copy under
  `H-sdirk4`).

**`solout`**, at every accepted step:

1. the scatter accumulator (`npRef*h`, markers at multiples of `dNscat`, weight `dNscat*npdot/npRef`);
2. `IamOut = .not. isPointInsideCell(y(1:3), vert, ...)`, `sectorOut` on a wedge, and `newGas` for the gas dual cell
   at `gas-order = 2` (`atGasBoundary` limits the geometry test to boundary dual cells);
3. the burnout test on the mass slot;
4. the NaN and tiny-step guard (`any(y /= y)`, `deltat < dtMin`);
5. `breakupEvent`, and `addChild` for a shed;
6. the solidification event functions;
7. the snapshot `oldLocal`, `oldStLocal`, `oldEvLocal` on a clean step.

Any interrupt returns `IRTRN = -2724`.

**Guards.** `nMaxCell = 10` segments in one cell; every `nMaxStall = 100` reflections the position is compared with
the last checkpoint, closer than `nMaxStall * stallStep` = 1e-4 m is a stall; `max-reflections`; `max-inner-iter`;
`maxIter`. The right-hand sides of models 2–8 replace a non-finite entry by `1e30`. A `deltat` outside
`[dtMin, huge)` is a hard error.

---

## Geometry

Module `IGLOO_RayFaceIntersection3D` (`geometry.f90`):

- `guide(6,4)`: the 4 nodes of each of the 6 faces,
  `reshape([1,2,3,4, 7,6,5,8, 6,2,1,5, 4,3,7,8, 5,1,4,8, 3,2,6,7], [6,4], order=[2,1])`.
- `isPointInsideCell`: the dispatcher. 2D: `isPointInsideQuadrilateral` (faces through `faceMap2D`). 3D with the cache:
  `isInsideHex12Fast`; without: `isPointInsideHexahedron12` (`chooseVector` picks a non-degenerate start node; a
  degenerate triangle returns `distance = -1`).
- `precomputeHex12`: at every cell entry and after a boundary event, caches `geoHexNorms(3,2,6)`,
  `geoHexCentroids(3,2,6)`, `geoHexDegen(2,6)` and the tolerance `geoHexTol` (`hexTol`); `gasHex*` for the dual cell.
- `obj_particles.f90`: `updateCell` (from `exitFace`), `findParticle` (after a boundary event), `advanceGasCell`
  (`gasExitFace`), `computeDs`, `computeDeltat`.
- `obj_bc.f90`: `bcDef` (the codes), `checkBoundary` and `closestFace` (the face crossed), `hopCrossing` (the impact),
  `axisymFold`, `periodicTransport`.

**Trajectory rows at a fold or a wrap.** With `print-dcell` set, a sector crossing writes no row and a periodic wrap is
a cell crossing, whose first row opens the new zone. Without it, the fold and the wrap each write two rows at one time,
the state before and after; the wrap's second row opens the new zone.

---

## Deposits

- `obj_sourceblock` (`obj_block.f90`): `sourceMass(nsm(), Nx, Ny, Nz)`, `sourceMom(3, ...)`, `sourceEn(...)`, and the
  coupling arrays `offered`, `absorbed`, `held`, `heldLast`. `variables.f90`: `nsm()` counts the rows, `rowCarry(m, s)`
  and `rowProc(m, s)` number the optional ones. `source.tec` writes the `nm` net rows only.
- `finalize` reduces the dual mesh to the cells with the weights of `precomputeDualWeights`;
  `assembleInterfaceDual` builds the vertex classes of 101 and 201 faces at setup, `combineInterfaceDual` sums them
  after the MPI reduction.
- `obj_eulerblock`: `density`, `np`, `velocity(3)`, `temperature`; `reduce`, then `normalize` (the Favre division;
  `cpVariable` inverts the enthalpy). `finalizeEUL` divides by the density.

## Mollification

`Lib_Mollify.f90`: `binomial_smooth` (one block), `mollifyBlocks` (blocks joined by 101, 201 and chimera faces),
`theta = 0.25`, `DEFAULT_MOLLIFY_PASSES = 8`. `obj_block.f90` calls them through `mollifyEulerBlocks` and
`mollifySourceBlocks`; without connected faces each block goes through `mollifyEuler` / `mollifySource` with the
scratch buffer `work`, reused over the 6 (Euler) or $n_m + 4$ (source) calls.

## Injection

`initialization.f90::pin_particles`, once per group in `obj_IGLOO%setup`. Method `FB` (faces) or `DB` (given points,
selected by `x`/`y`/`z`): `DB` copies `pos0`, `vel`, `temp`, `mdot`, `diam` into `group%particle`, the cell indices
`[0,0,0,0]` until `updateCell` finds them. `FB` uses `pin_particles_bc_center` when `ds = 0`, `mdotMax = 0` and no cell
carries its own spacing (`dsSwitch`), else `pin_particles_bc_ds`: `compute_ds` gives the spacing, `injViable` applies
`dsDegen`, the 2D sweep runs when `block(1)%Nz == 1`, the 3D front otherwise (the queue `Inj(:,:)`, `ind(:,:)`, `head`;
at most `alloc = 10^7` points; a uniform bucket grid finds the neighbours). `obj_IGLOO.f90::dropZeroFlow` removes the
parcels of zero-loading cells and prints
`[IGLOO] family <f>: <n> parcel(s) of inlet cells with a zero loading (field 1 = 0) not injected`.
`Lib_Statistics.f90::sampleDiameter` draws the diameters (`bc.txt` columns 6–8, `(fam, 6..8)`);
`obj_particles.f90::initializePart` reads the cell's columns.
