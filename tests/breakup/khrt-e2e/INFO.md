# INFO — breakup/khrt-e2e (KHRT KH-stripping validation, e2e, **GREEN, RT gate GREEN**)

## Reference paper
- **Title:** Modeling atomization processes in high-pressure vaporizing sprays (KHRT /
  Kelvin-Helmholtz–Rayleigh-Taylor hybrid)
- **Author(s):** R. D. Reitz (1987); Beale & Reitz (1999) hybrid
- **Tag:** `[Reitz87]` (`papers/Reitz-Atomization1987.pdf`)
- **Reproduced result (Tier V/P):** the **continuous KH-stripping rate**. KHRT (model 3)
  integrates a shed rate that reduces (mass-conservation, `d ∝ npdot^{-1/3}`) to
  ```
  dd/dt = (dStable − d)/τ_KH,   dStable = 2·B0·λ_KH,   τ_KH = 3.726·B1·r/(λ_KH·Ω_KH)
  ```
  with the Reitz-87 KH growth rate `Ω_KH` and wavelength `λ_KH` correlations. This is the
  Tier-V primary (`r(t)=r_s+(r0−r_s)e^{−t/τ_KH}`); gated by the same initial-rate
  method as [pilch-erdman-e2e](../pilch-erdman-e2e/INFO.md).

## Case construction
- `tools/make_pe_case.py --we-convention rad --u-gas 200` + 25-point radius-based sweep
  `We_r ∈ {30 … 1000}` (d 0.36–12 mm), slip = 100 m/s (`κ_v=0.5`), σ=0.072, μ_l=1e-3,
  water drops (ATLAS-GPB fixtures). Validity: Oh≈0.005<1 and ρg/ρp≈0.0012<0.1 (the KHRT
  hypotheses hold; no `[WARNING] Reitz-KHRT outside its hypotheses` line).
- **Weber convention:** radius-based (KHRT/Reitz and production `Lib_Breakup.f90` `ReitzKHRT`).

## KH-shed children
This case emits **241 KH-shed child parcels** in a single pass from 25 parents (deepest parcel
sheds 22 times; `run_out.txt` prints `number of children = 241` plus a `max N sheds from one
parcel` tripwire), so `outloc-A.dat` carries 266 exit records. The children integrate: they are
born inside the domain at the shed point, fly, and RT-shatter in their own right (`check_rt.py`
counts 187 shatters; 8 on a solver whose shed children do not fly).

Three pieces of the shed path make the children fly; without the shed-path gates below and
[khrt-stress](../khrt-stress/INFO.md), the suite is green on a solver missing any of them:

- **Cell-entry state.** `computeDeltat(vert)` runs on every outer-loop entry, and a child enters
  with `time /= 0`; `refreshCellEntry()` fills `vert`/`gasVert`/`gas` on every entry. A solver that
  fills them only inside the `part%time==0` init block sizes every child's first step off
  uninitialized stack and integrates that segment against uninitialized gas.
- **Child origin.** The child record is captured once, at the shed point, in the finalize block
  after `updatePart`, and pushed onto the parent's shed list. A capture in the outer loop guarded
  only by the *latched* `addChild` re-runs every iteration and walks the record forward to the
  parent's final state: the children are born on the exit plane with **zero** trajectory records.
- **Shed count.** Each parent keeps a growable `shedList`; with one slot per parent
  (`gr%child(ip)`) a parent sheds exactly once by construction.

⚠ **"all children exit at `x = 0.150000`" is NOT a signature of children born on the exit
plane.** Every parcel leaves through that same outflow plane — parents included — so the equality
holds just as firmly for correct output. What separates the two is whether a child leaves a
**trajectory** behind: on a solver without the three pieces, 0 records for all 25 children (one
per parent); here 2 to 61 records each (median 41). `check.py::check_children_integrate` gates that instead,
and goes red with the outer-loop capture reinstated (measured).

## Mass conservation
`check.py::check_mass_conservation` asserts total exit parcel mass-flow == injected
(3.09091e-01, tol 1e-4); measured 3.09091e-01 (6.4e-07). The event finalize re-derives
`mdot = npdot·ρ·d³·π/6`; a solver that takes both inputs stale (`d` from a once-per-call
`auxLocal`, `npdot` from segment start) **creates ~29 % mass** (exit 3.99106e-01 against the
injected 3.09091e-01), which the KH and RT gates do not see. **Do not remove this gate without
replacing it** — in this case no other check compares the exit mass flow with the injected flow.

## Shed-path gates
Three gates cover the KH-shed path, which the case's other gates cannot see:

- **`check.py::check_children_integrate`** (the child origin). A *population* test: ≥90 % of the
  shed children must leave ≥5 trajectory records (the birth row and four more) and cover ≥1e-3 in
  x between birth and exit. Observed 98.8 %; a solver whose children are born on the exit plane
  scores 0 %. It is a population test on purpose — a parent may legitimately shed in its last
  cells, giving a child that exits before the periodic trajectory print has fired four times
  (measured: 3 such children in 241, born in the last three cells). Requiring 100 % would be flaky.
  Deliberately **pairing-free**: `kid%ID` follows the drain order over shed lists, not the
  parent index, so the apparent `+25` parent↔child mapping is not a real invariant.
- **`check.py::check_multiple_sheds`** (the shed count). Pass-1 children (241) must exceed the
  25 injected parents, which is only possible if some parent sheds more than once. Without it,
  a one-shed cap would leave every other gate green.
- **`khrt-e2e-threads`** → `check_threads.py`. Runs the case at `OMP_NUM_THREADS` 1/2/4 and
  compares `trajectories-A.dat`/`outloc-A.dat` as **sorted multisets**, so record ordering
  (OMP-nondeterministic by design) is ignored while a *torn* record is not — those two files
  are written from inside the parallel region with no synchronisation; the `threads-*` gates
  catch interleaving on their own cases, this one on the KH-shed path. Shares
  `RESOURCE_LOCK khrt_output` with the other two.

Heavy shedding is stressed separately by [breakup/khrt-stress](../khrt-stress/INFO.md).

## Source budget (`check.py::check_source_budget`)
The case writes `source.tec`. A segment deposits `in − out` of the parcel's mass, momentum and energy flow,
so a parcel's deposits telescope to (first in − last out). At a KH shed the parent's outgoing flow drops by
the child's birth flow and the child starts from it: the parent hands that flow on to the child, it does not
deposit it. Over the whole field, with the flow along x,

    Σ Fx = P_in − Σ_exits ṁ |u_p|,    Σ E = E_in − Σ_exits ṁ (c_p T + |u_p|²/2),    Σ wdot = Σ Fy = Σ Fz = 0

injected state from the parents' injection rows and `MDOT_INJECTED`, exit state from outloc (columns 7, 5, 4).
Tolerances: `5e-6` of `P_in`/`E_in`, the E13.6 bound on the printed exit flows (measured **6.1e-7** and
**6.4e-7**); transverse and mass `1e-12` (exactly 0 measured). The cell-local form of the same statement,
with the exchange switched off, is [breakup/khrt-shed-noexchange](../khrt-shed-noexchange/INFO.md).

| solver | Σ Fx vs P_in − P_out | Σ E vs E_in − E_out |
|---|---|---|
| shed flow deposited in the gas | −5.265758 vs −13.017249 N: **+7.751 N** (0.25 of P_in) | +9.52172e4 vs −1.88588e3 W: **+9.710e4 W** (0.25 of E_in) |
| shed flow handed to the child | −13.017230 vs −13.017249 N (6.1e-7) | −1885.627 vs −1885.876 W (6.4e-7) |

The excess is the children's birth flux exactly: divided by their exit flow (0.0770868 kg/s) it is
100.555 m/s and c_p·300 + ½·100.555² J/kg.

## What `check.py` gates (KH stripping, GREEN)
IGLOO's initial `dd/dt` (leading-window LSQ, matched against the same-window RK4 of the
Reitz-87 KH rate — the decelerating rate's finite-window bias cancels on both sides) vs
the analytic rate at the injection We. **Scope: `We_r ≥ 340`, 14 drops, tol 2 %**
(observed worst 2e-4). Result: 14/14 within 0.03 %. The shed children do not enter it — it
selects parents by `We_r`, and the children are far below the cut.

## Why the KH gate is scoped to the initial rate at We_r≥340
The `check.py` KH-stripping gate uses the **initial rate only** (first ~3 % loss) at
`We_r≥340`, where the initial window is RT-free (below that the drop decelerates hard
enough that RT fires within the window). This keeps the KH-rate gate a clean measurement
of the continuous Kelvin-Helmholtz ODE, separate from the discrete RT/shed events handled
by `check_rt.py`.

**RT-freeness, measured.** Across all **14** gated drops the largest single-interval diameter
ratio inside the leading window is **1.001** — an RT shatter is ≥1.8× — so **zero** RT events
intrude: the `We_r≥340` threshold stands on that measurement.

Lockstep note: `check.py` codes the same `λ_KH/Ω_KH` correlation as production
(paper-anchored by the pointwise `test_breakup_khrt` unit test, whose oracle is written from
the Reitz-87 paper form, Eqs 4–5: the Ohnesorge `Z` on the liquid properties and the radius, the
Taylor `T = Z·√We_g` on the gas Weber), so this
validates the **integration path** (cell crossings, the model-3 npdot→d reduction)
end-to-end — not the Reitz-87 correlation itself.

**The independence gap, unit side.** Two unit tests in `tests/breakup/reitz-khrt/` cover it:
- **`test_kh_rayleigh_limit`** — the only KH reference that shares no constant with the
  fit. It recovers `λ_KH`/`Ω_KH` from `breakupOde` itself (the pure-KH rate is *affine in
  `B0`*, so two calls invert it exactly — no production change, no third copy of the
  correlation) and checks the `We_g→0, Z→0` corner against **Rayleigh's 1878** capillary
  dispersion relation `ω² = (σ/ρa³)·x(1−x²)I₁(x)/I₀(x)`, solved in-test: `λ/a → 9.0144`,
  `ω*→0.34334` vs the fit's `9.02`/`0.34` (6.2e-4 and 9.7e-3 — the fit's own accuracy,
  which is the tolerance floor). Reitz-87 p.318 states this reduction explicitly ("the
  maximum growth rate occurs at Λ = 9.02a").
- **`test_khrt_interaction`** — the KH↔RT path agreement: `ReitzKHRT` (ODE) and
  `ReitzKHRTevent` (event) each carry their own copy of the KH and RT blocks, and the test pins,
  at identical state, that the ODE's KH rate is exactly zero iff the event path takes the RT
  branch.

### Provenance of the `0.87` in Eq. (4)
Reitz-87 Eq. (4) prints `(1 + 0.87·We_g^1.67)^0.6`; **Beale-Reitz 1999** (`KHRT-hybrid.pdf`)
restates the same fit with **`0.865`**. Production follows the **origin paper (0.87)**. The
difference is ≈0.35 % on `Λ` at high `We_g` and vanishes as `We_g→0`, so the Rayleigh anchor
above does *not* discriminate it — recorded here so the choice is deliberate and traceable
rather than an accident of transcription.

## The RT gate (`check_rt.py`, ctest `khrt-e2e-rt`)
Asserts the RT/shed **event** fires, applies, and PERSISTS: ≥6 drops (187 here, parents and
children; 8, deterministically, on a solver whose shed children do not fly)
show a DISCONTINUOUS single-step `d`-collapse (≥1.8× in one recorded interval — the
continuous KH ODE can never do that) that reaches the fragment scale (`d≤0.5·d0`) early
and stays there, with per-droplet mass self-consistent (`m=ρd³·π/6`). This gates an
INTEGRATION property, the event surviving into the ODE state (the write-back below). The
`λ_RT` child diameter is pinned pointwise by `test_khrt_interaction` (the RT branch lands the
parcel on `λ_RT`, 1e-12); the `(dp/λ_RT)³` cubic child count (Beale-Reitz Eq. 11) has no
pointwise test (`test_breakup_khrt` covers the pure-KH path). A finite-difference `λ_RT`
reconstruction e2e is too noisy for a tight gate (`d_after/λ_RT` scatters ~1–2.4).

## The event write-back (what the RT gate sees)
On the event finalize path the solver pushes the event's npdot into the model-3 ODE state
`y(8)` and re-derives `mdot` from the event `(d, npdot)` — mass-conserving for both the RT
split and the KH shed — guarded on `eventFlag`, so a finalize without an event cannot freeze
the continuous stripping. The RT `told/tc` timer lives in the aux state, whose segment-start
snapshot `oldStLocal` (taken as `oldEvLocal` is, see [tab-e2e](../tab-e2e/INFO.md)) restores it
when a crossing cuts the step. A solver
without the write-back lets `updatePart` revert the event from the stale `stateVar(8)`
(trace: `part%npdot=4.23e7` (event) vs `stateVar(8)=6.76e5` (stale)): the event does not
survive into the ODE state, and the RT gate fails.

## Spray-level rejection
Diesel-spray penetration / dense-spray SMD (Reitz-KHRT's headline spray validation) is out
of scope by design: gas entrainment + two-way coupling + atomizing nozzle inlet, none in
IGLOO's steady one-way carrier.
