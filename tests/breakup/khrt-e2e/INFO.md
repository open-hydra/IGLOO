# INFO — breakup/khrt-e2e (B-VAL-6: KHRT KH-stripping validation, e2e, **GREEN, RT gate GREEN (A19 fixed)**)

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
  plan's Tier-V primary (`r(t)=r_s+(r0−r_s)e^{−t/τ_KH}`); gated by the same initial-rate
  method as PE (B-VAL-1).

## Case construction
- `tools/make_pe_case.py --we-convention rad --u-gas 200` + 25-point radius-based sweep
  `We_r ∈ {30 … 1000}` (d 0.36–12 mm), slip = 100 m/s (`κ_v=0.5`), σ=0.072, μ_l=1e-3,
  water drops (ATLAS-GPB fixtures). Validity: Oh≈0.005<1 and ρg/ρp≈0.0012<0.1 (the KHRT
  hypothesis holds; no `[WARNING] ReitzKHRT hypothesis falling`).
- **Weber convention:** radius-based (KHRT/Reitz and production `Lib_Breakup.f90:182`).

## KH-shed children (A23g/A23h **CLOSED** 2026-08-04, one-shed cap lifted)
This case emits **241 KH-shed child parcels** in a single pass from 25 parents (deepest parcel
sheds 22 times; `run_out.txt` prints `number of children = 241` plus a `max N sheds from one
parcel` tripwire), so `outloc-A.dat` carries 266 exit records. The children integrate: they are
born inside the domain at the shed point, fly, and RT-shatter in their own right (`check_rt.py`
counts 134 shatters, up from 8 when only parents flew).

Three defects were closed to get here, all previously invisible to a green suite:

- **A23g** — `vert`/`gasVert`/`gas` were filled only inside the `part%time==0` init block, yet
  `computeDeltat(vert)` runs on every outer-loop entry. Any parcel entering with `time /= 0` —
  i.e. every child — sized its first step off uninitialized stack and integrated that segment
  against uninitialized gas. Now refreshed unconditionally via `refreshCellEntry()`.
- **A23h** — the child origin capture sat in the outer loop guarded only by the *latched*
  `addChild`, so it re-ran every iteration and walked the record forward to the parent's final
  state. Children were born on the exit plane with **zero** trajectory records. Capture is now
  single-point, in the finalize block after `updatePart`.
- **the one-shed cap** — `gr%child(ip)` was one slot per parent, so a parent could shed exactly
  once by construction. Replaced by a growable per-parent `shedList`.

⚠ **"all children exit at `x = 0.150000`" is NOT a signature of A23h**, though it was recorded
as one while the bug was open. Every parcel leaves through that same outflow plane — parents
included — so the equality holds just as firmly for correct output. What actually separated the
broken state from the fixed one was whether a child left a **trajectory** behind: 0 records for
all 25, versus 58-59 records each now. `check.py::check_children_integrate` gates that instead,
and was confirmed to go red against the reinstated bug before being landed.

## Mass conservation (A23b+A23d, gated 2026-08-03)
`check.py::check_mass_conservation` asserts total exit parcel mass-flow == injected
(3.09091e-01, tol 1e-4). This gate exists because KHRT **created ~29 % mass** undetected for the
whole validation campaign: the A19 finalize re-derives `mdot = npdot·ρ·d³·π/6`, and both inputs
were stale in different ways (`d` from a once-per-call `auxLocal`; `npdot` from segment start).
Pre-fix exit 3.99106e-01 vs injected 3.09091e-01; post-fix 3.09090e-01 (3.3e-06). **Do not remove
this gate without replacing it** — nothing else in the suite looks at mass.

## Shed-path gates (A23/S4, added 2026-08-04)
Three gates cover the KH-shed path, none of which the older ones could see:

- **`check.py::check_children_integrate`** (A23h). A *population* test: ≥90 % of the shed
  children must leave ≥4 trajectory records and cover ≥1e-3 in x between birth and exit.
  Observed 98.3 %; the broken state scores 0 %. It is a population test on purpose — a parent
  may legitimately shed in its final segment, giving a child that exits before the periodic
  trajectory print fires (measured: 1 such child in 241). Requiring 100 % would be flaky.
  Deliberately **pairing-free**: `kid%ID` follows the drain order over shed lists, not the
  parent index, so the apparent `+25` parent↔child mapping is not a real invariant.
- **`check.py::check_multiple_sheds`** (the cap lift). Pass-1 children (241) must exceed the
  25 injected parents, which is only possible if some parent shed more than once. Without it,
  a silent return of the one-shed cap would leave every other gate green.
- **`khrt-e2e-threads`** → `check_threads.py`. Runs the case at `OMP_NUM_THREADS` 1/2/4 and
  compares `trajectories-A.dat`/`outloc-A.dat` as **sorted multisets**, so record ordering
  (OMP-nondeterministic by design) is ignored while a *torn* record is not — those two files
  are written from inside the parallel region with no synchronisation, and nothing else in the
  suite would catch interleaving. Shares `RESOURCE_LOCK khrt_output` with the other two.

Heavy shedding is stressed separately by [breakup/khrt-stress](../khrt-stress/INFO.md).

## What `check.py` gates (KH stripping, GREEN)
IGLOO's initial `dd/dt` (leading-window LSQ, matched against the same-window RK4 of the
Reitz-87 KH rate — the decelerating rate's finite-window bias cancels on both sides) vs
the analytic rate at the injection We. **Scope: `We_r ≥ 340`, 14 drops, tol 2 %**
(observed worst 2e-4). Result: 14/14 within 0.03 %. Unaffected by the shed work — it selects
parents by `We_r`, and the children are far below the cut.

## Why the KH gate is scoped to the initial rate at We_r≥340
The `check.py` KH-stripping gate uses the **initial rate only** (first ~3 % loss) at
`We_r≥340`, where the initial window is RT-free (below that the drop decelerates hard
enough that RT fires within the window). This keeps the KH-rate gate a clean measurement
of the continuous Kelvin-Helmholtz ODE, separate from the discrete RT/shed events handled
by `check_rt.py`. (Historically this scoping also made the KH gate A19-independent, before
A19 was fixed 2026-07-23 — see below.)

**RT-freeness re-measured 2026-08-03 (O7).** The `We_r≥340` cut was originally justified
by reasoning that assumed A19 was still open. With A19 fixed, that justification no longer
holds a priori, so the window was measured directly against the fixed code: across all
**14** gated drops the largest single-interval diameter ratio inside the leading window is
**1.001** — an RT shatter is ≥1.8× — so **zero** RT events intrude. The threshold stands on
evidence now, not on the bug.

Lockstep note: `check.py` codes the same `λ_KH/Ω_KH` correlation as production
(paper-anchored via A2/A13 + the pointwise `test_breakup_khrt` unit test), so this
validates the **integration path** (cell crossings, the model-3 npdot→d reduction)
end-to-end — not the Reitz-87 correlation itself.

**O7 — that independence gap is closed unit-side (2026-08-03).** Two new unit tests in
`tests/breakup/reitz-khrt/`:
- **`test_kh_rayleigh_limit`** — the only KH reference that shares no constant with the
  fit. It recovers `λ_KH`/`Ω_KH` from `breakupOde` itself (the pure-KH rate is *affine in
  `B0`*, so two calls invert it exactly — no production change, no third copy of the
  correlation) and checks the `We_g→0, Z→0` corner against **Rayleigh's 1878** capillary
  dispersion relation `ω² = (σ/ρa³)·x(1−x²)I₁(x)/I₀(x)`, solved in-test: `λ/a → 9.0144`,
  `ω*→0.34334` vs the fit's `9.02`/`0.34` (6.2e-4 and 9.7e-3 — the fit's own accuracy,
  which is the tolerance floor). Reitz-87 p.318 states this reduction explicitly ("the
  maximum growth rate occurs at Λ = 9.02a").
- **`test_khrt_interaction`** — the KH↔RT path agreement, O7 as filed in BUGS.md.

### Provenance of the `0.87` in Eq. (4)
Reitz-87 Eq. (4) prints `(1 + 0.87·We_g^1.67)^0.6`; **Beale-Reitz 1999** (`KHRT-hybrid.pdf`)
restates the same fit with **`0.865`**. Production follows the **origin paper (0.87)**. The
difference is ≈0.35 % on `Λ` at high `We_g` and vanishes as `We_g→0`, so the Rayleigh anchor
above does *not* discriminate it — recorded here so the choice is deliberate and traceable
rather than an accident of transcription.

## The RT gate (`check_rt.py`, ctest `khrt-e2e-rt`) — real gate since A19 fixed
Asserts the RT/shed **event** fires, applies, and PERSISTS: ≥6 drops (8 do, deterministically)
show a DISCONTINUOUS single-step `d`-collapse (≥1.8× in one recorded interval — the
continuous KH ODE can never do that) that reaches the fragment scale (`d≤0.5·d0`) early
and stays there, with per-droplet mass self-consistent (`m=ρd³·π/6`). This gates the
INTEGRATION property that A19 broke (the event surviving into the ODE state); the `λ_RT`
child diameter and the `(dp/λ_RT)³` cubic child count (Beale-Reitz Eq. 11, A14) are
verified pointwise by the unit test `test_breakup_khrt` — a finite-difference `λ_RT`
reconstruction e2e is too noisy for a tight gate (`d_after/λ_RT` scatters ~1–2.4).

## Bug A19 — FIXED 2026-07-23
The RT/shed event's npdot change was not propagated into the model-3 ODE state `y(8)`, so
`updatePart` reverted it from the stale `stateVar(8)` (trace: `part%npdot=4.23e7` (event)
vs `stateVar(8)=6.76e5` (stale)). Fixed on the event finalize path with three coupled
parts: push the event npdot into `y(8)` and re-derive `mdot` consistent with the event
`(d, npdot)` — mass-conserving for both RT split and KH shed — guarded on `eventFlag` so a
no-event finalize can't freeze the continuous stripping; plus the `oldStLocal` segment-start
init (the A18 defect-4 twin, for the RT `told/tc` timer). Byte-inert on all 16 non-event
e2e cases incl. the ord2 evap cases.

## Tier-2 (spray-level) documented rejection
Diesel-spray penetration / dense-spray SMD (Reitz-KHRT's headline spray validation) is out
of scope by design: gas entrainment + two-way coupling + atomizing nozzle inlet, none in
IGLOO's steady one-way carrier.
