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

## What `check.py` gates (KH stripping, GREEN)
IGLOO's initial `dd/dt` (leading-window LSQ, matched against the same-window RK4 of the
Reitz-87 KH rate — the decelerating rate's finite-window bias cancels on both sides) vs
the analytic rate at the injection We. **Scope: `We_r ≥ 340`, 14 drops, tol 2 %**
(observed worst 2e-4). Result: 14/14 within 0.03 %.

## Why the KH gate is scoped to the initial rate at We_r≥340
The `check.py` KH-stripping gate uses the **initial rate only** (first ~3 % loss) at
`We_r≥340`, where the initial window is RT-free (below that the drop decelerates hard
enough that RT fires within the window). This keeps the KH-rate gate a clean measurement
of the continuous Kelvin-Helmholtz ODE, separate from the discrete RT/shed events handled
by `check_rt.py`. (Historically this scoping also made the KH gate A19-independent, before
A19 was fixed 2026-07-23 — see below.)

Lockstep note: `check.py` codes the same `λ_KH/Ω_KH` correlation as production
(paper-anchored via A2/A13 + the pointwise `test_breakup_khrt` unit test), so this
validates the **integration path** (cell crossings, the model-3 npdot→d reduction)
end-to-end — not the Reitz-87 correlation itself.

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
