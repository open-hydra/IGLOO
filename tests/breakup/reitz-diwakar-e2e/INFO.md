# INFO — breakup/reitz-diwakar-e2e (B-VAL-4: Reitz-Diwakar validation, e2e, **GREEN**)

## Reference papers
- **[RD86]** Reitz, R. D.; Diwakar, R. *"Effect of Drop Breakup on Fuel Sprays."* SAE
  Technical Paper 860469, 1986 (SAE Trans. 95, pp. 218–227). The RD bag/stripping model.
- **[RD87]** Reitz, R. D.; Diwakar, R. *"Structure of High-Pressure Fuel Sprays."* SAE
  Technical Paper 870598, 1987 (SAE Trans. 96, pp. 492–509;
  `papers/Reitz-StructureHighPressureFuel-1987.pdf`). The calibrated form used here.
- **Reproduced result (Tier P):** the RD breakup rate. RD is model 3 (ODE-path, no
  event); the shed rate reduces (mass conservation, `d ∝ npdot^{-1/3}`) to
  `dd/dt = (dStable − d)/τ`, with the branch chosen per drop state.

## The model, verified against [RD87]
| | criterion | timescale τ | dStable |
|---|---|---|---|
| **Bag** | `We_r > 6` and `We_r ≤ 0.5√Re` (Eqs 5,7) | `π·√(ρ_l·r³/2σ)` (constant **D=π**) | `12σ/(ρ_g u²)` |
| **Stripping** | `We_r > 0.5√Re` (Eqs 6,8) | `C·(r/u)·√(ρ_l/ρ_g)`, **C=20** (curve-fit) | `σ²/(ρ_g u³ μ_g)` |

`We_r = ρ_g u² r/σ` (radius); `Re = ρ_g u d/μ_g` (gas, diameter). All four constants
(`WeBag=6`, `Cb=π`, `Cstrip=0.5`, `Cs=20`) and both stable sizes were checked against the
`[RD87]` PDF: Eq. 5 (`We>6`), Eq. 6 (`We/√Re>0.5`), Eq. 7 (bag, `D=π`), Eq. 8 (stripping),
with **C=20 obtained "by curve-fitting"** (p. 497 / Figs 10–11). `[RD86]`'s earlier
statement that the stripping constant `D₂` is "of order unity" was **refined to the
curve-fit C=20 in [RD87]**, which production uses — so the historical `Cs=20` vs `~O(1)`
flag is resolved: **no bug** (production is `[RD87]`-faithful).

## Case construction
- `tools/make_pe_case.py --we-convention rad --u-gas 200`, 25-point radius-based sweep
  `We_r ∈ {8 … 1000}`, slip = 100 m/s (`κ_v=0.5`), σ=0.072, μ_l=1e-3, water (ATLAS-GPB
  fixtures). At slip=100 the bag→stripping handoff (`We_r = 0.5√Re`) sits at **We_r=20**
  (`Re = 80·We_r`), so the sweep spans **bag** (We_r 8–18, 6 drops) and **stripping**
  (We_r 22–1000, 19 drops).
- **Weber convention:** radius-based (RD and production `Lib_Breakup.f90:154`).

## What `check.py` gates
Each drop's integrated **initial breakup rate** (leading-window LSQ, matched against the
same-window RK4 of the RD rate — the decelerating rate's finite-window bias cancels on
both sides) vs the RD closed form for that drop's regime, at the injection We. Tol 2 %
(observed worst 8.3e-4); requires ≥16 gated drops with ≥3 in each regime. RD is model 3
non-event, so there is no event-path scoping (unlike KHRT B-VAL-6): all breaking drops
are gated. Result: 25/25 within 0.1 % (bag 6, stripping 19), deterministic
(3× OMP=5 + OMP=1).

Lockstep note (as for KHRT): `check.py` codes the same RD correlation as production, so
this validates the **integration path** (cell crossings, the model-3 npdot→d reduction,
the bag/stripping branch selection) end-to-end; the pointwise RD formula (both branches +
the handoff at `We_r=0.5√Re`) is independently pinned by the unit test `test_breakup_rd`.

## Comparison plot (`verify.py` → `OUTPUT/reitz-diwakar-e2e.svg`, non-gating)
Initial `dd/dt` vs `We_r`: the RD closed-form rate curve and IGLOO's measured rates
coloured by regime (bag ○ / stripping □), with the bag→stripping handoff line.

## Tier-2 (spray-level) documented rejection
RD's headline validation in `[RD87]` is diesel-spray penetration / tip structure — out of
scope by design: gas entrainment + two-way coupling + an atomizing nozzle inlet, none of
which IGLOO's steady one-way carrier gas provides. Revisit only if a transient/two-way
mode is added.
