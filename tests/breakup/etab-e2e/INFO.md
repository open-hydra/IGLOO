# INFO — breakup/etab-e2e (ETAB cascade validation, e2e, **GREEN**)

## Reference paper
- **Title:** Liquid jet atomization and droplet breakup modeling of non-evaporating
  diesel fuel sprays (ETAB / Enhanced TAB cascade)
- **Author(s):** F. X. Tanner
- **Year:** 1997 (SAE 970050) / 1998 — **Tag:** `[Tan97/98]`
  (`papers/ETAB-Tanner1997.pdf`, `ETAB-Tanner1998.pdf`)
- **Reproduced result:** ETAB reuses the TAB damped oscillator (so onset `We_r=6` and
  the first-breakup time `t_bu` are the ORA87 closed form of [tab-e2e](../tab-e2e/INFO.md)), and replaces
  TAB's stochastic Rosin-Rammler product sizing with a **deterministic exponential
  cascade** (Tanner eq. 6/8):
  ```
  d_child/d_parent = exp( -(Kbr/ω)·acos(1 - 1/WeCr) ),  WeCr = We/WeCrit(=12)
  Kbr/ω = k1·(AWe·We⁴ + 1)   (bag,   We ≤ WeTrans=80)
        = k2·√We             (strip, We > WeTrans)
  AWe   = (k2/k1·√WeTrans - 1)/WeTrans⁴
  ```
  `ω` cancels → the ratio is a function of `We` alone. Defaults k1=k2=0.2222,
  WeCrit=6, WeTrans=80 (Tanner-1998 table, the reader's defaults in `Lib_INI.f90`); the
  stripping branch takes k2, continuous with the AWe smoother. The stripping branch asymptotes to
  `exp(-k2·√24) ≈ 0.337` — a size-independent ratio, the ETAB stripping signature.

## Case construction
- `tools/make_pe_case.py --we-convention rad --u-gas 50 --lx 0.6 --nx 240` +
  25-point **radius-based** sweep `We_r ∈ {3 … 105}` spanning the onset (6) and the
  bag→strip transition (80). Slip = 25 m/s (`κ_v=0.5`), σ=0.072, water drops (ATLAS-GPB
  `phase.txt`/`properties.dat`, ρ=1000/cp=4182). Box 0.6 m × 240 cells so `t_bu`
  (1.7–20 ms) is well sampled inside the ~24 ms residence.
- **Weber convention:** radius-based (ETAB/Tanner and production `Lib_Breakup.f90` `ETABmodel`).
- **`gas-order=1` (deliberate):** ETAB is model 1 (no auxState breakup fields), so the
  segment-start aux-state snapshot `oldStLocal` plays no part here.

## What `check.py` gates
- **A (onset):** `We_r ≤ 5.75` → `d` exactly constant (max deformation 2·We_r/12 < 1).
- **B (timing):** `We_r ≥ 6.5` → first `t_bu` inside the measured bracket (ORA87/TAB
  oscillator, slack 5 %) — the event-breakup path of [tab-e2e](../tab-e2e/INFO.md), for ETAB.
- **C (cascade size, Tier V):** first-break ratio `d1/d0` vs the Tanner formula
  evaluated at the **We at breakup** (slip just before the d-drop row — the drop
  decelerates, so injection-We ≠ break-We and the bag ratio is steeply We-dependent;
  the strip plateau is We-flat so it is drift-insensitive). Tol 1 %; observed worst
  **0.03 %** — the ETAB cascade reproduces Tanner to machine-ish precision.
- **D (kick transport):** the v⊥ product kick (Tanner eq. 8-10) reaches the parcel state:
  D1 some broken drop carries lateral velocity, D2 no unbroken drop does, D3
  `0 < max|v_lat|/max|U| < 0.10` (a sanity envelope; measured 0.0215). The rule is exact, with no
  tolerance: with the kick discarded, `V` and `W` are exactly zero in every row. 15 of the 20
  broken drops carry a kick (the low-We end is in the `A² ≤ 0` regime, which resizes without a
  kick), so D does not gate equality. The kick's magnitude and direction (azimuth-random) are
  verified at routine level, ET4a/b in [../etab](../etab/INFO.md), the no-kick branch by ET4c.

## Result (GREEN)
20/20 supercritical drops break; timing sub-percent; child-size ratio within 0.03 % of
Tanner across both branches (bag 0.55→0.79→0.34, strip plateau 0.337). Together with
`tab-e2e` this exercises the event-breakup path (ETAB, like TAB, is model 1).

## Spray-level rejection
Diesel-spray penetration / dense-spray SMD (Tanner's *other* headline validation) is
out of scope by design: gas entrainment + two-way coupling + atomizing nozzle inlet,
none available in IGLOO's steady one-way carrier. Revisit only with a transient/two-way
mode.
