# Breakup — Literature-Grounded Verification Tests

Scope: secondary breakup models in [Lib_Breakup.f90](../../src/lib/Lib_Breakup.f90).
Selector `assign_breakup` ([L28](../../src/lib/Lib_Breakup.f90#L28)):
Pilch-Erdman(1), Reitz-Diwakar(2), Reitz-KHRT(3), TAB(4), ETAB(5).
Defaults read in [Lib_INI.f90 L160-231](../../src/lib/Lib_INI.f90#L160).

Theory-only deliverable (no code run). "Coded ≠ paper" flagged as ⚠ CODE-VS-PAPER — NOT fixed here.

**UPDATE 2026-07-15 — all five models re-verified against the PRIMARY PDFs the user put in
`papers/` (full symbol-by-symbol tables: `plan-bucket/breakup-litcheck-report-2026-07-15.md`;
outcome ledger: A13/A14/A15):** Pilch-Erdman and TAB CLOSED as fully
verified (every previously-inferred constant confirmed, incl. PE87's 0.75/3 velocity poly,
4.5/1.2/1.64 Oh-override, and TAB's K=10/3 child energy balance). KHRT: KH curve fits and RT
structure confirmed; **A13** (Z/T definitions — wrong fluid/length/Weber, ~29× on T) and
**A14** (RT child count missing cube) found & FIXED. ETAB: oscillator/Kbr/size law confirmed;
**A15** (stripping branch k1→k2 + Tanner-1998 defaults 0.2222/We_t=80) FIXED; child
v⊥=A·ẋ spray-cone component ~~NOT-IN-CODE~~ **IMPLEMENTED 2026-07-16** (Tanner Eqs 8–10,
gated ET4a/b/c: closed-form magnitude, orthogonality, A²<0 no-kick). Reitz-Diwakar: criteria, forms,
stable sizes confirmed; time-scale defaults follow the Reitz-1987 lineage — its primary
("Structure of High-Pressure Fuel Sprays", 1987) is still missing (Cs=20 unverifiable).

---

## 0. DOI verification

| Tag | Paper | DOI | Resolves? |
|---|---|---|---|
| `[ORA87]` | O'Rourke & Amsden, TAB Method, SAE 872089 (1987) | 10.4271/872089 | **verified** — SAE page + OSTI biblio 6118786 confirm title/authors |
| `[RD86]` | Reitz & Diwakar, Effect of Drop Breakup on Fuel Sprays, SAE 860469 (1986) | 10.4271/860469 | **verified (title)** — SAE Mobilus lists 860469 w/ matching title; DOI format canonical, PDF not opened |
| `[PE87]` | Pilch & Erdman, IJMF 13(6) 741-757 (1987) | 10.1016/0301-9322(87)90063-2 | **verified** — ScienceDirect PII 0301932287900632 = this DOI, title/pages match |

Cross-check implementations used as secondary authorities (not citable primaries):
ANSYS/CFX theory (TAB Eq. 24, reproduced in Bhandarkar-Manna-Chakraborty, *Int. J. Fluid Mech. Res.* 27(1) 2017); OpenFOAM v7 `ReitzDiwakar.C` (RD reference form). Both confirmed independently below.

---

## 1. TAB — model (4) — STRONGEST ORACLE

### Paper / equation
[ORA87]. Distortion `y` (equator displacement `x` normalized: `y = x/(Cb·r)`, breakup at `y>1`) obeys a linear damped driven oscillator (ANSYS/CFX Eq. 24, = ORA87):

    ÿ = (Cf/Cb)·(ρg/ρl)·(u²/r²)  −  (Ck·σ/(ρl·r³))·y  −  (Cd·μl/(ρl·r²))·ẏ

ORA87 constants: **Cb=1/2, Ck=8, Cd=5, Cf=1/3, K=10/3**.

Grouped forms:
- damping   `1/td = (Cd/2)·μl/(ρl·r²)`
- frequency `ω² = Ck·σ/(ρl·r³) − 1/td²`
- equilibrium (steady, ÿ=ẏ=0) `y_eq = (Cf/(Cb·Ck))·We_r` with radius Weber `We_r = ρg·u²·r/σ`.
  With the constants `Cf/(Cb·Ck) = (1/3)/((1/2)·8) = 1/12` ⇒ **y_eq = We_r/12**.

> ⚠ Task-brief note: the hint `We_eff = Cf/Ck·(…)` omits the `1/Cb`. Correct grouping is `Cf/(Cb·Ck)=1/12`, NOT `Cf/Ck=1/24`. Code reproduces 1/12 (below), so the code is right and the hint was under-specified.

Closed form at constant forcing (heavy/fixed-velocity drop ⇒ `u`,`r` const ⇒ constant coeffs), with `We_eff ≡ y_eq`:

    y(t) = We_eff + e^(−t/td)·[ (y0−We_eff)·cos ωt + ( (y0−We_eff)/td + ẏ0 )·sin(ωt)/ω ]

### Code cross-check — MATCHES
[TABmodel L372-373](../../src/lib/Lib_Breakup.f90#L372):
`rdt = 0.5·bp(2)·mup/(rhop·radius²)` = `(Cd/2)μl/(ρl r²) = 1/td`, with `bp(2)=Cmu=5=Cd` ✓
`omega = bp(1)·sigma/(rhop·radius³) − rdt²` = `Ck σ/(ρl r³) − 1/td²`, with `bp(1)=Comega=8=Ck` ✓ (stored SQUARED).
[L377-378](../../src/lib/Lib_Breakup.f90#L377): `We = radius·vel²/sigma·rhog` = We_r ✓; `WeCr = We/bp(3)`.
[YupdateTAB L554-555](../../src/lib/Lib_Breakup.f90#L554):
`y = WeCr + exp(−dt·rdt)·(y1·cost + (yDot+y1·rdt)·sint/omega)` with `y1=y−WeCr`, `cost=cos(ω dt)`, `sint=sin(ω dt)`.
⇒ term-for-term identical to the closed form with `We_eff=WeCr`. **EXACT MATCH.** `yDot` line is the analytic derivative — also exact.

Product-drop size [L405](../../src/lib/Lib_Breakup.f90#L405): `k = 1 + 4/3 + ρl r³/(8σ)·ẏ²`, `rs=r/k`.
ORA87 energy balance: `r/r_new = 1 + (8K/20) + (ρl r³ ẏ²/σ)·(6K−5)/120`. With **K=10/3**: `8K/20=4/3` ✓ and `(6K−5)/120=1/8` ✓. **K=10/3 confirmed, MATCHES.**

`bp(3)` doubling — CORRECT DESIGN, not a bug. [Lib_INI.f90 L205-206](../../src/lib/Lib_INI.f90#L205): `WeCrit` default 6, then `bp(3)=2·bp(3)` (the `if(err) x=6; x=2*x` semicolon-guard makes the ×2 unconditional — INTENTIONAL). So `bp(3)=12` = the equilibrium denominator `Cb·Ck/Cf`. `WeCr = We_r/12 = y_eq` ✓. A drop from rest (y0=ẏ0=0) has amplitude `a=WeCr`, first peak `y_max = WeCr+a = 2·WeCr`; breakup gate `a+WeCr>1` ⇒ `We_r > bp(3)/2 = 6` = canonical TAB rest-drop critical Weber. User-facing `WeCrit` IS that rest threshold; ×2 maps it to the equilibrium factor. Self-consistent — no discrepancy.

### Verification test T-TAB (rank #1)
Non-tautological triangulation (YupdateTAB is literally the closed form, so code-vs-formula alone only checks transcription):
- **Inputs**: fixed slip `u` (heavy drop, drag→0 so `u`≈const), constant `r`, `μl`,`ρl`,`ρg`,`σ`. Choose `We_r` in [8,30] so `2·We_r/12 > 1` ⇒ oscillation reaches `y=1`. Seed `y0`,`ẏ0`.
- **Oracle A**: derive `ω, td, We_eff=We_r/12` *independently* from {Ck=8,Cd=5,Cf=1/3,Cb=1/2}; evaluate analytic `y(t)`.
- **Oracle B**: RK4-integrate the raw ODE `ÿ = (Cf/Cb)(ρg/ρl)(u²/r²) − (Ck σ/ρl r³)y − (Cd μl/ρl r²)ẏ` with small `h`. Confirms A *solves* the ODE (guards against a mistranscribed constant that would agree between code and A but not B).
- **Compare**: IGLOO YupdateTAB stepping vs A and B; and predicted breakup time `tb` (analytic first crossing `y=1`) vs code's `tb` at [L398](../../src/lib/Lib_Breakup.f90#L398).
- **Tolerance**: A-vs-code ≤ 1e-12 (arithmetic). B-vs-A ≤ O(h⁴)·(sim length) RK4 truncation. Constant-`u` approximation error bounded by `Δu/u` over the window (make drag negligible: large `ρl d²`).
- **Rest-drop threshold sub-test**: `y0=ẏ0=0`, sweep `We_r` across 6 ⇒ breakup fires iff `We_r>6` (checks the `bp(3)` design above).

---

## 2. ETAB — model (5)

[ORA87 §ETAB extension / Tanner 1997]. Same oscillator (`bp(5)=Comega=8`, `bp(6)=Cmu=5`, `bp(3)=WeCrit×2=12`), but product size from exponential mass-rate law instead of energy sampling.
[L472](../../src/lib/Lib_Breakup.f90#L472): `AWe = (k2/k1·√WeTrans − 1)/WeTrans⁴`, defaults `k1=k2=0.2`, `WeTrans=100`.
[L509-512](../../src/lib/Lib_Breakup.f90#L509): `Kbr = k1·ω·(AWe·We⁴+1)` for `We≤WeTrans`, `Kbr = k1·ω·√We` for `We>WeTrans`; `r_new = r·exp(−(Kbr/ω)·acos(1−1/WeCr))`.

### Verification test T-ETAB (rank #5)
**Continuity of `Kbr` at `We=WeTrans`.** Below: `Kbr = k1·ω·(AWe·WeTrans⁴+1) = k1·ω·(k2/k1·√WeTrans) = k2·ω·√WeTrans`. Above: `k1·ω·√WeTrans`. Continuous **iff k1=k2** (true at default). Test: sweep `We` through `WeTrans`, assert `Kbr` continuous to machine eps for k1=k2, and assert the designed jump `k2/k1` for k1≠k2. Pure unit check; oscillator half re-uses T-TAB. ⚠ minor: AWe branch is only continuous at default constants — document, don't fix.

---

## 3. Pilch-Erdman — model (1)

### Paper
[PE87]. Diameter Weber `We = ρg·u²·d0/σ`. Dimensionless total breakup time `T` (nondim by `t* = (d0/u)·√(ρl/ρg)`), piecewise:

| We range | T |
|---|---|
| 12 < We ≤ 18 | 6·(We−12)^(−0.25) |
| 18 < We ≤ 45 | 2.45·(We−12)^(0.25) |
| 45 < We ≤ 351 | 14.1·(We−12)^(0.25) |
| 351 < We ≤ 2670 | 0.766·(We−12)^(0.25) |
| We > 2670 | 5.5 |

Critical Weber (Oh-corrected): `Wc = 12·(1 + 1.077·Oh^1.6)`, `Oh = μl/√(ρl d σ)`.
Max stable fragment: `d_stable = Wc·σ / (ρg·((1−V)·u)²)` where `V` = accumulated dimensionless drop velocity at breakup.

### Code cross-check — T(We) & Wc(Oh) MATCH; d_stable structure MATCHES, V-coeffs UNVERIFIED
[PilchErdman L118-139](../../src/lib/Lib_Breakup.f90#L118): `We = dp·rhog·vel²/sigma` (diameter ✓, unlike the other models which use radius). `Wec = 12(1+1.077·Oh^1.6)` ✓ [L120]. Branch table [L125-137] = paper table **exactly** (5.5 / 0.766 / 14.1 / 2.45 / 6, with (We−12)^±0.25, breakpoints 2670/351/45/18/12) ✓.
`tau = TBT·dp/(vel·√(ρg/ρl))` [L145] = `T·t*` ✓.
Oh-override [L139]: `Oh>0.1 & We<228 ⇒ T = 4.5(1+1.2·Oh^1.64)` — this is the PE87 viscous-drop breakup-time correction; ⚠ coefficients (4.5, 1.2, 1.64) plausible but NOT re-read from PE87 PDF — mark **inferred**.
`VdV = √(ρg/ρl)·(0.75·Cd·T + 3·B·T²)` [L140], `dStable = Wec·σ/(ρg v²(1−VdV)²+toll)` [L141]. The `d_stable` *structure* matches PE87. The velocity polynomial coefficients (0.75, 3) with `Cd=bp(1)=1`, `B=bp(2)=0.116` are a truncated 2-term velocity-history fit; ⚠ **NOT verified against PE87** (PE87 gives a higher-order poly in normalized time). Treat as inferred.

### Verification tests
- **T-PE-time (rank #2)**: tabulate `T(We)` at We ∈ {13,20,50,400,3000} and just inside each break (17.9/18.1, 44.9/45.1, …); assert continuity is *not* required (paper table is discontinuous at 18/45 by design — 2.45 vs 14.1) but each branch value matches the formula to 1e-12. Confirms transcription + branch bounds (`>` vs `≥`). Oracle = the table above.
- **T-PE-Wec (rank #3)**: `Wc(Oh)` at Oh ∈ {0, 0.1, 1}: expect {12, 12·(1+1.077·0.1^1.6), 12·2.077}. 1e-12.
- **T-PE-dStable (rank #6, low)**: only after PE87 velocity-poly confirmed. Blocked on the ⚠ above.

---

## 4. Reitz-Diwakar — model (2)

### Paper / reference form
[RD86], reference implementation = OpenFOAM v7 `ReitzDiwakar.C`. OpenFOAM `We = 0.5·ρg·u²·d/σ` = **radius Weber** (= IGLOO convention). Two regimes:
- **Bag**: `We > Cbag` (Cbag=6). `d_stable = 2·Cbag·σ/(ρg u²)`. `τ_bag = Cb·d·√(ρl d/σ)`, OpenFOAM Cb=0.785.
- **Stripping**: `We > Cstrip·√Re` (Cstrip=0.5). `d_stable = (2 Cstrip σ)²/(ρg u³ μg)`. `τ_strip = Cs·d·√(ρl/ρg)/u`, OpenFOAM Cs=10.

### Code cross-check — MATCHES (constant repackaging)
[ReitzDiwakar L158-168](../../src/lib/Lib_Breakup.f90#L158): `We = radius·ρg·vel²/σ` = OpenFOAM's `0.5·ρg u² d/σ` ✓.
Bag: `We>bp(1)=6` ✓; `dStable = 2σ·bp(1)/(ρg vel²) = 12σ/(ρg u²)` ✓; `τ = bp(2)·dp·√(0.0625·ρl dp/σ) = (π/4)·dp·√(ρl dp/σ) = 0.785·dp√(ρl dp/σ)` (bp(2)=π, √0.0625=1/4 ⇒ π/4=0.785=OpenFOAM Cb) ✓.
Stripping: `We>bp(3)√Re`, bp(3)=0.5 ✓; `τ = bp(4)·0.5·dp·√(ρl/ρg)/u = 10·dp√(ρl/ρg)/u` (bp(4)=Cs=20, ×0.5 ⇒ 10=OpenFOAM Cs) ✓; `dStable=(2·bp(3)σ)²·Re/(dp ρg² u⁴)` = OpenFOAM `(2 Cstrip σ)²/(ρg u³ μg)` **iff `Re = ρg u dp/μg`** (gas Reynolds). ⚠ verify the `Re` passed into `breakupOde` is gas-based diameter Reynolds; if it is particle/other, `dStable` scales wrong. **All RD constants match once repackaged; ⚠ single open item = which Re.**

### Verification test T-RD (rank #4)
Unit checks: (i) bag regime `We_r∈(6, 0.5√Re)` ⇒ `dStable=12σ/(ρg u²)`, `τ` per formula; (ii) stripping `We_r>0.5√Re` ⇒ stripping branch; (iii) sub-critical `We_r<6` ⇒ `npdotDot=0`. Oracle = OpenFOAM formulas above. Confirm `Re` convention first (blocks dStable-strip check).

---

## 5. Reitz-KHRT — model (3) — ⚠ CONTAINS A DIMENSIONAL BUG

### Paper
[RD86 KH] Reitz 1987 KH linear stability: max growth rate and wavelength
`Ω_KH·√(ρl a³/σ) = (0.34 + 0.38 We_g^1.5)/((1+Oh)(1+1.4 T^0.6))` ⇒ `Ω_KH = […]·√(σ/(ρl a³))` [units 1/s];
`Λ_KH/a = 9.02(1+0.45 Oh^0.5)(1+0.4 T^0.7)/(1+0.87 We_g^1.67)^0.6`.
RT: `Ω_RT = √[2·(g_t(ρl−ρg))^1.5 / (3√(3σ)·(ρl+ρg))]`, `Λ_RT = 2π C_RT/√(g_t(ρl−ρg)/(3σ))`.

### Code cross-check
`lambdaKH` [L194](../../src/lib/Lib_Breakup.f90#L194): matches Reitz KH wavelength **exactly** ✓.
`omegaRT` [L201](../../src/lib/Lib_Breakup.f90#L201), `lambdaRT` [L202](../../src/lib/Lib_Breakup.f90#L202), `tauRT=bp(3)/omegaRT` (Ctau=1): match ✓.

> ⚠ **CODE-VS-PAPER DISCREPANCY (correctness-breaking, KHRT only).**
> [L192-193](../../src/lib/Lib_Breakup.f90#L192) (and duplicated [L280-281](../../src/lib/Lib_Breakup.f90#L280)):
> `omegaKH = (0.34+0.38·WeGas^1.5)/((1+Oh)(1+1.4·Tay^0.6)) * (sigma/rhop*radius**3)`
> The trailing factor parses as `(σ/ρl)·r³` → units **m⁶/s²**, NOT a frequency. The Reitz growth rate needs `sqrt(σ/(ρl·r³))` → units 1/s.
> So `omegaKH` is off by both a missing `sqrt` AND `r³` inverted (numerator vs denominator). Consequently `tauKH = 3.726·B1·r/(Λ_KH·omegaKH)` [L196] is dimensionally garbage, and `dStable=2·B0·Λ_KH` [L197] (which is dimensionally fine, ✓) is the only sound KH quantity. **Do NOT fix — flag only.** This is the direct analogue of the evaporation factor-2 the task warned about; it must be corrected before KHRT is trusted, in a separate change.

### Verification test T-KHRT (rank #7, BLOCKED)
A literature oracle for `Ω_KH`, `Λ_KH` exists (dimensionless plateau values at high We: `Ω_KH√(ρl a³/σ)→0.9`, etc.), but the coded `omegaKH` cannot pass any dimensional oracle as written. **Blocked pending the ⚠ fix.** `Λ_KH`, `Ω_RT`, `Λ_RT`, `Λ_RT<d` RT-trigger and child-drop cubic [L308-322] are independently testable and correct; defer to a KHRT-specific pass. Note KH-mass-shed branch [L337-346] is commented out.

---

## Ranked shortlist — implement first

1. **T-TAB (TAB closed-form oscillator)** — exact analytic oracle, all constants (Ck,Cd,Cf,Cb,K) verified vs [ORA87]/ANSYS Eq.24, triangulated code↔analytic↔RK4. Cleanest, strongest. **START HERE.**
2. **T-PE-time** — PE87 table verified exactly; pure branch/transcription unit check, trivial oracle.
3. **T-PE-Wec** — `Wc(Oh)` closed form, 3-point check.
4. **T-RD** — matches OpenFOAM after repackaging; blocks only on Re convention.
5. **T-ETAB** — Kbr continuity at WeTrans (exact at default k1=k2).
6. **T-PE-dStable** — LOW: depends on unverified PE87 velocity poly (0.75,3).
7. **T-KHRT** — BLOCKED by omegaKH ⚠ dimensional bug.

Confirmed: TAB closed-form is the right #1 (task's expected winner).

---

## Open questions for the user
1. **KHRT omegaKH** ⚠ (L192/L280): confirm you want this logged as CODE-VS-PAPER only (no fix), and whether KHRT is in active use — if so it warrants a follow-up fix task. It should be `sqrt(sigma/(rhop*radius**3))`.
2. **RD `Re`**: which Reynolds is passed into `breakupOde` (gas-based `ρg u d/μg`?)? Determines whether the stripping `dStable` matches [RD86]/OpenFOAM.
3. **PE87 velocity poly**: do you have the PE87 PDF? Need to confirm the `VdV` coefficients (0.75·Cd·T + 3·B·T²) and the viscous-time override (4.5,1.2,1.64) before T-PE-dStable is trustworthy.
4. Should these tests live under `tests/breakup/{tab,pilch-erdman,reitz-diwakar}/` (where they now live) with the harness in `test.sh`?
